##### required libraries
library(ggplot2)
library(dplyr)
library(gridExtra)
library(numbers)
library(parallel)
library(ggtern)
library(reshape2)
library(patchwork)
library(combinat)
library(viridis)
library(doParallel)

list.of.packages <- c("ggplot2", "dplyr","gridExtra","numbers","parallel","ggtern",
                      "reshape2","patchwork","combinat","viridis","doParallel")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

########################################################################################
################################Called functions########################################
########################################################################################
### This functions creates the intial frequencies of four genotypes. 
### This function requires the frequencies range (limits) and the step
### CONSIDER THAT THE FUCTION USES LIMITS AS arguments "from" and "to for seq() and step as argument"by".
create.gamete.frequencies<-function(x1.limits,x2.limits,x3.limits,x4.limits,step){
  a_=as.double(seq(x1.limits[1],x1.limits[2],step))
  b_=as.double(seq(x2.limits[1],x2.limits[2],step))
  c_=as.double(seq(x3.limits[1],x3.limits[2],step))
  d_=as.double(seq(x4.limits[1],x4.limits[2],step))
  ab<-data.frame(a=sort(rep(a_,length(b_))),b=b_)
  ### trimming excess of rows that sum more than 1 for x1 and x2
  ab<-ab[rowSums(ab)<=1,]
  rownames(ab)<-NULL
  ### trimming excess of rows that sum more than 1 for x3 and x4
  cd<-data.frame(c=sort(rep(c_,length(d_))),d=d_)
  cd<-cd[rowSums(cd)<=1,]
  rownames(cd)<-NULL
  sum.1.freqs<-data.frame()
  for(row.ab in 1:nrow(ab)){
    ### Adding columsn x3 and x4 to x1 and x2.
    temp.df<-data.frame(x1_init=ab[row.ab,"a"],x2_init=ab[row.ab,"b"],
                        x3_init=cd[,"c"],x4_init=cd[,"d"])
    ### Only rows that sum exactly 1.
    df<-temp.df[rowSums(temp.df)==1,]
    sum.1.freqs<-rbind(sum.1.freqs, df)
  }
  sum.1.freqs$x1_init<-as.double(sum.1.freqs$x1_init)
  sum.1.freqs$x2_init<-as.double(sum.1.freqs$x2_init)
  sum.1.freqs$x3_init<-as.double(sum.1.freqs$x3_init)
  sum.1.freqs$x4_init<-as.double(sum.1.freqs$x4_init)
  
  arrange(sum.1.freqs,x1_init,x2_init,x3_init)
  return(sum.1.freqs)
}

### This functions creates the transmission values of t.
### This function requires the frequencies range (limits) and the step
### CONSIDER THAT THE FUCTION USES LIMITS AS arguments "from" and "to for seq() and step as argument"by".
create.transmission.values<-function(tx.limits,ty.limits,step){
  t.values<-data.frame()
  for(tx in seq(tx.limits[1],tx.limits[2],step)){
    t.values<-rbind(t.values,data.frame(tx=tx,ty=seq(ty.limits[1],ty.limits[2],step)))
  }
  return(t.values)
}

### This functions merges the distinct initial frequencies and the transmission values 
merging.frequencies.and.t.values<-function(gamete.frequencies, transmision.values){
  list.frequency.values<-apply(gamete.frequencies,1,function(X){
    data.frame(x1_init=unname(unlist(X["x1_init"])),
               x2_init=unname(unlist(X["x2_init"])),
               x3_init=unname(unlist(X["x3_init"])),
               x4_init=unname(unlist(X["x4_init"])),
               tx=unname(unlist(transmision.values["tx"])),
               ty=unname(unlist(transmision.values["ty"])))
  })
  do.call("rbind",list.frequency.values)
}




##### Function to simulate two identical drivers
### This functions simulates the evolution of four the four possible genotypes computing.
### The initial genotype frequencies, "frequencies", e.g.=starting.genotype.frequencies
### The recombination frequency, r, "recombi.freq". e.g.=0.5
### The number of sexual generations, "generations", e.g.=10000.
### The driver strength, t, which is general for the two drivers.
### (Logic), Break simulations in case of being steady states lower than 10e-15 change between generation, "break.stable", e.g.=TRUE
### The decimal value at which frequencies are rounded to 1 or 0, "rounding.limit", e.g.=13
Simulate.TwoIdenticalDrivers<-function(frequencies,
                                       recombi.freq,
                                       t, generations,
                                       break.stable){
  #############################################################Function#######################################################
  t=as.double(t); 
  R=as.double(recombi.freq)
  listGenerations<-list()
  listFates<-list()
  AllRecFreq<-data.frame()
  for(r in R){
    x1=as.double(as.vector(unname((frequencies)["x1_init"])));
    x2=as.double(as.vector(unname((frequencies)["x2_init"])));
    x3=as.double(as.vector(unname((frequencies)["x3_init"])));
    x4=as.double(as.vector(unname((frequencies)["x4_init"])))
    W_=NULL
    for(i in 1:(generations)){
      #print(i)
      D = as.double(r*(x1[i]*x4[i]- x2[i]*x3[i]))
      W = as.double(1+ x4[i]*(-t*(2*x1[i]+x2[i]+x3[i])+x1[i]*t^2) + D*(2*t-t^2))
      W_= c(W_,W)
      if(break.stable==TRUE){
        if( near(x1[i],1,tol=1/(10^15)) | 
            near(x2[i],1,tol=1/(10^15)) |
            near(x3[i],1,tol=1/(10^15)) |
            near(x4[i],1,tol=1/(10^15)))
        {break}
        
        
      }
      # if(t==0){
      #   break
      # }else{
        x1=c(x1,as.double( (x1[i]-D)/W_[i] ))
        x2=c(x2,as.double( (x2[i]+D)/W_[i] ))
        x3=c(x3,as.double( (x3[i]+D)/W_[i] ))
        x4=c(x4,as.double( (x4[i]*(1-t*(2*x1[i]+x2[i]+x3[i])+x1[i]*t^2)-D*(1-t)^2)  /W_[i] ))
         if(i>1 & break.stable==TRUE){
        #   #### I had problems with the way R calculates small numbers,the simulation gave me rare numbers when I used a tolerance lower
        #   #### than -15. Calculations were not consistent with hand-made calculations.
           if(near(x1[i-1],x1[i],tol = 10e-15 ) & 
              near(x2[i-1],x2[i],tol = 10e-15 ) &
              near(x3[i-1],x3[i],tol = 10e-15 ) &
              near(x4[i-1],x4[i],tol = 10e-15 ) ){
             c=i+1
             W = as.double(1+ x4[c]*(-t*(2*x1[c]+x2[c]+x3[c])+x1[c]*t^2) + D*(2*t-t^2))
             W_= c (W_,W)  
             break
           }
           
         }
        if(i==(generations) & (i!=1)){
          c=i+1
          W = as.double(1+ x4[c]*(-t*(2*x1[c]+x2[c]+x3[c])+x1[c]*t^2)  + D*(2*t-t^2))
          W_= c (W_,W)
        }
        
    #}
    }
    GenerationsDataDriverDriver<-data.frame(F_xPlus_F_XPlus=x1,
                                            F_xPlus_F_XMinus=x2,
                                            F_xMinus_F_XPlus=x3,
                                            F_xMinus_F_XMinus=x4,
                                            Generations=c(0:(length(x4)-1)) ,
                                            x1_init=rep(unname(unlist(frequencies[1])),
                                                        length(1:length(x1))),
                                            x2_init=rep(unname(unlist(frequencies[2])),
                                                        length(1:length(x2))),
                                            x3_init=rep(unname(unlist(frequencies[3])),
                                                        length(1:length(x3))),
                                            x4_init=rep(unname(unlist(frequencies[4])),
                                                        length(1:length(x4))),
                                            recombFreq=rep(r,
                                                           length(1:length(x1))),
                                            tx=rep(t,
                                                   length(1:length(x1))),
                                            AbsFitness=W_,
                                            id=paste(c(frequencies,r, t),collapse="_"))
    AllRecFreq<-rbind(GenerationsDataDriverDriver,AllRecFreq)
  }
  return(AllRecFreq)
}




### This function simulates the evolution of four genotypes using parallel computing.
### The function requires the package parallel.
### The arguments to run in parallel are:
### Associated cluster, "X", e.g.=1:cores
### The cluster name "cl", e.g.=cl
### The recombination frequency, r, "recombi.freq". e.g.=0.5
### The numbers of cores, "cores", cores
### The initial genotype frequencies, "frequencies", e.g.=starting.genotype.frequencies
### The number of sexual generations, "generations", e.g.=10000.
### (Logic), Break simulations in case of being steady states lower than 10e-15 change between generation, "break.stable", e.g.=TRUE
### The decimal value at which frequencies are rounded to 1 or 0, "rounding.limit", e.g.=13

Simulate.Evolved.PickedPoints.parallel<-function(X,
                                                 cl,
                                                 recombi.freq, cores,
                                                 frequencies,
                                                 generations,
                                                 break.stable,
                                                 rounding.limit){
  ######## Function that calculates frequencies withot parallelizing.
  ######## Parametes are parsed.
  Simulate.DriverXvsDriverY.transmission.advantage.Simp<-function(frequencies,
                                                                  recombi.freq,
                                                                  tx,ty,
                                                                  generations,
                                                                  break.stable){
    #############################################################Function#######################################################
    t_x=as.double(tx); 
    t_y=as.double(ty); 
    R=recombi.freq
    listGenerations<-list()
    listFates<-list()
    for(r in R){
      x1=as.double(as.vector(unname((frequencies)["x1_init"])));x2=as.double(as.vector(unname((frequencies)["x2_init"])));
      x3=as.double(as.vector(unname((frequencies)["x3_init"])));x4=as.double(as.vector(unname((frequencies)["x4_init"])))
      W_=NULL
      for(i in 1:(generations)){
        #print(i)
        D=as.numeric(r*(x1[i]*x4[i]- x2[i]*x3[i]))
        # W = as.double(x1[i]^2+ x2[i]^2+ x3[i]^2+ x4[i]^2+
        #                 (x1[i]*x3[i]+x2[i]*x4[i])*(2-t_x)+ (x1[i]*x2[i]+x3[i]*x4[i])*(2-t_y)+
        #                 x1[i]*x4[i]*(t_x*t_y-t_x-t_y+2) + x2[i]*x3[i]*(-t_x-t_y+2) +
        #                 2*D*(-t_x*t_y)/2)
        W = as.double(1-(x2[i]+x4[i])*(x1[i]+x3[i])*t_y-(x3[i]+x4[i])*(x1[i]+x2[i])*t_x+
                        ((1-r)*x1[i]*x4[i]+r*(x2[i]*x3[i]))*(t_x*t_y))
        W_=c(W_, W)
        if(break.stable==TRUE){
          if( near(x1[i],1,tol=1/(10^15)) | near(x2[i],1,tol=1/(10^15)) |
              near(x3[i],1,tol=1/(10^15)) | near(x4[i],1,tol=1/(10^15)))
          {break}
        }
        
        if(t_x==0 & t_y==0){
          break
        }else{
          # x1=c(x1,as.double( (x1[i]^2+  x1[i]*x2[i]+ x1[i]*x3[i]+ (x1[i]*x4[i]-D))                                        /W ))
          # x2=c(x2,as.double( (x2[i]^2+  x1[i]*x2[i]*(1-t_y) + x2[i]*x4[i] + (x2[i]*x3[i]+D)*(1-t_y)) /W ))
          # x3=c(x3,as.double( (x3[i]^2+  x1[i]*x3[i]*(1-t_x)+  x3[i]*x4[i] + (x2[i]*x3[i]+D)*(1-t_x)) /W ))
          # x4=c(x4,as.double( (x4[i]^2+  x2[i]*x4[i]*(1-t_x)+  x3[i]*x4[i]*(1-t_y) + (x1[i]*x4[i]-D)*(1-t_x)*(1-t_y) )     /W ))

           x1=c(x1,as.double( (x1[i]-D)/W ))
           x2=c(x2,as.double( ((x2[i]*(1- t_y*(x1[i]+x3[i])))+D*(1-t_y)) /W ))
           x3=c(x3,as.double( ((x3[i]*(1- t_x*(x1[i]+x2[i])))+D*(1-t_x)) /W ))
           x4=c(x4,as.double( ((x4[i]*(1- t_x*(x1[i]+x2[i]) - t_y*(x1[i]+x3[i]) + x1[i]*t_x*t_y)- D*(1-t_x)*(1-t_y) ))     /W ))
          
          if(i!=1){
            #### I had problems with the way R calculates small numbers,the simulation gave me rare numbers when I used a tolerance lower
            #### than -13. Calculations were not consistent with hand-made calculations.
            if(near(x1[i-1],x1[i],tol = 10e-15 ) & 
               near(x2[i-1],x2[i],tol = 10e-15 ) &
               near(x3[i-1],x3[i],tol = 10e-15 ) &
               near(x4[i-1],x4[i],tol = 10e-15 ) ){
              c=i+1
              # W = as.double(x1[c]^2+ x2[c]^2+ x3[c]^2+ x4[c]^2+
              #                 (x1[c]*x3[c]+x2[c]*x4[c])*(2-t_x)+ (x1[c]*x2[c]+x3[c]*x4[c])*(2-t_y)+
              #                 x1[c]*x4[c]*(t_x*t_y-t_x-t_y+2) + x2[c]*x3[c]*(-t_x-t_y+2) +
              #                 2*D*(-t_x*t_y)/2)
              
              W = as.double(1-(x2[c]+x4[c])*(x1[c]+x3[c])*t_y-(x3[c]+x4[c])*(x1[c]+x2[c])*t_x+
                          ((1-r)*x1[c]*x4[c]+r*(x2[c]*x3[c]))*(t_x*t_y))
            
              W_= c (W_,W)  
              break}
          }
          
          if(i==(generations) & (i!=1)){
            c=i+1
            W =  as.double(1-(x2[c]+x4[c])*(x1[c]+x3[c])*t_y-(x3[c]+x4[c])*(x1[c]+x2[c])*t_x+
                        ((1-r)*x1[c]*x4[c]+r*(x2[c]*x3[c]))*(t_x*t_y))
            W_= c (W_,W)
          }
          
        }
      }
      GenerationsDataDriverDriver<-data.frame(F_xPlus_F_yPlus=x1,
                                              F_xPlus_F_yMinus=x2,
                                              F_xMinus_F_yPlus=x3,
                                              F_xMinus_F_yMinus=x4,
                                              Generations=c(1:length(x4)) ,
                                              x1_init=rep(unname(unlist(frequencies[1])),
                                                          length(1:length(x1))),
                                              x2_init=rep(unname(unlist(frequencies[2])),
                                                          length(1:length(x2))),
                                              x3_init=rep(unname(unlist(frequencies[3])),
                                                          length(1:length(x3))),
                                              x4_init=rep(unname(unlist(frequencies[4])),
                                                          length(1:length(x4))),
                                              recombFreq=rep(r,
                                                             length(1:length(x1))),
                                              tx=rep(t_x,
                                                     length(1:length(x1))),
                                              ty=rep(t_y,
                                                     length(1:length(x1))),
                                              AbsFitness=W_,
                                              id=paste(c(frequencies,r, t_x,t_y),collapse="_")
      )
      listGenerations<-c(listGenerations,list(GenerationsDataDriverDriver))
      listFates<-c(listFates,list(tail(GenerationsDataDriverDriver,1)))
    }
    do.call("rbind",listGenerations)->GenerationsDataDriverDriver.called
    do.call("rbind",listFates)->listFates.called
    list(GenerationsDataDriverDriver.called,listFates.called)
  }
  #################################################################
  ####### this section allocates each set of frequencies to each core
  #################################################################
  X->Position
  sectionsStart.pre<-base::seq(from=1,
                               to=nrow(frequencies)+round(nrow(frequencies)/cores),
                               by=round(nrow(frequencies)/cores))
  ####Reducing each section to -1 meaning previous position, and removing position 0. Shift limit
  sectionsEnd.pre<-(sectionsStart.pre-1)[-1]
  #### Removing the excess of position to start from
  sectionsStart<-sectionsStart.pre[-length(sectionsStart.pre)]
  #### Assigning the highest value for the last nrow
  sectionsEnd.pre[length(sectionsEnd.pre)]<-nrow(frequencies)
  ### Standardizing.
  sectionsEnd<-sectionsEnd.pre
  Start<-sectionsStart[Position]
  End<-sectionsEnd[Position]
  print(c(Start,End))
  flush.console()
  SelectedInterval<- frequencies[Start:End,]
  #################################################################
  #################################################################
  
  
  #################################################################
  #################################################################
  ##### Parsing each value from the selected interval of sequences to the specific core. 
  current.simulation.list<-apply(SelectedInterval,1,function(X){
    row.X<-unlist(X)
    print(X)
    #### Calling main function to each core
    simulation.values<-
      Simulate.DriverXvsDriverY.transmission.advantage.Simp(
        frequencies = row.X[c("x1_init","x2_init","x3_init","x4_init")],
        recombi.freq = recombi.freq,
        tx = row.X["tx"], ty =row.X["ty"],
        generations=generations,
        break.stable=break.stable)
    simulation.values<-simulation.values[[1]]
    ##### Determining genotype fate stablished by frequency and "rounding.limit"
    rounded.new.values<-tail(simulation.values[,c("F_xPlus_F_yPlus","F_xPlus_F_yMinus","F_xMinus_F_yPlus","F_xMinus_F_yMinus")],1)
    rounded.new.values[near(rounded.new.values,1,tol=1/10^rounding.limit)]<-1
    rounded.new.values[near(rounded.new.values,0,tol=1/10^rounding.limit)]<-0
    fate.desc<-NULL
    if(near(rounded.new.values[,"F_xPlus_F_yPlus"],1,tol=1/10^rounding.limit)){
      fate.desc<-c(fate.desc,"BothDriversFixed")
    }
    if(near(rounded.new.values[,"F_xPlus_F_yMinus"],1,tol=1/10^rounding.limit)){
      fate.desc<-c(fate.desc,"DriverXFixed")
    }
    if(near(rounded.new.values[,"F_xMinus_F_yPlus"],1,tol=1/10^rounding.limit)){
      fate.desc<-c(fate.desc,"DriverYFixed")
    }
    if(near(rounded.new.values[,"F_xMinus_F_yMinus"],1,tol=1/10^rounding.limit)){
      fate.desc<-c(fate.desc,"BothDriversExtinct")
    }
    if( !(near(rounded.new.values[,"F_xPlus_F_yPlus"]  , 1 ,tol=1/10^rounding.limit) |
          near(rounded.new.values[,"F_xPlus_F_yMinus"]  , 1 ,tol=1/10^rounding.limit) |
          near(rounded.new.values[,"F_xMinus_F_yPlus"]  , 1 ,tol=1/10^rounding.limit) |
          near(rounded.new.values[,"F_xMinus_F_yMinus"] , 1 ,tol=1/10^rounding.limit) )){
      fate.desc<-c(fate.desc,"Polymorphic")
    }
    simulation.values$DriversFate<-fate.desc
    return(simulation.values)
  })
  #################################################################
  #################################################################
  current.simulation.df<-do.call("rbind",current.simulation.list)
  return(current.simulation.df)
}


