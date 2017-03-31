###### CODE FOR GITHUB FOR PRODUCING PAPER RESULTS
##### Lee LE, Burke Ó, Foster C, Doherty D
##### "Dynamic time warping: a novel approach for assessment of inter-rater reliability of categorical time series"

###### Updated: 12-2-16
## Note: Throughout, "level" refers to prefix length used


library(dtw)
library(irr)
library(TraMineR)
library(stringr)
library(lubridate)
library(RColorBrewer)
library(xtable)
library(lsr)
library(plyr)


### SET WORKING DIRECTORY ------------

setwd("~/Desktop/Output")



### CREATE VARIABLES -------
##build list of inputs for all available cases
Participants=data.frame(Number=c('001','002','003','004','005','006',
                                 '002','007','009','010','011','013'),
                        Type=c(rep('Main',6),rep('Pilot',6)),
                        Rater1Name=rep("LindsayLee",12),
                        Rater2Name=rep("AngelWong",12),
                        Taxonomy=c(rep("taxonomy for main data set.csv",6),
                                   "taxonomy pilot 2 - P002.csv",
                                   "taxonomy pilot 3 - P007.csv",
                                   "taxonomy pilot 6 - P009.csv",
                                   "taxonomy pilot 5 - P010.csv",
                                   "taxonomy pilot 4 - P011.csv",
                                   "taxonomy pilot 1 - P013.csv"),
                        Iteration=c(7:12,2,3,6,5,4,1),stringsAsFactors=FALSE)

#levels 5 and 6 same for all participants (max seq length is 5)
#level 4 same for all participants except P002, P004 in Main and P002 and P013 in Pilot (max seq length is <=4 for all others)
#levels 1,2,3 different for all participants (max seq length 3 only for Angel's P003Main)

short=Participants[c(2,4,7,12),] #only participants where level 4 is different than 5

##create list of step patterns used for each participant
steps=c('symmetric2','asymmetric',"rabinerJuangStepPattern(3,'c',smoothed=FALSE)","rabinerJuangStepPattern(5,'d',smoothed=FALSE)")

##nice names for steps for legends and titles
prettystep=c('symmetric2','asymmetric','R-J3','R-J5')

##generate colors for plots for steps
#'RdYlBu' has max 11 colors, need to change palette if # of participants or step patterns is higher than this
colstep=brewer.pal(length(steps),"RdYlBu")

##generate colors for plots for all participants - 'Paired' has max 12 colors, need to change palette if # of participants or step patterns is higher than this
partcols=brewer.pal(length(Participants$Number),"Paired")

#order data were calculated
or=c(12,7,8,11,10,9,1:6)

#pretty names for participants
prettypart=paste("Iteration",Participants$Iteration)


### DTW ALIGNMENT FUNCTIONS ---------

### read2: read annotations and level desired and compute Cohen's kappa and Krippendorff's alpha ------
read2=function(p,type,rater1,rater2,lev){
  setwd(paste("~/Documents/Rhodes/MSc Applied Stats/Dissertation/Annotation/",type,sep=''))
  
  ##read annotations
  note1=read.csv(paste('P',p,'_',rater1,'.csv',sep=''),na.strings=c("",NA),quote="")
  note2=read.csv(paste('P',p,'_',rater2,'.csv',sep=''),na.strings=c("",NA),quote="")
  
  ##remove NA rows (eg for photos not on main day of wear)
  note1=note1[!is.na(note1$annotation),]
  note2=note2[!is.na(note2$annotation),]
  
  ##convert time to POSIXlt type (times are same for each rater, only need one vector of time)
  time=strptime(as.character(note1$image_time),format="%Y-%m-%d %H:%M:%S")
  
  ##convert time to decimal hour
  t.ymd <- ymd_hms(time)
  dectime <- hour(t.ymd) + minute(t.ymd)/60 + second(t.ymd)/3600
  
  ##store annotations as characters
  Note1=as.character(note1$annotation)
  Note2=as.character(note2$annotation)
  
  ##make a copy of full annotation
  fullNote1=Note1
  fullNote2=Note2
  
  ##split Note1 and Note2
  Note1=strsplit(Note1,';')
  Note2=strsplit(Note2,';')
  
  ##make each entry in Note1 and Note2 same length as lev, and paste back together
  for (k in 1:length(Note1)) {
    Note1[[k]]=Note1[[k]][1:lev]
    Note1[[k]]=Note1[[k]][!is.na(Note1[[k]])]
    Note2[[k]]=Note2[[k]][1:lev]
    Note2[[k]]=Note2[[k]][!is.na(Note2[[k]])]
    
    Note1[[k]]=paste(Note1[[k]],collapse=';')
    Note2[[k]]=paste(Note2[[k]],collapse=';')
    
  }
  
  Note1=as.character(Note1)
  Note2=as.character(Note2)
  
  ##give new names to annotations for easy reading of confusion matrix
  allval=union(unique(Note1),unique(Note2))
  
  newnames=as.character(1:length(allval))
  
  ##change number codes to new names (ie integers) for ease of reading
  newnameNote1=sapply(Note1,function(i) return(newnames[allval==i][1]))
  newnameNote2=sapply(Note2,function(i) return(newnames[allval==i][1]))
  
  
  ##compute cohen's kappa
  kappascore=kappa2(cbind(newnameNote1,newnameNote2),weight='unweighted')$value
  
  ##compute Krippendorff's alpha
  ka=kripp.alpha(rbind(newnameNote1,newnameNote2))$value
  
  ##compute overall proportion of agreement
  ag=agree(cbind(newnameNote1,newnameNote2))$value/100
  
  ##compile annotations
  compare=data.frame(Time=time,DecTime=dectime,
                     Rater1=Note1,Rater2=Note2,
                     fullRater1=fullNote1,fullRater2=fullNote2,
                     shortRater1=newnameNote1,shortRater2=newnameNote2,
                     stringsAsFactors=FALSE)
  
  ##output the data, kappa score, Krippendorff's alpha, and overall proportion of agreement
  return(list(Data=compare,Kappa=kappascore,KrippAlpha=ka,JPA=ag))
  
}

### mydtw: compute alignments for all participants and step patterns for given level, and complete new simulations---------
mydtw=function(participants,lev){
  
  ##count number of participants
  num=dim(participants)[1]
  
  ##list of names of Results list
  charname=paste0(participants$Number,participants$Type)
  
  ##list of different numerical results we want for each participant and step pattern
  output=c('Kappa',"KrippAlpha","JPA",steps)
  
  ##Results list of normalized distances for each participant for each step pattern
  res <- sapply(participants[,1],function(x) NULL) #empty list for each participant
  Results=lapply(1:num,function(j) {
    res[[j]]=sapply(output,function(x) NULL,simplify='array')}) #list of empty lists
  names(Results)=charname
  
  ##list of coordinates for warping functions for each participant for each step pattern
  lin <- sapply(participants[,1],function(x) NULL) #empty list for each participant
  warps=lapply(1:num,function(j) {
    lin[[j]]=sapply(steps,function(x) NULL,simplify='array')}) #list of empty lists
  names(warps)=charname
  
  ##list of random simulated annotations for each participant
  Randsim.Decom <- sapply(participants[,1],function(x) NULL) #empty list for each participant
  names(Randsim.Decom)=charname
  
  ##list of weighted simulation annotations for each participant
  Wrongsim.Decom <- sapply(participants[,1],function(x) NULL) #empty list for each participant
  names(Wrongsim.Decom)=charname
  
  ##list of data frames for time of each image
  time <- sapply(participants[,1],function(x) NULL) #empty list for each participant
  names(time)=charname
  
  ##intialize storage of sequence lengths
  N=numeric(dim(participants)[1]) #for annotation with certain prefix length
  M=numeric(dim(participants)[1])
  fN=numeric(dim(participants)[1]) #for full annotation (max length prefix)
  fM=numeric(dim(participants)[1])
  
  #####do the following for each participant
  for (i in 1:num) {
    dat=read2(participants$Number[i],   
             participants$Type[i],
             participants$Rater1Name[i],
             participants$Rater2Name[i],
             lev)
    
    ##store Cohen's Kappa and Krippendorff's Alpha - should calculate using level to match distance measures
    Results[[charname[i]]][['Kappa']]=dat$Kappa
    Results[[charname[i]]][['KrippAlpha']]=dat$KrippAlpha
    Results[[charname[i]]][['JPA']]=dat$JPA
    
    ##store timings for each image
    time[[charname[[i]]]]=data.frame(Time=dat$Data$Time,DecTime=dat$Data$DecTime,stringsAsFactors=FALSE)
    
    
    ##data frame of annotations only
    datnote=data.frame(Rater1=dat$Data$Rater1,Rater2=dat$Data$Rater2,
                       fullRater1=dat$Data$fullRater1,fullRater2=dat$Data$fullRater2,stringsAsFactors=FALSE)
    
    ##### USE TraMineR PACKAGE TO CREATE COST MATRIX (DISTANCE MATRIX) FOR CATEGORICAL DATA
    
    ##read in taxonomy
    tax=read.csv(paste('~/Documents/Rhodes/MSc Applied Stats/Dissertation/Annotation/',
                       participants$Type[i],'/',participants$Taxonomy[i],sep=''),header=FALSE)
    tax=as.character(tax$V1)
    
    ##remove commas at end of taxonomy
    c=sapply(tax,function(x) str_sub(x,start=-1))
    commas=which(c==',')
    for (k in commas) {
      tax[k]=substr(tax[k], 1, nchar(tax[k])-1)
    }
    
    suppressMessages(tax<-seqdecomp(tax,sep=';'))
    suppressMessages(tax.seq<-seqdef(tax))
    taxcol=ncol(tax.seq) #see number of columns Rater1 and Rater2 must have
    
    ##convert annotations to matrix representation of sequences
    suppressMessages(Rater1<-seqdecomp(datnote$Rater1,sep=';'))
    suppressMessages(Rater2<-seqdecomp(datnote$Rater2,sep=';'))
    suppressMessages(fullRater1<-seqdecomp(datnote$fullRater1,sep=';'))
    suppressMessages(fullRater2<-seqdecomp(datnote$fullRater2,sep=';'))
    
    ##make number of columns the same=max possible length of sequence by adding missing-value columns
    Rater1col=ncol(Rater1)
    Rater2col=ncol(Rater2)
    fRater1col=ncol(fullRater1)
    fRater2col=ncol(fullRater2)
    
    #print info about hierarchy
    cat('P',charname[i],participants$Rater1Name[i],'max seq length:',Rater1col,
        '\nP',charname[i],participants$Rater2Name[i],'max seq length:',Rater2col,
        '\nTaxonomy max seq length:',taxcol,
        '\nLevel of variation chosen:',lev,'\n')
    
    N[i]=dim(Rater1)[1] #necessary to store string lengths when plotting multiple warping functions together
    M[i]=dim(Rater2)[1]
    fN[i]=dim(fullRater1)[1]
    fM[i]=dim(fullRater2)[1]
    
    if (Rater1col!=taxcol) {
      while (Rater1col<taxcol) {
        Rater1=cbind(Rater1,rep(NA,N[i]))
        Rater1col=ncol(Rater1)
        dimnames(Rater1)[[2]][Rater1col]=paste0('[',Rater1col,']')
      }
    }
    
    if (Rater2col!=taxcol) {
      while (Rater2col<taxcol) {
        Rater2=cbind(Rater2,rep(NA,M[i]))
        Rater2col=ncol(Rater2)
        dimnames(Rater2)[[2]][Rater2col]=paste0('[',Rater2col,']')
      }
    }
    
    if (fRater1col!=taxcol) {
      while (fRater1col<taxcol) {
        fullRater1=cbind(fullRater1,rep(NA,fN[i]))
        fRater1col=ncol(fullRater1)
        dimnames(fullRater1)[[2]][fRater1col]=paste0('[',fRater1col,']')
      }
    }
    
    if (fRater2col!=taxcol) {
      while (fRater2col<taxcol) {
        fullRater2=cbind(fullRater2,rep(NA,fM[i]))
        fRater2col=ncol(fullRater2)
        dimnames(fullRater2)[[2]][fRater2col]=paste0('[',fRater2col,']')
      }
    }
    
    ##make sure level isn't higher than taxonomy level
    if (lev>taxcol) lev<-taxcol
   
    ##combine Rater1 and Rater2 into one long matrix of sequences for reading into seqdist
    Seq=rbind(Rater1,Rater2)[,1:lev]
    suppressMessages(Seq<-seqdef(as.matrix(Seq))) #only look at sequence up to certain categorical level
   
    ##compute normalized Longest Common Prefix distance
    suppressMessages(cost<-seqdist(Seq,method="LCP",norm=TRUE))
    cost=cost[1:N[i],(N[i]+1):(N[i]+M[i])] #extract only Rater1.seq (query) vs. Rater2.seq (ref) part
   
    randsim.decom=randsim(tax,fullRater2) #single matrix of randomly simulated sequence
    wrongsim.decomlist=wrongsim(tax,fullRater2,charname[i]) #list of matrices of weighted simulation sequences, list length equals number of probablity weights
   
    Randsim.Decom[[i]]=randsim.decom
    Wrongsim.Decom[[i]]=wrongsim.decomlist
   
    ##### INPUT DISTANCE MATRIX INTO DTW
   
    ##output message
    cat('Starting alignments for ',charname[i],'\n')
   
    for (j in steps) {
      ##compute the alignment for each step pattern
      align=dtw(cost,keep=TRUE,step=eval(parse(text=j)))
      
      ##do simulation just for level specified - simulations used in paper are computed using simresults, which uses the same simulated annotations at each level
      randsimdist=randsimalign(randsim.decom,fullRater2,j,lev,charname[i]) #single value
      wrongsimdist=wrongsimalign(wrongsim.decomlist,fullRater2,j,lev,charname[i]) #vector of values for each probability weight
      
      ##store warp function for each step pattern
      warps[[charname[i]]][[j]]=data.frame(index1=align$index1,index2=align$index2) #store coordinates for warping function
      
      ##print pdf of warp function
      pdf(paste('level',lev,"_warpP",charname[i],'_',prettystep[which(steps==j)],'.pdf',sep=''))
      plot(align,type='alignment', xlab="rater 1 (query) index", ylab="rater 2 (reference) index",
           main=paste("Warping function for ",charname[i],' prefix length ',lev,'\nand step ',prettystep[which(steps==j)],sep=''))
      abline(0,1,col=2,lty=2)
      dev.off()
      
      ##print pdf of density plot
      pdf(paste('level',lev,"_densityP",charname[i],'_',prettystep[which(steps==j)],'.pdf',sep='')) #slow, note from documentation: The density plot is more colorful than useful.
      dtwPlotDensity(align, normalize=TRUE,
                     xlab="rater 1 (query) index", ylab="rater 2 (reference) index",
                     main=paste("Contour map & warping function for ",charname[i],' prefix length ',lev,'\nand step ',prettystep[which(steps==j)],sep='')) 
      dev.off()
      
      ##store normalized distance as measures of fit quality and simulation normalized distances
      Results[[charname[i]]][[j]]=list(NormDist=align$normalizedDistance,RandSim=randsimdist,WrongSim=wrongsimdist)

    }
   
   ##output message
   cat('Finished all alignments for ',charname[i],'\n')
   
  }
  
  ##return list of all numerical results, the warping functions' coordinates, the random sequences computed, the weighted sequences computed
  return(list(Results=Results,warps=warps,Randsim.Decom=Randsim.Decom,Wrongsim.Decom=Wrongsim.Decom))
}


### SIMULATION FUNCTIONS ------------

### randsim: simulate random annotation ----------
randsim=function(tax,rater) {
  
  n=nrow(as.matrix(rater))
  
  ##build random simulation from same taxonomy as participant
  randsim.decom=as.matrix(tax)[sample(nrow(as.matrix(tax)),n,replace=TRUE),]
  
  ##return normalized distance
  return(randsim.decom)
}

### randsimalign: alignment for random simulation ---------
randsimalign=function(randsim.decom,rater,st,lev,part) {
  
  n=nrow(as.matrix(rater))
  
  ##bind random simulation and reference sequence together, then cut level
  Seqrandsim=rbind(as.matrix(randsim.decom),as.matrix(rater))[,1:lev]
  suppressMessages(Seqrandsim<-seqdef(as.matrix(Seqrandsim)))
  
  ##compute LCP cost matrix
  suppressMessages(costrandsim<-seqdist(Seqrandsim,method="LCP",norm=TRUE))
  costrandsim=costrandsim[1:n,(n+1):(2*n)]
  
  ##compute alignment
  alignrandsim=dtw(costrandsim,keep=TRUE,step=eval(parse(text=st)))
  
  ##output message
  cat('Finished randomized simulation alignment for level ',lev,' ',part,' ',st,'\n')
  
  return(alignrandsim$normalizedDistance)
}

probweights=seq(0,1,by=0.05) #create vector of probabilities of match at each image

### wrongsim: simulate weighted annotation--------
wrongsim=function(tax,rater,part) {
  
  n=nrow(as.matrix(rater)) #number of annotations in reference sequence
  m=nrow(as.matrix(tax))
  wrongsim.decom=matrix(0,nrow=n,ncol=ncol(as.matrix(rater))) #intialize matrix to store annotations
  
  ##initialize list of decomposed simulated sequence matrices wrongsim.decom for each prob weight
  decomlist <- sapply(as.character(probweights),function(x) NULL) #empty list for each probweight
  
  ##create decomposed simulated sequence matrices
  for (p in probweights) {
    for (a in 1:n){
      index=which(apply(tax,1,function(x) identical(x,as.matrix(rater)[a,]))) #look for taxonomy index with same annotation as rater annotation
      probnote=rep((1-p)/(m-1),m) #make probdist equal probabilities for all annotations in taxonomy...
      probnote[index]=p #...except for the one that is the same as rater
      wrongsim.decom[a,]=as.matrix(tax)[sample(m,1,replace=TRUE,prob=probnote),]
    }
    
    decomlist[[as.character(p)]]=wrongsim.decom
    
    ##output message
    cat('Finished weighted simulation for ',p,'\n')
  }
  
  ##output message
  cat('Finished all weighted simulations for ',part,'\n')
  
  return(decomlist) 
}

### wrongsimalign: compute alignments for weighted simulation for one participant, one step type (all probability weights) -------
wrongsimalign=function(decomlist,rater,st,lev,part) {
  
  n=nrow(as.matrix(rater)) #number of annotations in reference sequence
  dlen=length(decomlist) #same length as probability weights vector
  probdist=numeric(length(probweights)) #intialize vector to store normalize distances
  
  for (d in 1:dlen) {
    ##bind random simulation and reference sequence together, then cut level
    Seqwrongsim=rbind(as.matrix(decomlist[[d]]),as.matrix(rater))[,1:lev]
    suppressMessages(Seqwrongsim<-seqdef(as.matrix(Seqwrongsim)))
    
    ##compute LCP cost matrix
    suppressMessages(costwrongsim<-seqdist(Seqwrongsim,method="LCP",norm=TRUE))
    costwrongsim=costwrongsim[1:n,(n+1):(2*n)]
    
    ##compute alignment
    alignwrongsim=dtw(costwrongsim,keep=TRUE,step=eval(parse(text=st)))
    probdist[d]=alignwrongsim$normalizedDistance
  }
  
  ##output message
  cat('Finished all weighted simulation alignments for level ',lev,' ',part,' ',st,'\n')
  
  return(probdist)
}

### simresults: create simulation results for all levels based on same simulated sequence at each prob weight---------
simresults=function(Randsim.Decom,Wrongsim.Decom,maxlev,participants) {
  
  ##list of names of Results list
  charname=paste0(participants$Number,participants$Type)
  num=length(participants$Number)
  
  ##intialize storage of randsimdists
  rand <- sapply(participants[,1],function(x) NULL) #empty list for each participant
  randsimdists=lapply(1:num,function(j) {
    rand[[j]]=sapply(steps,function(x) numeric(maxlev),simplify=FALSE)}) #list of empty lists
  names(randsimdists)=charname
  
  ##intialize storage of wrongsimdists
  wrong <- sapply(participants[,1],function(x) NULL) #empty list for each participant
  wrongsimdists=lapply(1:num,function(j) {
    wrong[[j]]=sapply(steps,function(x) matrix(0,nrow=maxlev,ncol=length(probweights),dimnames=list(Levels=as.character(1:maxlev),Weights=as.character(probweights))),simplify=FALSE)}) #list of empty lists
  names(wrongsimdists)=charname
  
  for (part in charname) {
    i=which(charname==part)
    dat=read2(participants$Number[i],   
              participants$Type[i],
              participants$Rater1Name[i],
              participants$Rater2Name[i],
              maxlev)
    
    fullRater2=dat$Data$fullRater2
    
    ##read in taxonomy to create alphabet
    tax=read.csv(paste('~/Documents/Rhodes/MSc Applied Stats/Dissertation/Annotation/',
                       participants$Type[i],'/',participants$Taxonomy[i],sep=''),header=FALSE)
    tax=as.character(tax$V1)
    
    ##remove commas at end of taxonomy
    c=sapply(tax,function(x) str_sub(x,start=-1))
    commas=which(c==',')
    for (k in commas) {
      tax[k]=substr(tax[k], 1, nchar(tax[k])-1)
    }
    
    suppressMessages(tax<-seqdecomp(tax,sep=';'))
    suppressMessages(tax.seq<-seqdef(tax))
    taxcol=ncol(tax.seq) #see number of columns Rater1 and Rater2 must have
    
    ##make sure reference sequence has same number of columns as taxonomy
    suppressMessages(fullRater2<-seqdecomp(fullRater2,sep=';'))
    fRater2col=ncol(fullRater2)
    fM=dim(fullRater2)[1]
    
    if (fRater2col!=taxcol) {
      while (fRater2col<taxcol) {
        fullRater2=cbind(fullRater2,rep(NA,fM))
        fRater2col=ncol(fullRater2)
        dimnames(fullRater2)[[2]][fRater2col]=paste0('[',fRater2col,']')
      }
    }
    
    for (st in steps)  {
      for (l in 1:maxlev) {
        randsimdists[[part]][[st]][l]=randsimalign(Randsim.Decom[[part]],fullRater2,st,l,part)
        wrongsimdists[[part]][[st]][l,]=wrongsimalign(Wrongsim.Decom[[part]],fullRater2,st,l,part)
      }
    }
  }
  return(list(Randsimdists=randsimdists,Wrongsimdists=wrongsimdists))
}


### simKresults: create simulation of kappa score results for all levels based on same simulated sequence at each prob weight--------------
simKresults=function(Randsim.Decom,Wrongsim.Decom,maxlev,participants) {
  
  ##list of names of Results list
  charname=paste0(participants$Number,participants$Type)
  num=length(participants$Number)
  
  ##intialize storage of randsimkappa
  rand <- sapply(participants[,1],function(x) NULL) #empty list for each participant
  randsimkappa=lapply(1:num,function(j) {rand[[j]]=numeric(maxlev)}) #length 5 0 vector for each participant
  names(randsimkappa)=charname
  
  ##intialize storage of wrongsimkappa
  wrong <- sapply(Participants[,1],function(x) NULL) #empty list for each participant
  wrongsimkappa=lapply(1:num,function(j) {
    wrong[[j]]=matrix(0,nrow=maxlev,ncol=length(probweights),
                      dimnames=list(Levels=as.character(1:maxlev),Weights=as.character(probweights)))}) #list of empty matrices
  names(wrongsimkappa)=charname
  
  for (part in charname) {
    i=which(charname==part)
    dat=read2(participants$Number[i],   
              participants$Type[i],
              participants$Rater1Name[i],
              participants$Rater2Name[i],
              maxlev)
    
    fullRater2=dat$Data$fullRater2
    
    ##read in taxonomy to create alphabet
    tax=read.csv(paste('~/Documents/Rhodes/MSc Applied Stats/Dissertation/Annotation/',
                       participants$Type[i],'/',participants$Taxonomy[i],sep=''),header=FALSE)
    tax=as.character(tax$V1)
    
    ##remove commas at end of taxonomy
    c=sapply(tax,function(x) str_sub(x,start=-1))
    commas=which(c==',')
    for (k in commas) {
      tax[k]=substr(tax[k], 1, nchar(tax[k])-1)
    }
    
    suppressMessages(tax<-seqdecomp(tax,sep=';'))
    suppressMessages(tax.seq<-seqdef(tax))
    taxcol=ncol(tax.seq) #see number of columns Rater1 and Rater2 must have
    
    ##make sure reference sequence has same number of columns as taxonomy
    suppressMessages(fullRater2<-seqdecomp(fullRater2,sep=';'))
    fRater2col=ncol(fullRater2)
    fM=dim(fullRater2)[1]
    
    if (fRater2col!=taxcol) {
      while (fRater2col<taxcol) {
        fullRater2=cbind(fullRater2,rep(NA,fM))
        fRater2col=ncol(fullRater2)
        dimnames(fullRater2)[[2]][fRater2col]=paste0('[',fRater2col,']')
      }
    }
    
    
    for (l in 1:maxlev) {
      #FIX: APPLY CAN'T WORK FOR LEVEL=1
      
      Rand=apply(as.matrix(Randsim.Decom[[part]][,1:l]),1,paste,collapse=";")
      Orig=apply(as.matrix(fullRater2[,1:l]),1,paste,collapse=";")
      
      randsimkappa[[part]][l]=kappa2(cbind(Rand,Orig),weight="unweighted")$value
      
      for (p in 1:length(probweights)) {
        
        Wrong=apply(as.matrix(Wrongsim.Decom[[part]][[p]][,1:l]),1,paste,collapse=";")
        wrongsimkappa[[part]][l,p]=kappa2(cbind(Wrong,Orig),weight="unweighted")$value
        
      }
    }
    
  }
  return(list(Randsimkappa=randsimkappa,Wrongsimkappa=wrongsimkappa))
}



##### PLOT FUNCTIONS -------------------
### samepartallstep: create plots of same participant with all step patterns----------
samepartallstep=function(warps,lev){
  
  ##store participants used
  charname=names(warps)
  num=length(charname)
  
  ##do for each participant
  for (i in 1:num) {
    setwd(paste0("~/Documents/Rhodes/MSc Applied Stats/Dissertation/Annotation/",str_sub(charname[i],start=4)))
    
    N=round_any(max(warps[[i]][[1]][[1]]),500,ceiling) #xlim and ylim as closest multiple of 500 to maximum length of image time series
    
    ##initialize plot
    pdf(paste('level',lev,'_warps_sameP',charname[i],'.pdf',sep=''))
    mar.default <- c(5,4,4,2) + 0.1
    par(mar=mar.default+c(0,1,0,0))
    plot(1,type='n',xlim=c(0,N),ylim=c(0,N),main=charname[i],cex.lab=1.75,cex.axis=1.75,cex.main=1.75,
         xlab="Rater 1 (query) index",ylab="Rater 2 (reference) index")
    abline(0,1,col=1,lty=2)
    
    ##print warps for each step pattern
    for (j in steps) {
      lines(warps[[charname[i]]][[j]]$index1,warps[[charname[i]]][[j]]$index2,col=colstep[which(steps==j)])
    }
    
    legend("topleft",legend=prettystep,col=colstep,lty=1,cex=1.75)
    par(mar=mar.default)
    dev.off() 
  }
}

### allpartsamestepTIME: create plots of all participants warping functions with same step pattern ALIGNED BY TIME ---------
allpartsamestepTIME=function(warps,participants,lev){

  ##store participants used
  charname=names(warps)
  num=length(charname)
  
  ##generate colors for plots for participants used for the function
  colpart=partcols[match(charname,paste0(Participants$Number,Participants$Type))]
  
  ##find half point for participants for formatting legend
  h=ceiling(length(charname)/2)
  
  ##create empty list of image times for each participant
  time <- sapply(charname,function(x) NULL)
  names(time)=charname
  
  ##do for each step pattern
  for (j in steps) {
    
    ##intialize plot
    pdf(paste('level',lev,'_warps_samestepAlign_',prettystep[which(steps==j)],'.pdf',sep=''))
    mar.default <- c(5,4,4,2) + 0.1
    par(mar=mar.default+c(0,1,0,0))
    plot(1,type='n',xlim=c(0,24),ylim=c(0,24),cex.axis=1.75,cex.lab=1.75,cex.main=1.75, #main=prettystep[which(steps==j)]
        xlab="Rater 1 (query) time (hours)",ylab="Rater 2 (reference) time (hours)")
  
    ##plot warp functions aligned by time for each participants
    for (i in 1:num) {
      ##read in times for each image
      dat=read2(participants$Number[i],   
                participants$Type[i],
                participants$Rater1Name[i],
                participants$Rater2Name[i],
                lev)
      
      ##store timings for each image
      time[[charname[[i]]]]=data.frame(Time=dat$Data$Time,DecTime=dat$Data$DecTime,stringsAsFactors=FALSE)
      
      ##store warping paths
      w1=warps[[charname[i]]][[j]]$index1 #warping index 1
      w2=warps[[charname[i]]][[j]]$index2 #warping index 2
    
      ##plot warping paths transformed by time
      points(time[[charname[[i]]]]$DecTime[w1],time[[charname[[i]]]]$DecTime[w2],col=colpart[i],cex=0.5,pch='.')
    }
    
    ##split legend into two pieces for readability
    L=legend("topleft",legend=prettypart[1:8],col=colpart[1:8],pch='.',cex=1.75,pt.cex=7,bty='n')
    legend(x=L$rect$left+L$rect$w,y=L$rect$top,legend=prettypart[(9):num],col=colpart[(9):num],pch='.',cex=1.75,pt.cex=7,bty='n') #was x="bottomright"
    par(mar=mar.default)
    dev.off()
  }
}

### allpartsamestep2: create plots of all participants warping functions with same step pattern ---------
#index doesn't refer to same time of day across participants
#PUB: COLOR DOES NOT MATCH OTHER PLOTS AND LEGEND FORMATTED FOR SMALL NUMBER OF PARTICIPANTS
allpartsamestep2=function(warps,lev){
  
  ##store participants used
  charname=names(warps)
  num=length(charname)
  
  ##generate colors for plots for participants used for the function
  colpart=partcols[c(2,6)]
  
  ##get nice participant names for plot
  partnames=prettypart[which(paste0(Participants$Number,Participants$Type) %in% names(newlist))]
  
  ##do for each step pattern
  for (j in steps) {
    
    ##intialize plot
    pdf(paste('level',lev,'_warps_samestep_',prettystep[which(steps==j)],'.pdf',sep=''))
    mar.default <- c(5,4,4,2) + 0.1
    par(mar=mar.default+c(0,1,0,0))
    plot(1,type='n',xlim=c(0,3000),ylim=c(0,3000),cex.axis=1.75,cex.main=1.75,cex.lab=1.75, #main=prettystep[which(steps==j)]
         xlab="Rater 1 (query) index",ylab="Rater 2 (reference) index")
    abline(0,1,col=1,lty=2)
    
    ##plot warp functions aligned by index for each participants
    for (i in 1:num) {
      lines(warps[[charname[i]]][[j]]$index1,warps[[charname[i]]][[j]]$index2,col=colpart[i])
    }
    
    ##split legend into two pieces for readability
    legend("topleft",legend=partnames,col=colpart,lty=1,cex=1.75)
    par(mar=mar.default)
    dev.off()
  }
}

###distbar2: barplot of normalized distances----------------
#PUB: COLOR DOES NOT MATCH OTHER PLOTS AND ONLY FOR RJ-5, puts in order by iteration
distbar2= function(results,lev){
  
  ##create empty dataframe for normalized distances
  df <- data.frame(matrix(ncol = length(names(results)), nrow = length(steps)))
  colnames(df) <- names(results)
  rownames(df) <- steps
  
  ##extract normalized distances and put into data frame
  for (i in 1:length(names(results))) {
    d=unlist(lapply(results[[i]][steps],'[[','NormDist')) 
    df[,i]=d
  }
  
  #use only RJ-5
  df=df[which(rownames(df)=="rabinerJuangStepPattern(5,'d',smoothed=FALSE)"),]
  
  #put in order that data was gathered
  df=df[,or]
  
  ##print pdf of barplot
  pdf(paste0('level',lev,'_barplotdist.pdf'))
  mar.default <- c(5,4,4,2) + 0.1
  par(mar=mar.default+c(1,0,0,0),xpd=TRUE)
  mp=barplot(as.matrix(df),beside=TRUE,ylim=c(0,0.1+max(df)), col=partcols[2], xaxt='n', space=c(0,2),
             ylab="Normalized distance",xlab='',main='')
  text(x=colMeans(mp),y=-0.01 , srt = 45, adj = 1,labels = names(results)[or], xpd = TRUE)
  mtext("Participant",1,4)
  par(mar=mar.default,xpd=FALSE)
  dev.off()   
  
  ##print pdf of line graph - ONLY FOR ONE STEP PATTERN
  pdf(paste0('level',lev,'_lineplotdist.pdf'))
  mar.default <- c(5,4,4,2) + 0.1
  par(mar=mar.default+c(0,1,0,0),xpd=TRUE)
  plot(as.matrix(df)[1,],ylim=c(0,0.1+max(df)),main='',ylab="Normalized distance",xlab="Iteration",
       type="l",cex.axis=1.75,cex.main=1.75,cex.lab=1.75,lwd=2)
  par(mar=mar.default,xpd=FALSE)
  dev.off()
  
  
  return(df)
}

#se: compute standard error-----------------
se <- function(x) sqrt(var(x)/length(x))

### simplotCI2: plot mean simulation results outputed by simresults given 1-alpha CI using either z or t score -------------
#PUB: COLOR DOES NOT MATCH OTHER PLOTS, DOES NOT OUTPUT ALL LEVELS AND ONLY FOR RJ-5
simplotCI2=function(resultssim,maxlev,resultssimK,score.type,alpha) {
  
  ##vector of participants used
  charname=names(resultssim[['Randsimdists']])
  
  ##vector of step functions used (OVERRIDES steps in global environment): 
  steps="rabinerJuangStepPattern(5,'d',smoothed=FALSE)"
  
  ##calculate critical value
  if (score.type=="z") {crit <- qnorm(1-alpha/2)}
  else if (score.type=="t") {crit <- qt(1-alpha/2,df=(length(charname)-1))}
  
  for (lev in c(5)) {
    
    ##set up plot
    pdf(paste0('level',lev,'_simplotCI_',score.type,'.pdf'))
    mar.default <- c(5,4,4,2) + 0.1
    par(mar=mar.default+c(0,0.5,2.75,0),xpd=TRUE)
    plot(1,type='n', ylab="Normalized distance or -kappa + 1",xlab='Probability weight',xlim=c(0,1),ylim=c(0,1),main='',cex.lab=1.75,cex.axis=1.75)
    legend("topright",c("DTW","Kappa"),col=c(1,2),pch=1)
    
    randd=numeric(length(charname))
    wrongd=matrix(0,length(charname),length(probweights))
    rownames(wrongd)=charname
    colnames(wrongd)=as.character(probweights)
    
    randdK=numeric(length(charname))
    wrongdK=matrix(0,length(charname),length(probweights))
    rownames(wrongdK)=charname
    colnames(wrongdK)=as.character(probweights)
    
    for (r in 1:length(charname)) {
      
      ##extract normalized distances from simulation and put into data frame
      randd[r]=resultssim[['Randsimdists']][[r]][[steps]][lev] #extract random simulation normalized distance
      wrongd[r,]=resultssim[['Wrongsimdists']][[r]][[steps]][lev,] #extract weighted simulation normalized distances
      
      randdK[r]=resultssimK[['Randsimkappa']][[r]][lev] #extract random simulation normalized distance
      wrongdK[r,]=resultssimK[['Wrongsimkappa']][[r]][lev,] #extract weighted simulation normalized distances
      
    }
    
    #mean and CI for normalized distance
    randmean=mean(randd) #mean of random simulation
    wrongmean=apply(wrongd,2,mean) #vector of means of weighted sim across participants for each prob weight
    se.rand=sqrt(var(randd)/length(randd)) #SE of random simulation
    se.wrong=apply(wrongd,2,se)
    
    ciL.rand=randmean-crit*se.rand #CI for rand sim
    ciU.rand=randmean+crit*se.rand
    
    ciL.wrong=wrongmean-crit*se.wrong #CI for wrong sim
    ciU.wrong=wrongmean+crit*se.wrong
    
    
    #mean and CI for kappa
    randmeanK=mean(randdK) #mean of random simulation
    wrongmeanK=apply(wrongdK,2,mean) #vector of means of weighted sim across participants for each prob weight
    se.randK=sqrt(var(randdK)/length(randdK)) #SE of random simulation
    se.wrongK=apply(wrongdK,2,se)
    
    ciL.randK=randmeanK-crit*se.randK #CI for rand sim
    ciU.randK=randmeanK+crit*se.randK
    
    ciL.wrongK=wrongmeanK-crit*se.wrongK #CI for wrong sim
    ciU.wrongK=wrongmeanK+crit*se.wrongK
    
    #plot CI and mean lines
    polygon(c(probweights,rev(probweights)),c(ciL.wrong,rev(ciU.wrong)),col="gray",border="gray")
    polygon(c(probweights,rev(probweights)),c((-ciL.wrongK+1),rev(-ciU.wrongK+1)),col="pink",border="pink")
    
    lines(x=probweights,y=wrongmean,col=1,type='b') #plot points of mean of weighted sim kappa (flipped vertically +1)
    lines(x=probweights,y=(-wrongmeanK+1),col=2,type='b') #plot points of mean of weighted sim ND
    
    par(mar=mar.default,xpd=FALSE)
    dev.off()
  }
  
  return(rbind(se.wrong,ciU.wrong,wrongmean,ciL.wrong,se.wrongK,ciU.wrongK,wrongmeanK,ciL.wrongK))
  
}

### RUN IT, PRINT IT -----------------------

f=file("results_level5_DATE.txt")
sink(f)
sink(f,type='message')
system.time(results5decom <- mydtw(Participants,5))
results5decom[['Results']]
newlist=results5decom[['warps']][c(6,12)] #extract warping functions for 006Main and 013Pilot
samepartallstep(results5decom[['warps']],5)
allpartsamestep2(results5decom[['warps']],5)
allpartsamestepTIME(results5decom[['warps']],Participants,5)
distbar2(results5decom[['Results']],5)
results5sim=simresults(results5decom[['Randsim.Decom']],results5decom[['Wrongsim.Decom']],5,Participants)
results5simK=simKresults(results5decom[['Randsim.Decom']],results5decom[['Wrongsim.Decom']],5,Participants)
ci.plot.new.t=simplotCI2(results5sim,5,results5simK,alpha = 0.05,score.type = "t") #output sim plot with mean and 95% CI using t-score
sink()
sink(type='message')


### Save workspace
save.image("~/Desktop/Output/resultsDATE.RData")

