##################################################
# Information content of stepped wedge and other 
# centrosymmetric longitudinal cluster randomised 
# trial designs with unequal cluster-period sizes
# J Kasza, AB Forbes
#
# Functions to aid in replicating figures in paper
# J Kasza, 2019-10-01
# jessica.kasza@monash.edu
##################################################

#Generate a standard Stepped Wedge design matrix
SWdesmat <- function(T) {
  Xsw <- matrix(data=0, ncol = T, nrow = (T-1))
  for(i in 1:(T-1)) {
    Xsw[i,(i+1):T] <- 1
  }
  return(Xsw)
}

#Take all of the permutations of a set of clusters and sort within
#blocks of mym (i.e. for the 8-cluster, 4-sequence example, although 
#there are 40320 ways to arrange the 8 clusters, the arrangement within
#sequences does not matter (i.e. 1,2 assigned to seq 1 is the same as
#2,1 assigned to seq 1)). sortbym returns only the distinct SW designs.
sortbym <- function(myperms, mym){
  #The number of sequences in the SW:
  numberseqs <- ncol(myperms)/mym
  #want to sort each mym-tuple by size
  sorted <- matrix(data=NA, nrow=nrow(myperms), ncol=ncol(myperms))
  for(i in 1:numberseqs){
    sorted[,((i-1)*mym+1):(i*mym)] <- t(apply(myperms[,((i-1)*mym+1):(i*mym)], MARGIN=1, sort))
  }
  #return(sorted)
  #Delete duplicate rows:
  return(distinct(as.data.frame(sorted)))
  
}


#Calculate the variance of the CRT when clusters have 
#different cluster-period sizes
CRTVar_mMAT <- function(Xmat, Mmat, rho0, r, type) {
  
  totalvar <- 1
  sig2CP <- rho0*totalvar
  sig2E <- totalvar - sig2CP
  sig2 <- sig2E/Mmat #this is now a matrix. 
  
  T <- ncol(Xmat)
  K <- nrow(Xmat)
  
  
  Xvec <- as.vector(t(Xmat))
  stackI <- matrix(rep(diag(1,T)), nrow=K*T, ncol=T, byrow=TRUE)
  Zmat <- cbind(stackI[!is.na(Xvec),], Xvec[!is.na(Xvec)])
  
  #Variance matrix for first cluster, with decay in correlation over time
  #Vi <- diag(sig2,T) + matrix(sig2CP,nrow=T, ncol=T)
  #Constant decay var if type==0
  if(type==0) { 
    Vi <-diag(sig2[1,] +(1-r)*sig2CP, T) + matrix(data=sig2CP*r, nrow=T, ncol=T)
  }
  #exponential decay structure
  if(type==1) { 
    Vi <- diag(sig2[1,],T) + sig2CP*(r^abs(matrix(1:T,nrow=T, ncol=T, byrow=FALSE) - matrix(1:T,nrow=T, ncol=T, byrow=TRUE)))
  }
  
  for(k in 2:K){
    if(type==0) { 
      Vi <- bdiag(Vi, diag(sig2[k,] +(1-r)*sig2CP, T) + matrix(data=sig2CP*r, nrow=T, ncol=T))
    }
    #exponential decay structure
    if(type==1) { 
      Vi <- bdiag(Vi, diag(sig2[k,],T) + sig2CP*(r^abs(matrix(1:T,nrow=T, ncol=T, byrow=FALSE) - matrix(1:T,nrow=T, ncol=T, byrow=TRUE))))
    }
  }
  
  #Variance matrix for all clusters
  Vall <- Vi[!is.na(Xvec),!is.na(Xvec)]
  
  
  #there will be problems if Zmat is not of full column rank
  if(rankMatrix(Zmat)[1] < ncol(Zmat)) return(NA)
  else return(solve((t(Zmat)%*%solve(Vall)%*%Zmat))[ncol(Zmat),ncol(Zmat)])
  
}

CRTVar_mMAT_omitcol <- function(ocol, Xmat, Mmat, rho0, r, type) {
  #ocol is the column to be omitted i.e. the omitted time period
  totalvar <- 1
  sig2CP <- rho0*totalvar
  sig2E <- totalvar - sig2CP
  sig2 <- sig2E/Mmat #this is now a matrix. 
  
  T <- ncol(Xmat[,-ocol])
  K <- nrow(Xmat)
  Xvec <- as.vector(t(Xmat[,-ocol]))
  stackI <- matrix(rep(diag(1,T)), nrow=K*T, ncol=T, byrow=TRUE)
  Zmat <- cbind(stackI[!is.na(Xvec),], Xvec[!is.na(Xvec)])
  
  #Variance matrix for one cluster, with decay in correlation over time
  #Vi <- diag(sig2,T) + matrix(sig2CP,nrow=T, ncol=T)
  #Constant decay var if type==0
  if(type==0) { 
    Vi <-diag(sig2[1,] +(1-r)*sig2CP, T+1) + matrix(data=sig2CP*r, nrow=(T+1), ncol=(T+1))
    Vi <- Vi[-ocol, -ocol]
  }
  #exponential decay structure
  if(type==1) { 
    Vi <- diag(sig2[1,],T+1) + sig2CP*(r^abs(matrix(1:(T+1),nrow=(T+1), ncol=(T+1), byrow=FALSE) - matrix(1:(T+1),nrow=(T+1), ncol=(T+1), byrow=TRUE)))
    Vi <- Vi[-ocol, -ocol] 
  }
  
  for(k in 2:K){
    if(type==0) { 
      nextVi <- diag(sig2[k,] +(1-r)*sig2CP, T+1) + matrix(data=sig2CP*r, nrow=(T+1), ncol=(T+1))
      nextVi <- nextVi[-ocol, -ocol]
      Vi <- bdiag(Vi, nextVi)
    }
    #exponential decay structure
    if(type==1) { 
      nextVi <-  diag(sig2[k,],T+1) + sig2CP*(r^abs(matrix(1:(T+1),nrow=(T+1), ncol=(T+1), byrow=FALSE) - matrix(1:(T+1),nrow=(T+1), ncol=(T+1), byrow=TRUE)))
      nextVi <- nextVi[-ocol, -ocol]
      Vi <- bdiag(Vi, nextVi)
    }
  }
  
  
  Vall <- Vi
  
  vartheta <- solve((t(Zmat)%*%solve(Vall)%*%Zmat))[ncol(Zmat),ncol(Zmat)]
  #[T+1, T+1]
  
  return(vartheta)
  
}



#Information content of cluster-periods for one design:
ICclusterperiod <- function(Xmat, Mmat, rho0, r, type){
  
  Xdes_colsum <- colSums(Xmat, na.rm=TRUE)
  Xdes_nasums <- colSums(is.na(Xmat) == FALSE)
  
  varmat_excl<-matrix(data=NA, nrow=nrow(Xmat), ncol=ncol(Xmat))
  varmatall <- CRTVar_mMAT(Xmat, Mmat, rho0, r, type)
  
  for(i in 1:nrow(Xmat)) {
    for(j in 1:ncol(Xmat)) {
      if(is.na(Xmat[i,j])==TRUE)  varmat_excl[i,j] <- NA
      else if(is.na(Xmat[i,j])==FALSE) {
        if(Xdes_colsum[j] == 1 & Xdes_nasums[j]==2) varmat_excl[i,j] <- 2.2772
        else   {
          Xdesij <- Xmat
          Xdesij[i,j] <- NA
          Mmatij <- Mmat
          Mmatij[i,j] <- NA
          varmat_excl[i,j] <- CRTVar_mMAT(Xdesij, Mmatij, rho0, r, type)/varmatall
        }
      }
    }
  }
  
  return(varmat_excl)
  
} 

#Information content of sequence-periods for one design:
ICseqperiod <- function(Xmat, Mmat, emat, rho0, r, type){
  
  #emat defines which rows of the Xmat define the 
  #sequences: it is an adjacency matrix with 
  ## of clusters rows and #seqs cols
  # 1 in (k,s) when k assigned to sequence s
  
  seqmat <- as.matrix(distinct(as.data.frame(Xmat)))
  seqmat_colsum <- colSums(seqmat, na.rm=TRUE)
  seqmat_nasums <- colSums(is.na(seqmat) == FALSE)
  
  varmat_excl<-matrix(data=NA, nrow=nrow(seqmat), ncol=ncol(seqmat))
  varmatall <- CRTVar_mMAT(Xmat, Mmat, rho0, r, type)
  
  for(i in 1:nrow(seqmat)) {
    for(j in 1:ncol(seqmat)) {
      if(is.na(seqmat[i,j])==TRUE)  seqmat_excl[i,j] <- NA
      else if(is.na(seqmat[i,j])==FALSE) {
        if(seqmat_colsum[j] == 1 & seqmat_nasums[j]==2) seqmatt_excl[i,j] <- 2.2772
        else   {
          Xdesij <- Xmat
          #remove all clusters assigned to this sequence
          Xdesij[emat[,i]== 1,j] <- NA
          Mmatij <- Mmat
          Mmatij[emat[,i]== 1,j] <- NA
          varmat_excl[i,j] <- CRTVar_mMAT(Xdesij, Mmatij, rho0, r, type)/varmatall
        }
      }
    }
  }
  
  return(varmat_excl)
  
}

#Information content of cluster-periods averaged over a set of designs
ICclusterperiod_ave <- function(Xmat, Mmat, Emat, rho0, r, type){
  
  #Emat is the set of cluster orderings
  #Xmat designates the design 
  #Note that this function does not yet allow for the investigation
  #of designs where the number of clusters assigned to each sequence 
  #changes
  
  icclusterperiods <- matrix(data = 0, nrow= nrow(Mmat), ncol=ncol(Mmat))
  
  for(i in 1:nrow(Emat)){
    desmat <- Xmat[as.numeric(levels(reorder(seq(1,nrow(Xmat),1),Emat[i,]))), ]
    
    ICcp_thisorder <- ICclusterperiod(desmat, Mmat,  rho0, r, type)    
    #Need to arrange these rows so that cluster 1 is first etc.
    #Otherwise won't be averaging over orderings within clusters
    #(will instead be getting averaging within sequences)
    
    # print(icclusterperiods)
    icclusterperiods <- icclusterperiods + ICcp_thisorder#[Emat[i,], ]
    
  }
  
  return(icclusterperiods/nrow(Emat))
  
} 

#Information content of sequence-periods averaged over a set of designs
ICseqperiod_ave <- function(Xmat, Mmat, Emat, emat, rho0, r, type){
  
  #emat defines the number of clusters assigned to each sequence
  #Emat is the set of cluster orderings
  #Xmat designates the design 
  #Note that this function does not yet allow for the investigation
  #of designs where the number of clusters assigned to each sequence 
  #changes
  
  seqmat <- as.matrix(distinct(as.data.frame(Xmat)))

  icseqperiods <- matrix(data = 0, nrow= nrow(seqmat), ncol=ncol(seqmat))
  
  for(i in 1:nrow(Emat)){
    desmat <- Xmat[as.numeric(levels(reorder(seq(1,nrow(Xmat),1),Emat[i,]))), ]
    reord_emat <- emat[as.numeric(levels(reorder(seq(1,nrow(Xmat),1),Emat[i,]))), ]
    
    icseqperiods <- icseqperiods + ICseqperiod(desmat, Mmat, reord_emat,  rho0, r, type)
    
    
  }
  
  return(icseqperiods/nrow(Emat))
  
} 


#Information content of entire clusters for one design:
ICcluster <- function(Xmat, Mmat, rho0, r, type){
  
  varmat_excl<-matrix(data=NA, nrow=nrow(Xmat), ncol=ncol(Xmat))
  
  varmatall <- CRTVar_mMAT(Xmat, Mmat, rho0, r, type)
  
  for(i in 1:nrow(Xmat)) {
    Xdesij <- Xmat
    Xdesij[i,] <- NA
    Mmatij <- Mmat
    Mmatij[i,] <- NA
    varmat_excl[i,] <- CRTVar_mMAT(Xdesij, Mmatij, rho0, r, type)/varmatall
  }
  
  return(varmat_excl)
  
} 

#IC of clusters averaged over a set of designs
ICcluster_ave <- function(Xmat, Mmat, Emat, rho0, r, type){
  
  #Emat is the set of cluster orderings
  #Xmat designates the set of sequences
  #Note that this function does not yet allow for the investigation
  #of designs where the number of clusters assigned to each sequence 
  #changes
  
  icclusterperiods <- matrix(data = 0, nrow= nrow(Mmat), ncol=ncol(Mmat))

  for(i in 1:nrow(Emat)){
    desmat <- Xmat[as.numeric(levels(reorder(seq(1,nrow(Xmat),1),Emat[i,]))), ]
    
    
    ICcp_thisorder <- ICcluster(desmat, Mmat,  rho0, r, type)
    #Need to arrange these rows so that cluster 1 is first etc.
    #Otherwise won't be averaging over orderings within clusters
    #(will instead be getting averaging within sequences)
    
    icclusterperiods <- icclusterperiods + ICcp_thisorder 
    
  }
  
  return(icclusterperiods/nrow(Emat))
  
} 

#Information content of sequences for one design:
ICseq <- function(Xmat, Mmat, emat, rho0, r, type){
  
  #emat defines which rows of the Xmat define the 
  #sequences: it is an adjacency matrix with 
  ## of clusters rows and #seqs cols
  # 1 in (k,s) when k assigned to sequence s
  
  seqmat <- as.matrix(distinct(as.data.frame(Xmat)))
  
  #If there are any periods with only 2 treatment sequences, with 
  #differential exposure, cannot calculate the information content 
  #of either of these two sequences.
  #Need to flag that the cells in these periods MUST be included and
  #thus do not have an information content.
  seqmat_colsum <- colSums(seqmat, na.rm=TRUE)
  seqmat_nasums <- colSums(is.na(seqmat) == FALSE)
  
  
  varmat_excl<-matrix(data=NA, nrow=nrow(seqmat), ncol=ncol(seqmat))
  
  varmatall <- CRTVar_mMAT(Xmat, Mmat, rho0, r, type)
  
  #Code in for loops, but re-do this later
  for(i in 1:nrow(seqmat)) {
    Xdesij <- Xmat
    #remove all clusters assigned to this sequence
    Xdesij[emat[,i]== 1,] <- NA
    Mmatij <- Mmat
    Mmatij[emat[,i]== 1,] <- NA
    varmat_excl[i,] <- CRTVar_mMAT(Xdesij, Mmatij, rho0, r, type)/varmatall
  }
  
  
  return(varmat_excl)
  
}

#Information content of sequences averaged over a set of designs
ICseq_ave <- function(Xmat, Mmat, Emat, emat, rho0, r, type){
  
  #emat defines the number of clusters assigned to each sequence
  #Emat is the set of cluster orderings
  #Xmat designates the design 
  #Note that this function does not yet allow for the investigation
  #of designs where the number of clusters assigned to each sequence 
  #changes
  
  
  seqmat <- as.matrix(distinct(as.data.frame(Xmat)))
  #need to not code in for loops!
  icseqs <- matrix(data = 0, nrow= nrow(seqmat), ncol=ncol(seqmat))
  
  
  for(i in 1:nrow(Emat)){
    desmat <- Xmat[as.numeric(levels(reorder(seq(1,nrow(Xmat),1),Emat[i,]))), ]
    reord_emat <- emat[as.numeric(levels(reorder(seq(1,nrow(Xmat),1),Emat[i,]))), ]
    icseqs <- icseqs + ICseq(desmat, Mmat, reord_emat,  rho0, r, type)
    
    
  }
  
  return(icseqs/nrow(Emat))
  
} 

#Information content of entire periods for one design:
ICperiod <- function(Xmat, Mmat, rho0, r, type){
  
  varmat_excl<-matrix(data=NA, nrow=nrow(Xmat), ncol=ncol(Xmat))
  
  varmatall <- CRTVar_mMAT(Xmat, Mmat, rho0, r, type)
  
  #Code in for loops, but re-do this later
  for(i in 1:ncol(Xmat)) {
    
    varmat_excl[,i] <- CRTVar_mMAT_omitcol(i, Xmat, Mmat, rho0, r, type)/varmatall
    
  }
  
  return(varmat_excl)
  
} 

#IC of periods averaged over a set of designs
ICperiod_ave <- function(Xmat, Mmat, Emat, rho0, r, type){
  
  #Emat is the set of cluster orderings
  #Xmat designates the design 
  #Note that this function does not yet allow for the investigation
  #of designs where the number of clusters assigned to each sequence 
  #changes
  
  icclusterperiods <- matrix(data = 0, nrow= nrow(Mmat), ncol=ncol(Mmat))
  
  
  for(i in 1:nrow(Emat)){
    desmat <- Xmat[as.numeric(levels(reorder(seq(1,nrow(Xmat),1),Emat[i,]))), ]
    ICcp_thisorder <- ICperiod(desmat, Mmat,  rho0, r, type)
    #Need to arrange these rows so that cluster 1 is first etc.
    #Otherwise won't be averaging over orderings within clusters
    #(will instead be getting averaging within sequences)
    
    icclusterperiods <- icclusterperiods + ICcp_thisorder
    
  }
  
  return(icclusterperiods/nrow(Emat))
  
} 


#information content plotting function
#takes the information contents in a matrix
#and plots them in the usual style
ICplot <- function(ICmatrix, myylab = "Sequence"){
  
  ICmatrix<-round(ICmatrix, 4)
  melted_varmatexcl <- melt(ICmatrix)
  names(melted_varmatexcl)[names(melted_varmatexcl)=="Var1"] <- "Sequence"
  names(melted_varmatexcl)[names(melted_varmatexcl)=="Var2"] <- "Period"
  
  color_palette <-colorRampPalette(c( "yellow", "red"))(length(table(ICmatrix)))
  if(sum(melted_varmatexcl$value==2.2772, na.rm=TRUE) > 0)    color_palette[length(table(melted_varmatexcl$value))]<- "#000000"
  
  T <- ncol(ICmatrix)
  K <- nrow(ICmatrix)
  
  ggplot(data = melted_varmatexcl, aes(x=Period, y=Sequence, fill = factor(value))) + 
    geom_tile( colour = "grey50") +
    scale_y_reverse(breaks=c(1:K)) +
    ylab(myylab) +
    scale_x_continuous(breaks=c(1:T)) +
    theme(panel.grid.minor = element_blank()) +      
    geom_text(aes(Period, Sequence, label = round(value,4)), color = "black", size = 4) +
    scale_fill_manual(values = color_palette, breaks=levels(melted_varmatexcl$value)[seq(90, 150, by=5)]) 
  
  
}

##Plot to include cluster-period sizes on the IC heat plot
ICplot_clustersizes <- function(ICmatrix, sizematrix, myylab = "Sequence"){
  
  ICmatrix<-round(ICmatrix, 4)
  melted_varmatexcl <- melt(ICmatrix)
  names(melted_varmatexcl)[names(melted_varmatexcl)=="Var1"] <- "Sequence"
  names(melted_varmatexcl)[names(melted_varmatexcl)=="Var2"] <- "Period"
  
  melted_size <- melt(sizematrix)
  
  
  color_palette <-colorRampPalette(c( "yellow", "red"))(length(table(ICmatrix)))
  if(sum(melted_varmatexcl$value==2.2772, na.rm=TRUE) > 0)    color_palette[length(table(melted_varmatexcl$value))]<- "#000000"
  
  T <- ncol(ICmatrix)
  K <- nrow(ICmatrix)
  
  ggplot(data = melted_varmatexcl, aes(x=Period, y=Sequence, fill = factor(value))) + 
    geom_tile( colour = "grey50") +
    scale_y_reverse(breaks=c(1:K)) +
    ylab(myylab) +
    scale_x_continuous(breaks=c(1:T)) +
    theme(panel.grid.minor = element_blank()) +      
    geom_text(aes(Period, Sequence, label = round(value,4)), color = "black", size = 4) +
    geom_text(aes(melted_size$Var2, melted_size$Var1, label = paste("m =", melted_size$value)), 
              color = "black", size = 2, nudge_y=-0.35) +
    scale_fill_manual(values = color_palette, breaks=levels(melted_varmatexcl$value)[seq(90, 150, by=5)]) 
  
  
}


#Profile plots for entire clusters, periods, sequences:
Profileplot <- function(ICmatrix, dim,  myxlab = "Sequence"){
  
  #if dim = 1, first row considered (i.e. periods)
  #if dim = 2, first column considered (i.e. periods)
  if(dim == 1) mydata <- ICmatrix[1,]
  if(dim == 2) mydata <- ICmatrix[,1]
  
  myseq <- seq(1, length(mydata)) 
  mydataset <- cbind(mydata, myseq)
  mydataset <- as.data.frame(mydataset)
  
  ggplot(mydataset, aes(x=myseq, y=mydata)) +
    geom_line() + geom_point(size=2) +
    labs(x=myxlab, y="Information content", colour="", linetype="") + 
    scale_x_continuous(breaks=seq(1,length(mydata),1)) 
}

#Function that performs one simulation of cluster sizes
#and produces the cluster-period ICs for that simulation
ICseqperiod_sim <- function(myalpha, mybeta, Totsize, Xmat, emat, rho0, r, type){
  
  #myalpha, mybeta: parameters of the gamma distribution
  #Totsize: total size of the study
  #Xmat: design matrix
  
  sim01 <- rgamma(n=nrow(Xmat), shape = myalpha, scale = mybeta)
  Mmat <- matrix(data =rep(sim01*Totsize/sum(sim01), each = ncol(Xmat)), 
                 nrow=nrow(Xmat), ncol=ncol(Xmat), byrow=T)
  
  return(ICseqperiod(Xmat, Mmat, emat, rho0, r, type))
}

#Variances of all possible allocations
var_alldesigns <- function(Xmat, Mmat, Emat, emat, rho0, r, type){
  
  #emat defines the number of clusters assigned to each sequence
  #Emat is the set of cluster orderings
  #Xmat designates the design 
  #Note that this function does not yet allow for the investigation
  #of designs where the number of clusters assigned to each sequence 
  #changes
  
  
  seqmat <- as.matrix(distinct(as.data.frame(Xmat)))
  #need to not code in for loops!
  icseqperiods <- matrix(data = 0, nrow= nrow(seqmat), ncol=ncol(seqmat))
  
  myvars <- matrix(data=NA, nrow=nrow(Emat), ncol=1)
  
  for(i in 1:nrow(Emat)){
    desmat <- Xmat[as.numeric(levels(reorder(seq(1,nrow(Xmat),1),Emat[i,]))), ]
    reord_emat <- emat[as.numeric(levels(reorder(seq(1,nrow(Xmat),1),Emat[i,]))), ]
    
    myvars[i,1] <- CRTVar_mMAT(desmat, Mmat, rho0, r, type)
    
    
  }
  
  return(myvars)
  
} 
