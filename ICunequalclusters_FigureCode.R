##################################################
# Information content of stepped wedge and other 
# centrosymmetric longitudinal cluster randomised 
# trial designs with unequal cluster-period sizes
# J Kasza, AB Forbes
#
# Code to replicate figures in paper
# J Kasza, 2019-10-01
# jessica.kasza@monash.edu
##################################################

#Required libraries:
library(Matrix)
library(partitions)
library(dplyr)
library(multicool)
library(reshape2)
library(ggplot2)
library(gridExtra)

source("ICunequalclusters_Functions.R", local=TRUE)

#################################################
# Hill et al. (2015) example: a stepped wedge
# trial assessing an intervention on fall rates 
# in rehabilitation units
#################################################

#Input the cluster sizes:
clustersizes <- matrix(data=NA, nrow = 8, ncol= 5)
clustersizes[1,] <- c(66,52,112,119,117)
clustersizes[2,] <- c(36,37,42,45,33)
clustersizes[3,] <- c(42, 50, 65, 73, 67)
clustersizes[4,] <- c(71, 65, 73, 80, 68)
clustersizes[5,] <- c(97, 103, 121, 138, 124)
clustersizes[6,] <- c(81, 88, 82, 58, 73)
clustersizes[7,] <- c(84, 88, 149, 98, 88)
clustersizes[8,] <- c(128, 184, 189, 161, 159)

#Sum over clusters assigned to each sequence
clustersizes_total <- matrix(data=NA, nrow=4, ncol=5)
clustersizes_total[1,] <- clustersizes[1,] + clustersizes[2,]
clustersizes_total[2,] <- clustersizes[3,] + clustersizes[4,]
clustersizes_total[3,] <- clustersizes[5,] + clustersizes[6,]
clustersizes_total[4,] <- clustersizes[7,] + clustersizes[8,]

#Row and column mean cluster-period sizes:
rowMeans(clustersizes)
colMeans(clustersizes)

#The design matrix: two clusters assigned to each
# of the four sequences of a 5-period stepped wedge
desmat <- matrix(rep(SWdesmat(5),each=2),nrow=2*4)

#emat: The assignment of clusters to sequences
emat <- matrix(data=0, nrow= 8, ncol=4)
emat[1,1] <- 1
emat[2,1] <- 1
emat[3,2] <- 1
emat[4,2] <- 1
emat[5,3] <- 1
emat[6,3] <- 1
emat[7,4] <- 1
emat[8,4] <- 1

##Design parameters:
rho0 <- 0.14
myr <- 1
decaytype <- 0 #indicates a Hussey and Hughes correlation structure

##How many possible allocations of clusters to 
##sequences are there? Each sequence has two 
##clusters assigned to it.
  #Initialise the permutations of clusters:
  m1 <- initMC(seq(1,8,1))
  myperms <- allPerm(m1)
  #Within each pair of clusters, order from smallest to largest 
  #and then get rid of duplicate rows. 
  uniqueSWs <- sortbym(myperms, 2)
  uniqueSWs <- as.matrix(uniqueSWs)
  dim(uniqueSWs)
  #[1] 2520    8

#################################################################
#Figures 3 and 4

#Information content of design components for the actual trial:
  #cluster-periods
  temp <- ICclusterperiod(desmat, clustersizes, rho0, myr, decaytype)
  Hill_cp <- ICplot_clustersizes(temp, clustersizes, myylab = "Cluster")
  #sequence-periods
  temp <- ICseqperiod(desmat, clustersizes, emat, rho0, myr, decaytype)
  Hill_seqp <- ICplot_clustersizes(temp, clustersizes_total)
  #clusters
  temp <- ICcluster(desmat, clustersizes, rho0, myr, decaytype)
  Hill_c <- ICplot(temp, myylab = "Cluster")
  Hill_c_profile <- Profileplot(temp, 2, myxlab="Cluster")
  #periods
  temp <- ICperiod(desmat, clustersizes, rho0, myr, decaytype)
  Hill_p <- ICplot(temp, myylab = "Cluster")
  Hill_p_profile <- Profileplot(temp, 1, myxlab="Period")
  #sequences
  temp <- ICseq(desmat, clustersizes, emat, rho0, myr, decaytype)
  Hill_s <- ICplot(temp, myylab = "Sequence")
  Hill_s_profile <- Profileplot(temp, 2, myxlab="Sequence")

  #Save the plots:
  Figure3 <- grid.arrange(arrangeGrob(Hill_cp, 
                                       Hill_seqp, nrow=2, ncol=1))
  ggsave("Figure3.jpeg", plot=Figure3, width= 14, height=21, units="cm", dpi=600)
  ggsave("Figure3.pdf", plot=Figure3, width= 14, height=21, units="cm", dpi=600)

  Figure4 <- grid.arrange(arrangeGrob(Hill_c_profile, Hill_s_profile,
                                     Hill_p_profile, nrow=1, ncol=3))
  ggsave("Figure4.jpeg", plot=Figure4, width= 14, height=10, units="cm", dpi=600)
  ggsave("Figure4.pdf", plot=Figure4, width= 14, height=10, units="cm", dpi=600)

#################################################################

#################################################################  
#Figures 5 and 6: averaging over all valid arrangements:
##NOTE: THESE COMMANDS CAN TAKE SOME TIME TO RUN  
  
  #Cluster-periods
  CPaverages <- ICclusterperiod_ave(desmat, clustersizes, uniqueSWs, rho0, myr, decaytype)
  CPaveplot <- ICplot_clustersizes(CPaverages, clustersizes,myylab = "Cluster") 
  #sequence-periods
  SeqPaverages <- ICseqperiod_ave(desmat, clustersizes, uniqueSWs, emat, rho0, myr, decaytype) 
  SEQaveplot <- ICplot(SeqPaverages,  myylab = "Sequence")

  Figure5 <- grid.arrange(arrangeGrob(CPaveplot, 
                                          SEQaveplot, nrow=2, ncol=1))
  ggsave("Figure5.jpeg", plot=Figure5, width= 14, height=21, units="cm", dpi=600)
  ggsave("Figure5.pdf", plot=Figure5, width= 14, height=21, units="cm", dpi=600)

  #sequences, clusters, and periods
  Seqaverages <- ICseq_ave(desmat, clustersizes, uniqueSWs, emat, rho0, myr, decaytype)
  Clusteraverages <- ICcluster_ave(desmat, clustersizes, uniqueSWs,  rho0, myr, decaytype)
  Periodaverages <- ICperiod_ave(desmat, clustersizes, uniqueSWs,  rho0, myr, decaytype)

  Hill_s_profileAVE <- Profileplot(Seqaverages, 2, myxlab="Sequence")
  Hill_c_profileAVE <- Profileplot(Clusteraverages, 2, myxlab="Cluster")
  Hill_p_profileAVE <- Profileplot(Periodaverages, 1, myxlab="Period")

  Figure6 <- grid.arrange(arrangeGrob(Hill_c_profileAVE, Hill_s_profileAVE,
                                        Hill_p_profileAVE, nrow=1, ncol=3))
  ggsave("Figure6.jpeg", plot=Figure6, width= 14, height=10, units="cm", dpi=600)
  ggsave("Figure6.pdf", plot=Figure6, width= 14, height=10, units="cm", dpi=600)
#################################################################  
  
#################################################################  
#Figure 7: simulation study
# Use the distribution for sample sizes considered in 
# Martin et al (2019) BMC Med Res Methodol
#Generate cluster-period sizes from a Gamma distribution
#assume equality of period sizes within each cluster
#Need to scale to ensure total number of subjects the same in each 
#simulation (a la Martin et al)
  
#Use the same design as the Hill et al trial
#4 sequence 5 period design, with 8 clusters and 
#2 per sequence.
  
  #Total study size: 
  sum(clustersizes)
  #3606
  mean(clustersizes)
  #[1] 90.15
  sd(clustersizes)
  #[1] 40.55863
  
  #Coefficient of variation: 
  mycv <- sd(clustersizes)/mean(clustersizes)
  #0.45
  
  #Using formulas in Martin et al
  myalpha <- (mean(clustersizes)^2)/(sd(clustersizes)^2)
  mybeta <- (sd(clustersizes)^2)/mean(clustersizes)
  Totsize <- 3600

set.seed(142934)
SPreplicates <- replicate(1000, ICseqperiod_sim(myalpha, mybeta, Totsize, desmat, emat, rho0, myr, decaytype), simplify=FALSE)

SPreps_red <- Reduce('+', SPreplicates)/length(SPreplicates)

#How does this compare to the scenario when all cluster-periods of same size
#assume each cp cell of size 90  (90*5*8 = 3600)
SPreps_ave <- ICseqperiod(desmat, matrix(data=90, nrow=8, ncol=5), emat, rho0, myr, decaytype)

SP_simsplot <- ICplot(SPreps_red, myylab = "Sequence")
SP_aveplot <- ICplot(SPreps_ave, myylab = "Sequence")

Figure7 <- grid.arrange(arrangeGrob(SP_simsplot, 
                                     SP_aveplot, nrow=2, ncol=1))
ggsave("Figure7.jpeg", plot=Figure7, width= 14, height=21, units="cm", dpi=600)
ggsave("Figure7.pdf", plot=Figure7, width= 14, height=21, units="cm", dpi=600)
#################################################################  

#################################################################  
#Figure 9: Incomplete designs for the Hill et al. (2015) trial

#The design matrices for the incomplete designs:
  incompleteA <-  matrix(rep(SWdesmat(5),each=2),nrow=2*4)
  incompleteA[1:2, 3:5] <- NA
  incompleteA[3:4, 1] <- NA
  incompleteA[3:4, 4:5] <- NA
  incompleteA[5:6, 1:2] <- NA
  incompleteA[5:6, 5] <- NA
  incompleteA[7:8, 1:3] <- NA

  incompleteB <- incompleteA
  incompleteB[1:2,5] <- 1
  incompleteB[7:8,1] <- 0

  incompleteC <- incompleteB
  incompleteC[3:4,4] <- 1
  incompleteC[5:6,2] <- 0

  incompleteD <- incompleteC
  incompleteD[1:2,5] <- NA
  incompleteD[1:2,3] <- 1
  incompleteD[7:8,1] <- NA
  incompleteD[7:8,3] <- 0

#Variances for the given allocation for these incomplete designs
completevar <- CRTVar_mMAT(matrix(rep(SWdesmat(5),each=2),nrow=2*4), clustersizes, rho0, myr, decaytype)
incompleteAvar <- CRTVar_mMAT(incompleteA, clustersizes, rho0, myr, decaytype)
incompleteBvar <- CRTVar_mMAT(incompleteB, clustersizes, rho0, myr, decaytype)
incompleteCvar <- CRTVar_mMAT(incompleteC, clustersizes, rho0, myr, decaytype)
incompleteDvar <- CRTVar_mMAT(incompleteD, clustersizes, rho0, myr, decaytype)

varratios<- c(incompleteAvar/completevar, incompleteBvar/completevar,
              incompleteCvar/completevar, incompleteDvar/completevar)

#Averages over all possible allocations of clusters to sequences
complete_allvar <- var_alldesigns(matrix(rep(SWdesmat(5),each=2),nrow=2*4), clustersizes, uniqueSWs, emat, rho0, myr, decaytype)
incompleteA_allvar <- var_alldesigns(incompleteA, clustersizes, uniqueSWs, emat, rho0, myr, decaytype)
incompleteB_allvar <- var_alldesigns(incompleteB, clustersizes, uniqueSWs, emat, rho0, myr, decaytype)
incompleteC_allvar <- var_alldesigns(incompleteC, clustersizes, uniqueSWs, emat, rho0, myr, decaytype)
incompleteD_allvar <- var_alldesigns(incompleteD, clustersizes, uniqueSWs, emat, rho0, myr, decaytype)

varratio_allvar_1 <- incompleteA_allvar/complete_allvar
varratio_allvar_2 <- incompleteB_allvar/complete_allvar
varratio_allvar_3 <- incompleteC_allvar/complete_allvar
varratio_allvar_4 <- incompleteD_allvar/complete_allvar
design <- rep(c(1,2,3,4), each=2520)
myratios <- as.data.frame(cbind(design, rbind(varratio_allvar_1, varratio_allvar_2,
                                              varratio_allvar_3, varratio_allvar_4)))

#Now plot the data....
Figure9<-   ggplot(myratios, aes(y=V2, x=as.factor(design))) +
  geom_boxplot() + theme_minimal() +
  scale_x_discrete(labels=c("A","B","C","D")) +
  labs(y="Ratio of variances: incomplete/complete", x="Design") + 
  theme(legend.position ="bottom") +
  theme(legend.key.size = unit(0.75, "cm")) +
  theme(legend.text=element_text(size=10)) 


ggsave("Figure9.jpeg", plot=Figure9, width= 10, height=14, units="cm", dpi=600)
ggsave("Figure9.pdf", plot=Figure9, width= 10, height=14, units="cm", dpi=600)

#Alternative way to plot data: overlapping histograms
#Figure9<-   ggplot(myratios, aes(V2, fill=as.factor(design))) +
#  geom_histogram(binwidth= 0.025) + theme_minimal() +
#  scale_fill_manual(values = c("#BDB8AD","#382119","#C6D4E1","#44749D"), labels=c("A", "B", "C", "D"), 
#                    name="Design: ") +
#  labs(y="Count", x="Ratio of variances: incomplete/complete", colour="", linetype="") + 
#  theme(legend.position ="bottom") +
#  theme(legend.key.size = unit(0.75, "cm")) +
#  theme(legend.text=element_text(size=10)) 