##################################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################################################
# Start
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
##################################################################################################
# NanoString Classifiation Project
# Dr Reza Rafiee, Research Associate at NICR, Newcastle University, 2014-2016
# Aim: Reliably classify the samples into four subgroups using 22 Genes
# Because of the copyright law, You should not redistribute this code in any way.
##################################################################################################
# Input
# 36 features (gene expression) of medulloblastoma (brain tumour) samples (exctracted from RNA materials)
#
# Reading NanoString test data  =====>    matrix Sample_test: 
#                       
#                         Sample name 
#
#                      NMB147 ... NMB154     
#   Gene name     _                      _
#    Gene 1      |      2915  ...  5875   |
#    Gene 2      |                        |
#    Gene 3      |                        |
#      ...       |                        | 
#      ...       |           ...          |
#      ...       |                        |  
#    Gene 36     |_                      _|
#
#      

# WNT genes (5)  _
# 1: DKK2         |
# 2: EMX2         |
# 3: GAD1         |
# 4: TNC          |
# 5: WIF1        _|

# ----------------     
# SHH genes (4)  _
# 6: ATOH1        |
# 7: EYA1         |
# 8: HHIP         |
# 9: SFRP1       _|

# ----------------     
# Grp3 genes (5) _
# 10: GABRA5      |
# 11: IMPG2       |
# 12: MAB21L2     |
# 13: NPR3        |
# 14: NRL        _|

# ----------------     
# Grp4 genes (5) _
# 15: EOMES       |
# 16: KCNA1       |
# 17: KHDRBS2     |
# 18: RBM24       |
# 19: UNC5D      _|

# ----------------     
# Negative control genes (8)
# 20: NEG_A
# 21: NEG_B
# 22: NEG_C
# 23: NEG_D
# 24: NEG_E
# 25: NEG_F
# 26: NEG_G
# 27: NEG_H
# ----------------     
# Positive control genes (6)
# 28: POS_A
# 29: POS_B
# 30: POS_C
# 31: POS_D
# 32: POS_E
# 33: POS_F
# ----------------     
# Housekeeping genes (3)
# 34: GAPDH
# 35: ACTB
# 36: LDHA

# ==============================================================================
# Input csv file, raw data of nanoString gene expression
# The column header (column names) must be sample names, row names must be 36 gene names 
###############################################################################################
# Output
# Classifier confidence and subgroup labels for the samples
################################################################################################
#STARTSTARTSTARTSTARTSTARTSTARTSTARTSTARTSTARTSTARTSTARTSTARTSTARTSTARTSTARTSTARTSTARTSTARTSTART
################################################################################################
################################## Normalisation procedure #####################################
setwd("~/NanoString")

# load libraries
library(NanoStringNorm) # normalising nanoString raw data
library(e1071)          # Support verctor machine (SVM) classifier
library(parallel)       # mclapply speeds up probability estimation 


# Reading an input csv file including raw gene expression data of 19 genes as well as conrtol genes
NanostringSamplesGenes <- as.data.frame(read.csv("~/NanoString/81RawNanoString08112016_TestSample.csv",header=T))

# Subgroup medulloblastoma genes 
# WNT_genes <- c("DKK2","EMX2","GAD1","TNC","WIF1")
# SHH_genes <- c("ATOH1","EYA1","HHIP","SFRP1")
# Grp3_genes <- c("GABRA5","IMPG2","MAB21L2","NPR3","NRL")
# Grp4_genes <- c("EOMES","KCNA1","KHDRBS2","RBM24","UNC5D")

All_genes_ordered <- c("DKK2","EMX2","GAD1","TNC","WIF1","ATOH1","EYA1","HHIP","SFRP1","GABRA5","IMPG2","MAB21L2","NPR3","NRL","EOMES","KCNA1","KHDRBS2","RBM24","UNC5D",
                  "NEG_A", "NEG_B", "NEG_C", "NEG_D", "NEG_E", "NEG_F", "NEG_G", "NEG_H", "POS_A", "POS_B", "POS_C", "POS_D", "POS_E", "POS_F", "GAPDH", "ACTB", "LDHA")

# Order the subgroup-based genes of the input file
rownames(NanostringSamplesGenes) <- NanostringSamplesGenes$Name
NanostringSamplesGenes <- NanostringSamplesGenes[All_genes_ordered,]

# Separating data from annotation
NanostringSamples.anno <- NanostringSamplesGenes[,c(1:3)]
NanostringSamples.data <- NanostringSamplesGenes[,-c(1:3)]

# Normalising raw data
NanoStriang_mRNA_norm_tech_bio <- NanoStringNorm(
  x = NanostringSamples.data, 
  anno = NanostringSamples.anno,
  CodeCount = 'geo.mean',
  Background = 'mean.2sd',
  SampleContent = 'housekeeping.geo.mean',
  round.values = TRUE,
  take.log = TRUE,  #log2 transformed
  return.matrix.of.endogenous.probes = TRUE
);

# # Chnage return.matrix.of.endogenous.probes value into FALSE
# Plot.NanoStringNorm(
#   x = NanoStriang_mRNA_norm_tech_bio,
#   label.best.guess = TRUE,
#   #plot.type = c('volcano'),
#   plot.type = c('cv', 'mean.sd'),
#   title = FALSE
# )


# NanoString.mRNA.norm.tech.bio # results
# Normalising gene expression data into range of [0 1]
Nano1_log2_Normalised_4 <- (NanoStriang_mRNA_norm_tech_bio-min(NanoStriang_mRNA_norm_tech_bio))/(max(NanoStriang_mRNA_norm_tech_bio)-min(NanoStriang_mRNA_norm_tech_bio))
Nano1_log2_Normalised_4 <- as.data.frame(Nano1_log2_Normalised_4)

Sample_test <- Nano1_log2_Normalised_4

################################################################################################
############################### End of normalisation procedure #################################
################################################################################################
# 22 January 2016
GoldCohort.Threshold <- 0.6847321
#GoldCohort.Threshold <- 0.1
# Fitting of the distribution ' beta ' by matching moments 
# Parameters : 
#   estimate
# shape1 4.3858188
# shape2 0.3586241

##########################################################################
# Each col is a Nanostring sample, Row is a Gene value
# if already normalised
# new normalising test, 20th May 2016 #####################################################################################
#Sample_test <- as.data.frame(read.csv("112TorontoNanoStringCohort_NewNormalising_20thMay2016.csv",header=T,row.names=1)) # 20th May 2016
#RefSub_RNA_Nano <- as.integer(c(rep.int(1,16),rep.int(2,31),rep.int(3,20),rep.int(4,45)))# WNT=16, SHH=31, Grp4=20, Grp3=45
# new normalising test, 20th May 2016 #####################################################################################
#rn_Sample_test <- colnames(Sample_test)
#Sample_test <- t(apply(Sample_test,1,as.numeric))
#colnames(Sample_test) <- rn_Sample_test

NanoString_test_19 <- Sample_test
Total.No.of.Samples <- ncol(NanoString_test_19)
#############################################################
#############################################################
#############################################################
# Training set, n=101 (readcount log2 transformed), 

TrainingsetRNA_seq19 <- as.data.frame(as.matrix(read.delim("101RNA_seq19gene_log2_normalised01_5thApril2016.csv",sep=",",header=T)))
sam_names <- as.character(TrainingsetRNA_seq19[,1])
TrainingsetRNA_seq19 <- TrainingsetRNA_seq19[,-1]
rownames(TrainingsetRNA_seq19) <- sam_names

subgroup_labels <- as.integer(c(rep.int(1,11),rep.int(2,24),rep.int(3,33),rep.int(4,33)))# WNT=11, SHH=24, Grp3=33, Grp4=33
y_training <- factor(subgroup_labels)
# #############################################################
#------------------------------------
# Adding the data, time, the number of samples in "file.name" variable
var_date <- date()
#var_date <- gsub(" ", "_",var_date) 
var_date <- gsub(":","",var_date)
file.name <- paste(Total.No.of.Samples," samples, ") # 
#file.name <- gsub(" ", "",file.name)
#------------------------------------
#############################################################

### n=101, training set, 13th May 2016
cost_1 <- 4.4 
gamma_1 <- 0.01 
kernel_1 <- "radial"

x <- 1000 ## Number of iterations

# for checking the performance, 27th June 2016

amount <- round(0.9*nrow(TrainingsetRNA_seq19))

sel2<- lapply(1:x, function(i) {
  set.seed(i)
  sample(1:nrow(TrainingsetRNA_seq19), amount, replace=F)
})


## MB this bit causes a delay
#incProgress(0.10, detail = paste("Assessing confidence of subgroup calls: stage 2"))
Radial_SVM <- mclapply(1:x, 
                       mc.cores=16,
                       function(i)  svm(x = TrainingsetRNA_seq19[sel2[[i]],],
                                        y = y_training[sel2[[i]]], scale = F, 
                                        tolerance = 0.00001, type = "C-classification",
                                        kernel = kernel_1,cost = cost_1,
                                        gamma=gamma_1, probability = T,
                                        seed=i)
)


## prediction 
#incProgress(0.10, detail = paste("Assessing confidence of subgroup calls: stage 3"))
Radial_test <- mclapply(1:x,
                        mc.cores=16,
                        function(i) predict(Radial_SVM[[i]],
                                            newdata=t(NanoString_test_19), 
                                            decision.values = T,
                                            probability = T)
)

#incProgress(0.10, detail = paste("Assessing confidence of subgroup calls: stage 5"))
prob.test <- (lapply(1:x,
                     function(i) attr(Radial_test[[i]], "probabilities")))

####################################### Creating Pobes2 #############################################
k <- FALSE
for (j in 1:x) # the number of iterations
{
  for (i in 1:4) # 4 subgroups
  {
    
    predProbTemp <-prob.test[[j]] # j iteration
    predProbTemp <- predProbTemp[,c("1", "2", "3","4")] # order the matrix based on the subgroup orders
    
  }
  if (k == FALSE) # Making defult tables
  {
    
    predProbabilities <- matrix(ncol = 4, nrow =nrow(predProbTemp)*x, 0.0)
    predProbabilities <- predProbTemp 
    colnames(predProbabilities) <- c("WNT","SHH", "Grp3", "Grp4")
    
    k <- TRUE
  }
  else
  {
    #Adding other iteration probabilities to the created table in the ordered columns
    predProbabilities <- rbind(predProbabilities,predProbTemp)
  }
}

probs2 <- matrix(ncol=nrow(predProbTemp),nrow=x,0.0)
colnames(probs2) <- rownames(predProbTemp)

for (ttt in 1:nrow(predProbTemp)) # number of samples
{
  mmm <- matrix(ncol = 4, nrow =x, 0.0)
  colnames(mmm) <- c("WNT","SHH", "Grp3", "Grp4")
  gg <- 0
  for (fftt in 1:x)
  {
    gg <- gg + 1
    mmm[gg,] <- predProbabilities[ttt+nrow(predProbTemp)*(fftt-1),]
    #predProbabilities[1+3*0,1+3*1,1+3*2,1+nrow(predProbTemp)*(x-1)] # for the first sample, n=3
  }
  
  ProbSubgroup <- apply(mmm[,1:4],1,max)
  probs2[,ttt] <- ProbSubgroup
  
}
####################################### creating pobes2 #############################################

# Optimising the Gamma 
#I would suggest the following theoretical guidance. 
#When you are using Gaussian RBF kernel, your separating surface will be based on a combination of bell-shaped surfaces centered
#at each support vector. The width of each bell-shaped surface will be inversely proportional to γ. 
#If this width is smaller than the minimum pair-wise distance for your data, you essentially have overfitting.
#If this width is larger than the maximum pair-wise distance for your data, all your points fall into one class and 
#you don't have good performance either. So the optimal width should be somewhere between these two extremes

#it is essentially data dependent. Grid search (over log-transformed hyper-parameters) is a very good method if you only 
#have a small number of hyper-parameters to tune, but don't make the grid resolution too fine or
#you are likely to over-fit the tuning criterion. For problems with a larger number of kernel parameters,
#I find the Nelder-Mead simplex method works well.

# Optimising the Gamma and C (parameters of the SVM classifier)
#--------------------------------------------------------

i=1234

model <- svm(TrainingsetRNA_seq19,y_training,scale = F, tolerance = 0.00001, type = "C-classification", kernel = kernel_1,cost = cost_1, gamma=gamma_1, probability = T, seed=i)  

test.pred <- predict(object=model, newdata=t(NanoString_test_19), probability=TRUE)
prob.test <- signif(attr(test.pred, "probabilities"), digits=2)
maxProbs <- apply(prob.test,1,max)

###########################################################################################################
#maxProbsWhich <- factor(test.pred[1:nrow(prob.test)],levels=c("1", "2", "3", "4"))

Threshold.Min <- GoldCohort.Threshold # based on the new threshold for Gold cohort (n=101), 8th of October 2015.
maxProbsWhich <- factor(test.pred[1:nrow(prob.test)],levels=c("1", "2", "3", "4", "NC"))

for (i in 1:nrow(prob.test))
{
  if (maxProbs[i] <= Threshold.Min)
    maxProbsWhich[i] <- "NC"      
}
##
###########################################################################################################

maxProbsCol <- ifelse(maxProbsWhich==1,"blue",ifelse(maxProbsWhich==2,"brown1",
                                                     ifelse(maxProbsWhich==3,"khaki",
                                                            ifelse(maxProbsWhich==4,"SeaGreen","SlateGray4")))) # else: "NC"

maxProbsCol2 <- ifelse(maxProbsCol=="blue","#0000FF66", ifelse(maxProbsCol=="brown1","#FF000066",
                                                               ifelse(maxProbsCol=="khaki","#F0E68C", 
                                                                      ifelse(maxProbsCol=="SeaGreen","#2E8B57","#6C7B8B"))))

############################################################
plot.new()
par(mfrow=c(1,1))
par(mar=c(6,4,4,5) + 0.1)
par(cex.axis=0.8)

heading = paste("NanoString subgrouping, ",file.name,"19 genes, classifier ver3.2.1")

boxplot(yaxt="n",xlab="",main=heading,ylab="Probability",probs2[,order(maxProbsCol, maxProbs)],outpch=NA,ylim=c(0,1),las=2, notch=FALSE,
        col=maxProbsCol2[order(maxProbsCol,maxProbs)] )

#maxProbs[order(maxProbsCol,maxProbs)] >=0.75
#abline(col="red",lty=1, v=)
abline(col="grey",lty=1, h=GoldCohort.Threshold)
#abline(col="grey",lty=1, h=0.70)
tmp <- table(maxProbsCol)

grp.sum <- cumsum(c(tmp[1], tmp[2], tmp[3], tmp[4]))

abline(v=c(grp.sum[1] + 0.5, grp.sum[2] +0.5, grp.sum[3]+0.5, grp.sum[4]+0.5))
#lines(col="black",lwd=2,maxProbs[order(maxProbsCol,maxProbs)])
points(col=maxProbsCol[order(maxProbsCol,maxProbs)],pch=19, maxProbs[order(maxProbsCol,maxProbs)])
#points(col=maxProbsCol[order(maxProbsCol,maxProbs)],pch=19, maxProbs[order(maxProbsCol,maxProbs)] >= 0.75)
axis(2, las=2)

legend("bottomleft", legend = c("WNT", "SHH", "Grp3", "Grp4","NC"), col=c("blue", "red", "yellow2", "darkgreen","grey"), pch=19)
axis(2, las=2)

colnames(prob.test) <- c("WNT","SHH","Grp3","Grp4") 
write.csv(prob.test, file = "NanoString_Subgroups_Probabilities.csv") 

#table(as.character(test.pred),as.character(t(RefSub_RNA_Nano)))  # confusion matrix (NanoString vs. NanoStringToronto  call)
# Which samples are discordant?
#hc2 <- as.character(test.pred)
#discor_calls <- which(RefSub_RNA_Nano != hc2) 
#rownames(prob.test)[discor_calls]

#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# End
##################################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################################################
