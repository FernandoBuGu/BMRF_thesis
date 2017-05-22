#From general (filtered byDF), P, EES
setwd("/home/bueno002/From_gitHub/large_coex")
general<-read.table("GO_general.txt",header=T)
general = general[general$evid %in% c('EXP', 'IDA', 'IEP', 'IMP', 'IPI', 'IGI'),]
general = general[general$term == "P",]
general$evid<-NULL
general$term<-NULL
general<-unique(general)
general <- general[c("label", "go")]
write.table(general,file="general_P_EES_ignoreGO50",sep="\t", col.names = F, row.names = F, quote = F)
FILEnet="../large_coex/big_topless.tsv"
FILEann="../large_coex/general_P_EES_ignoreGO50"
FILEclust="../large_coex/big_ipr_labels.tsv"
minGOsize=20
minDFsize=20
maxGOsize=0.8
maxDFsize=0.8
library(ROCR);
library(Matrix);
library(brglm);
library(mvtnorm); 
library(glmnet);
setwd("/home/bueno002/From_gitHub/reference_implementation_yiannis_kourmpetis")
source("bmrf_functions.R");
A = loadSingleNet(FILEnet);
L = loadAnnFixedProteins(FILEann, A);
maxGOsize = (maxGOsize*nrow(A));
Ls = colSums(L);
sel = which(Ls >= minGOsize & Ls <= maxGOsize);
L = L[,sel];
set.seed(888)
GO50<-sample(colnames(L),25)
setwd("/home/bueno002/From_gitHub/large_coex")
write.table(GO50,file="GO50",sep="\t", col.names = F, row.names = F, quote = F)
