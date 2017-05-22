setwd("/home/bueno002/From_gitHub/reference_implementation_yiannis_kourmpetis")
library(ROCR);
library(Matrix);
library(brglm);
library(mvtnorm); 
library(glmnet);
library(caret);
source("bmrf_functions.R");


#network and domains files
FILEnet="../large_coex/big_topless.tsv"		
FILEclust="../large_coex/big_ipr_labels.tsv"
A = loadSingleNet(FILEnet);
D_m = loadAnnFixedProteins(FILEclust, A);
minDFsize=30
maxDFsize=0.8
minClustsize = minDFsize; 
maxClustsize = (maxDFsize*nrow(A));
Ds = colSums(D_m);
sel = which(Ds >= minClustsize & Ds <= maxClustsize);
D_m = D_m[,sel];


#GO file
minGOsize = 10
maxGOsize = (0.8*nrow(A));
FILEann="../large_coex/go_file.uppropagated_extended_ready"
L_m = loadAnnFixedProteins(FILEann, A);
Ls = colSums(L_m);
sel = which(Ls >= minGOsize & Ls <= maxGOsize);
L_m = L_m[,sel];



u = rowSums(L_m);				
U = which(u == 0);


#inside each GO
Lsingle2<-L_m[,colnames(L_m) == "GO:0006888"]
Lsingle2[U]<--1
glmnetpred_m = try(glmnetDalpha(Y = Lsingle2, X = D_m, MAXVAR = (ncol(D_m)-1)), silent=FALSE);
posteriors = try(BMRFz(A, Lsingle2, glmnetpred_m, burnin = 20, niter = 20), silent=FALSE);	
posteriors_ordered<-posteriors[order(-posteriors)]		#This shows the highest predictions. Only 4 above 0.5			
library(AUC)
AUC_including_train=try(auc(roc(posteriors,as.factor(a))))	#this gives 82%
AUC_including_train
posteriors_ordered

