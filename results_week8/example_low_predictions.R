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
a<-L_m[,colnames(L_m) == "GO:0006888"]
a[U] = -1;
a_check<-a[a!=-1] 
a_check2<-a[a==0]
a_check<-append(a_check,a_check2)
a1s<-a[a==1]
a1s<-a1s[names(a1s) %in% names(a_check)]	
a0s<-a[a==0]
am1s<-a[a==-1]
flds1 <- createFolds(as.factor(a1s), k = 5)							
flds0 <- createFolds(as.factor(a0s), k = 5)					
fldsm1 <- createFolds(as.factor(am1s), k = 5)
AUCs<-c()
prots_test<-c()					
for (j in 1:length(unlist(flds1[1]))) {
	oo<-unlist(flds1[1])[j]
	prots_test<-append(prots_test,a1s[oo])
}
for (j in 1:length(unlist(flds0[1]))) {
	oo2<-unlist(flds0[1])[j]
	prots_test<-append(prots_test,a0s[oo2])
}
for (j in 1:length(unlist(fldsm1[1]))) {
	oo3<-unlist(fldsm1[1])[j]
	prots_test<-append(prots_test,am1s[oo3])
}
prots_test <- sample(prots_test)
name<-names(prots_test)
a_test<-rep(1,length(prots_test))	
names(a_test)<-name
train<-a[!names(a) %in% names(prots_test)]	
Lsingle2<-append(train,prots_test)
Lsingle2<-Lsingle2[order(names(Lsingle2))]

glmnetpred_m = try(glmnetDalpha(Y = Lsingle2, X = D_m, MAXVAR = (ncol(D_m)-1)), silent=FALSE);
posteriors = try(BMRFz(A, Lsingle2, glmnetpred_m, burnin = 20, niter = 20), silent=FALSE);	
posteriors_ordered<-posteriors[order(-posteriors)]		#This shows the highest predictions
posterios_check<-posteriors[names(posteriors) %in% names(a_check)]			
check_probs<-posterios_check[names(posterios_check) %in% names(prots_test)]		
check_labels<-a_check[names(a_check) %in% names(a_test)]					
library(AUC)
AUC_including_train=try(auc(roc(posteriors,as.factor(a))))	#this gives 82%
AUC_only_test=try(auc(roc(check_probs,as.factor(check_labels))))		#this below 0.5
AUC_including_train
posteriors_ordered
AUC_only_test
