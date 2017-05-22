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
#minGOsize = 10
#maxGOsize = (0.8*nrow(A));
FILEann="../large_coex/go_file.uppropagated_extended_ready"
L_m = loadAnnFixedProteins(FILEann, A);
#Ls = colSums(L_m);
#sel = which(Ls >= minGOsize & Ls <= maxGOsize);
#L_m = L_m[,sel];
GO50<-read.table("../large_coex/GO50")
L_m<-L_m[,colnames(L_m) %in% as.character(GO50$V1)]


u = rowSums(L_m);				
U = which(u == 0);


#inside each GO
a<-L_m[,colnames(L_m) == "GO:0030437"]    #"GO:0016192"
a[U] = -1;
a_check<-a[a!=-1] 
a1s<-a[a==1]	
a0s<-a[a==0]
flds1 <- createFolds(as.factor(a1s), k = 15)							
flds0 <- createFolds(as.factor(a0s), k = 15)					
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
prots_test <- sample(prots_test)
name<-names(prots_test)
a_test<-rep(-1,length(prots_test))	
names(a_test)<-name
train<-a[!names(a) %in% names(prots_test)]	
Lsingle2<-append(train,a_test)
Lsingle2<-Lsingle2[order(names(Lsingle2))]


glmnetpred_m = try(glmnetDalpha(Y = Lsingle2, X = D_m, MAXVAR = (ncol(D_m)-1)), silent=FALSE);
posteriors = try(BMRFz(A, Lsingle2, glmnetpred_m, burnin = 4000, niter = 4000), silent=FALSE);	

posteriors_ordered<-posteriors[order(-posteriors)]		#This shows the highest predictionsç
posteriors_ordered[1:20]

posterios_check<-posteriors[names(posteriors) %in% names(a_test)]	
posterios_check_ordered<-posterios_check[order(-posterios_check)]	
posterios_check_ordered[1:10]	

check_labels<-a_check[names(a_check) %in% names(a_test)]					
library(AUC)
AUC_including_train=try(auc(roc(posteriors,as.factor(a))))	#this gives 82%
roc_including_train=roc(posteriors,as.factor(a))
AUC_only_test=try(auc(roc(posterios_check,as.factor(check_labels))))
roc_only_test=roc(posterios_check,as.factor(check_labels))


library(ROCR)
pred <- prediction( posterios_check,as.factor(check_labels) )



#go from  -1    0    1 
#		2367 3075  318
to
		3500  1700 318


Reason 2) in previous email is the one
with the following proportions in the train
  -1    0    1 
3766 1614  380

I achieve AUC:75 for GO term x with some True Positives, for  a test set that was set to "-1"
So, by reducing GO file, we can achieve good reasults. The only problem is that with so many -1s, a large number of gibbs iterations are reuired (~3000, which takes ~10 minutes per fold and per GO. So, we will have to reduce #folds as much as possible).

knowing this I will try to fill the tables again. This time including results of FP,FN,TP,TN
