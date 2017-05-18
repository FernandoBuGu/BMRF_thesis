setwd("/home/bueno002/From_gitHub/reference_implementation_yiannis_kourmpetis")
FILEnet="../large_coex/big_topless.tsv"		#all coexpression data from yeasy combined
FILEann="../large_coex/big_go.tsv"
FILEclust="../large_coex/big_ipr_labels.tsv"
out="../large_coex/ouT"
minGOsize=20	#as in test from Kourmpetis
minDFsize=20
maxGOsize=0.9
maxDFsize=0.9
library(ROCR);
library(Matrix);
library(brglm);
library(mvtnorm); 
library(glmnet);
source("bmrf_functions.R");

# this is as in bmrf.R
A = loadSingleNet(FILEnet);
L = loadAnnFixedProteins(FILEann, A);
D = loadAnnFixedProteins(FILEclust, A);
minGOsize = minGOsize;
maxGOsize = (maxGOsize*nrow(A));
Ls = colSums(L);
sel = which(Ls >= minGOsize & Ls <= maxGOsize);
L = L[,sel];
minClustsize = minDFsize; 
maxClustsize = (maxDFsize*nrow(A));
Ds = colSums(D);
sel = which(Ds >= minClustsize & Ds <= maxClustsize);
D = D[,sel];
u = rowSums(L);
U = which(u == 0);


#sample 50 go
set.seed(123)
GO50<-sample(colnames(L),50)
L<-L[,colnames(L) %in% GO50]



AUCs<-function(L) {
	df<-as.data.frame(0)
	for(N in 1:50) {
		library(caret)
		a<-L[,N]
		a[U] = -1;	#as in bmrf.R. -1 if the protein has 0 GO terms assigned
		a_check<-a[a!=-1]  #Me: I will not consider these ones for AUC

		#I create 3 folds for the 3 types of labels, so that they are perfectly balanced. consider than only around 20 proteins of 5070 have a "1". Even with this precaution, there are about 5% of the folds for which it gives error: "not enough disticnt predictions to compute AUC". I did not find the reason why this happens sometimes 
		a1s<-a[a==1]
		a0s<-a[a==0]
		am1s<-a[a==-1]
		flds1 <- createFolds(as.factor(a1s), k = 5)					
		flds0 <- createFolds(as.factor(a0s), k = 5)					
		fldsm1 <- createFolds(as.factor(am1s), k = 5)

		AUCs<-c()
		for(i in 1:length(flds1)){
			prots_test<-c()							#for each fold, append to prots_test each element of the fold
			for (j in 1:length(unlist(flds1[i]))) {
				oo<-unlist(flds1[i])[j]
				prots_test<-append(prots_test,a1s[oo])
			}
			for (j in 1:length(unlist(flds0[i]))) {
				oo<-unlist(flds0[i])[j]
				prots_test<-append(prots_test,a0s[oo])
			}
			for (j in 1:length(unlist(fldsm1[i]))) {
				oo<-unlist(fldsm1[i])[j]
				prots_test<-append(prots_test,am1s[oo])
			}
			prots_test<-sample(prots_test)
			name<-names(prots_test)
			a_test<-rep(-1,length(prots_test))	#Proteins in test will enter the model with a -1 (they have to enter the model somehow, otherwise no predictions)
			names(a_test)<-name
			train<-a[!names(a) %in% names(prots_test)]	#train wit all proteins except those in test		
			all<-append(train,a_test) #so, train and test (with -1) enter the model
			Lsingle2 = all
			glmnetpred = try(glmnetDalpha(Y = Lsingle2, X = D, MAXVAR = (ncol(D)-1)), silent=FALSE);			#as in bmrf.R. I had to increase maxit from 1000 to 9000
			if(class(glmnetpred) == 'try-error')
			{
			  next;
			}
			posteriors = try(BMRFz(A=A, L = Lsingle2, D = glmnetpred, burnin = 100, niter = 100), silent=FALSE);	#as in bmrf.R
			posterios_check<-posteriors[names(posteriors) %in% names(a_check)]			#look at those probabilities for which at least one GO is known
			check_probs<-posterios_check[names(posterios_check) %in% names(a_test)]		#look at those probabilities that are in test (and for which at least one GO is known)
			check_labels<-a_check[names(a_check) %in% names(a_test)]					# and the labels
			library(AUC)
			AUC=try(auc(roc(check_probs,as.factor(check_labels))))	#this requires a vector of proteins and a vector of labels
			if(class(AUC) == 'try-error')
			{
			  next;
			}
			AUCs<-append(AUCs,AUC)
		}
		df<-rbind(df,mean(AUCs))	#mean of AUC of all foldes. Get a single AUC per GO (per replicate)
		df<-rbind(df,sd(AUCs))
	}
	return(df)
}


library(pbapply)
w<-pbreplicate(15, AUCs(L)) #call the function 15 times
