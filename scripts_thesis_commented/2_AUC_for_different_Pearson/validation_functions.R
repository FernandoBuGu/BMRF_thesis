#Script with the functions that are requred to run "runAUC.R" in the same folder.

    #MSc thesis bioinfomratics WUR. Protein function prediction for poorly annotated species.
    #Author: Fernando Bueno Gutierrez
    #email1: fernando.buenogutierrez@wur.nl
    #email2: fernando.bueno.gutie@gmail.com


AUCf<-function(go,k,all,noit_GS=30,PR=0,RMO=0,domainMAX=0.9,domainMIN=20)
{
#Computes AUC for a given GO term given some parameters.

#In:
    #go: int, index of the GO term in the only-GO file
    #k: int, number of folds
    #all:logical, true if the non-valid associations are also considered in the trainning set.
    #noit_GS: nº iterations in the Gibb sampling in the Yiannis code
    #PR (percentage remove), numeric from 0 to 1. percentage of rows to remove randomly from the network file
    #RMO: numeric from 0 to 1. percentage of rows to remove randomly from the GO-file

#Out:
    #mean(AUCs,na.rm = TRUE). The mean AUC of the k folds.


    #print some information
	print("myGOLm")
	print(dim(myGO$L_m)) #dimension of the GO file considering only the rows of the target GO term 
	print("myGOA")
	print(dim(myGO$A)) #dimension of the network file
	print("myGOminus_ones")
	print(length(myGO$minus_ones)) #number of "minus_ones". Genes associated with 0 GO terms
	print("GOfile") 
	print(dim(myGO$GOfile)) # dimension of the GO file considering all GO terms 
	print("dimD_m")
	print(dim(myGO$D_m)) #dimension of the domain file




	library(Matrix);
	myGO<-reduce_myGO(go,reduce,PR,RMO,domainMAX,domainMIN)

	go_number<-go	
	GOfile=myGO$GOfile
	names<-GOfile[as.character(GOfile$go) == as.character(go_number),] # names (df): associations of the target GO term
	valid_positive_GO<-names[names$reliab == "valid",] 
	print("dim1_valid_positive_GO:#positive cases valid")
	print(dim(valid_positive_GO)[1])



	data_fold<-create_Lsingle2(k,all=all,myGO)
	AUCs<-c()
	for(i in 1:length(data_fold$dat)){
		glmnetpred = try(glmnetDalpha(Y = unlist(data_fold$dat[i]), X = myGO$D_m, MAXVAR = (ncol(myGO$D_m)-1)), silent=FALSE)
		if(class(glmnetpred) == 'try-error')
		{
		  return(NA);
		}
		posteriors = try(BMRFz(myGO$A, unlist(data_fold$dat[i]), glmnetpred, burnin = noit_GS, niter = noit_GS), silent=FALSE)
		if(class(posteriors) == 'try-error')
		{
		  return(NA);
		}
		posteriors = round(calibrate(posteriors),5)
		posteriors_test<-posteriors[names(posteriors) %in% names(unlist(data_fold$test[i]))]	
		posteriors_test<-posteriors_test[order(names(posteriors_test))]
		labels_test<-unlist(data_fold$test[i])
		labels_test<-labels_test[order(names(labels_test))]
		library(AUC)	
		AUCc=try(auc(roc(posteriors_test,as.factor(labels_test))));
		if(class(AUCc) == 'try-error')
		{
		  return(NA);
		}
		AUCs<-append(AUCs,AUCc)
	}
	#
	l<-list("dim_myGO$L_m"=dim(myGO$L_m),"sum_myGO$L_m"=sum(myGO$L_m),"AUC"=mean(AUCs,na.rm = TRUE))
	return(mean(AUCs,na.rm = TRUE))
}


create_Lsingle2<- function(k=10,all,myGO)
{

# given noGP, and k-cross-validation, it returns Lsingle2

#In: 
    # k: integer that refers to the number of folds in the cross validation
    # all: logical that if true transforms an integer referring to the (index) of the go-term in L_m into a go_number (or "GO ID")
    # myGO: list with the data frames that are required to run BMRF (A, L_m and D_m) after removing PR % of the rows

#Out:
    # lista: list with the objects "test", "dat" and "labels".

	go_number<-go
	GOfile=myGO$GOfile

	labels<-myGO$L_m[,colnames(myGO$L_m) == as.character(go_number)] #labels: a binary table with: (1) names of all genes; (2) value "1" if assoc. with target GO
	print("sum_labels:#positive cases")
	print(sum(labels))
	

	names<-GOfile[as.character(GOfile$go) == as.character(go_number),] # names (df): associations of the target GO term
	valid_positive_GO<-names[names$reliab == "valid",]  #valid_positive_GO: a df(Zx3) with all the valid associations for the target GO
	NONvalid_positive_GO<-names[names$reliab == "NONvalid",]

	labels_nega<-labels[!names(labels) %in% valid_positive_GO$label & !names(labels) %in% NONvalid_positive_GO$label] #labels_nega: a binary table with: (1) genes without assoc. with target GO; (2) value "0" in all cases

	labels[names(labels) %in% myGO$minus_ones] = -1; # myGO$minus_ones. a list with the names of all genes that are not associated with any GO term and that have not been removed from the network 
	print("minusOnes")
	print(length(myGO$minus_ones)) 
	labels_check<-labels[labels!=-1] #Check among labels that do not have minusones. (although not all of them, as seen in next 2 lines)
	labels_check_posi<-labels_check[names(labels_check) %in% valid_positive_GO$label] #labels_check_posi: a binary table with: (1) genes assoc. with target GO and reliab="valid"; (2) value "1" in all cases
	labels_check_nega<-labels_check[names(labels_check) %in% names(labels_nega)] #labels_check_nega: a binary table with: (1) genes WITHOUT assoc. with target GO that are not "minusones"; (2) value "0" in all cases

	folds_1s<-folds(labels_check_posi,k)
	folds_0s<-folds(labels_check_nega,k) 
	Lsingle2_L<-list()
	prots_test_L<-list()
	for(i in 1:k){
		prots_test<-c()
		prots_test<-append(prots_test,unlist(folds_1s[i]))
		prots_test<-append(prots_test,unlist(folds_0s[i]))
		name<-names(prots_test)
		test<-rep(-1,length(prots_test))
		names(test)<-name
		prots_test_L<-append(prots_test_L,list(prots_test))
		train<-labels[!names(labels) %in% names(test)]	
		Lsingle2<-append(train,test)
		Lsingle2<-Lsingle2[order(names(Lsingle2))]
		Lsingle2_L<-append(Lsingle2_L,list(Lsingle2))
	}
	labels<-myGO$L_m[,colnames(myGO$L_m) == as.character(go_number)]
	lista<-list("test"=prots_test_L,"dat"=Lsingle2_L,"labels"=labels)
	return(lista)
}


folds<-function(v,k=10){
# Given a vector of 1s or 0s, it returns k folds
	all<-c()
	list_test<-list()
	for(i in 1:k) {
		samplable<-v[!names(v) %in% all]
		take<-floor(1/k*length(v))
		if(i == 1){take_also<-(length(v)-take*k)}
		if(take_also != 0){
			s<-sample(samplable,take+1)
			take_also<-take_also-1
		} else {
			s<-sample(samplable,take)
		}
		list_test<-append(list_test,list(s))
		all<-append(all,names(s))
		i=i+1
	}
	return(list_test)
}



return_L_m<-function(subset=F,minGOsize=20,maxGOsize=0.1,only_EES_BP=F){

# Returns a data_frame L_m with GOs as columns and Genes as rows, where all GOs have passed the filters and are in subset_network.

# subset: an integer that refers to the folder of network data from which network data will be extracted. All folders of networks are in 
# 		/home/WUR/bueno002/humans/run/run_subsets/network_subsets. If false, all network data will be used.
# minGOsize and maxGOsize: Criterias such as in https://github.com/jwbargsten/bmrf
# only_EES_BP: a logical that if true, only associations with Experimental evidence scores and Biological process will be used to train the model
#
# L_m: a data_frame with GOs as columns and Genes as rows, where all GOs have passed the filters and are in subset_network.
# GOfile: a data_frame corresponding to the GO file, where all GOs have passed the filters and are in subset_network.d 
# network: a data_frame corresponding to the network file as loaded

	GOfile<-read.table("/home/bueno002/pearsons_ready/rbind_upvalid_upNONvalid_02" ,header=F)
	network<-read.table("/home/bueno002/pearsons_ready/big_topless02_updated.tsv" ,header=F,fill=T)
    GOfile<-GOfile[GOfile$V1 %in% network$V1 | GOfile$V1 %in% network$V2,]
    #network<-network[complete.cases(network), ]
#    N<-network[nchar(as.character(network[, "V1"])) >2,] Not required because big_topless02.tsv has been updated
#    N2<-N[nchar(as.character(N[, "V2"])) >2,]
    
#    network<-network[network$V1 %in% GOfile$V1,] 
#    network<-network[network$V2 %in% GOfile$V1,] 
	A = loadSingleNet(network)


	#options(warn=0)
	#subset network
	#I skip part of "if subset" for chickens


	#Go file: apply filter only on valid
	valid_positive<-GOfile[GOfile$V3 == "valid",]  #valid_positive: the associations in GOfile that can be used in the validation(EES and BP)
	NONvalid_positive<-GOfile[GOfile$V3 == "NONvalid",]
	valid_positive$V1<-NULL
	valid_positive$V3<-NULL
	table<-table(valid_positive$V2) #table: a table that counts the number of asssoc (usable for validation) per GO
	minGOsize=minGOsize
	maxGOsize=maxGOsize*nrow(A)
	table<-table[table >= minGOsize & table <= maxGOsize]
	valid_positive<-GOfile[GOfile$V3 == "valid",]	
	valid_positive<-valid_positive[valid_positive$V2 %in% names(table),]

	#unless "only_EES_BP", add the info from NONvalid for those GOs that passed the filter 
	NONvalid_positive<-NONvalid_positive[NONvalid_positive$V2 %in% valid_positive$V2,]
	GOfile<-rbind(valid_positive,NONvalid_positive)
	if(only_EES_BP){	
		GOfile<-valid_positive
	}
	    prots<-as.vector(unique(network$V1))
	DF<-append(prots,as.vector(unique(network$V2)))
    GOfile<-GOfile[GOfile$V1 %in% DF,]
	L_m = loadAnnFixedProteins(GOfile, A); #i.e. for minGOsize:50, maxGOsize:0.9 and only_EES_BP, L_m:9998x1093
 #i.e. for minGOsize:50, maxGOsize:0.9 and only_EES_BP, L_m:9998x1093


	Ret<-list("L_m" = L_m, "GOfile"=GOfile, "network"=network)
	return(Ret)
}



reduce_myGO<-function(go,reduce=F,PR=0,RMO=0,domainMAX,domainMIN){

# Returns a list with the data frames that are required to run BMRF (A, L_m and D_m) after removing PR % of the rows. Also returns 

# go: an integer that refers to the index of the GO file in L_m (from the "return_L_m" function) for which we want to measure prediction 		#			performance. The indexes range from 1 to dim(loaded$L_m)[2]
# reduce: a character that specifies what type of reduction we want to do: 
# 			amg,oa,epmg,oe   associationsof my go; other associations; edges of proteins of my go; other edges	
# PR: a decimal (from 0 to 1) that species the portion of rows to be removed (either from the GOfile or from the network file 
#			(depending on "reduce")


	GOfile<-loaded$GOfile
	colnames(GOfile)<-c("label","go","reliab")
	network<-loaded$network
	L_m<-loaded$L_m


	go_number<-go	

	proteins_myGO<-GOfile[GOfile$go == go_number,]
	proteins_myGO_V<-proteins_myGO$label[proteins_myGO$reliab=="valid"]
	proteins_NOTmyGO<-unique(GOfile$label[!GOfile$label %in% proteins_myGO])

	network_proteins_myGO_pp<-network[network$V1 %in% proteins_myGO_V & network$V2 %in% proteins_myGO_V,]
	network_proteins_myGO_pn1<-network[network$V1 %in% proteins_myGO_V & !network$V2 %in% proteins_myGO_V,]
	network_proteins_myGO_pn2<-network[!network$V1 %in% proteins_myGO_V & network$V2 %in% proteins_myGO_V,]
	network_proteins_myGO_pn<-sample(rbind(network_proteins_myGO_pn1,network_proteins_myGO_pn2))
	non_nn<-rownames(network_proteins_myGO_pp)
	non_nn<-append(non_nn,rownames(network_proteins_myGO_pn))
	network_proteins_myGO_nn<-network[!rownames(network) %in% non_nn,]

	assoc_myGO<-GOfile[GOfile$go == go_number,]											#The associations of the target GO term
	assoc_myGO_Val<-assoc_myGO[assoc_myGO$reliab == "valid",]							#The validated associations of the target GO term
	other_assoc<-GOfile[GOfile$go != go_number,]										#associations that are not form the target GO term
		#If other_assoc becomes very small, n decreases, it starts to fail, cause of lack of ceros in validation

	#reduce whatever
	if(reduce=="amg"){
		out_GOfile<-assoc_myGO_Val[sample(nrow(assoc_myGO_Val), dim(assoc_myGO_Val)[1]*PR), ]
		GOfile<-GOfile[!rownames(GOfile) %in% rownames(out_GOfile),] 
	}
	if(reduce=="oa"){
		out_GOfile<-other_assoc[sample(nrow(other_assoc), dim(other_assoc)[1]*PR), ]
		GOfile<-GOfile[!rownames(GOfile) %in% rownames(out_GOfile),] 
	}
	if(reduce=="epp"){
		out_network<-network_proteins_myGO_pp[sample(nrow(network_proteins_myGO_pp), dim(network_proteins_myGO_pp)[1]*PR), ]
		network<-network[!rownames(network) %in% rownames(out_network),]
	}
	if(reduce=="epn"){
		out_network<-network_proteins_myGO_pn[sample(nrow(network_proteins_myGO_pn), dim(network_proteins_myGO_pn)[1]*PR), ]
		network<-network[!rownames(network) %in% rownames(out_network),]
	}
	if(reduce=="enn"){
		out_network<-network_proteins_myGO_nn[sample(nrow(network_proteins_myGO_nn), dim(network_proteins_myGO_nn)[1]*PR), ]
		network<-network[!rownames(network) %in% rownames(out_network),] 
	}


	#restrict domains file to network.
	A = loadSingleNet(network);
	prots<-as.vector(unique(network$V1))
	DF<-append(prots,as.vector(unique(network$V2)))
	DF<-unique(DF)  #DF: a list that contains all genes in the network. It may have been modified with respect to previous
	print("lengthDF")
	print(length(DF))
	domains<-read.table("/home/bueno002/pearsons_ready/big_ipr_labels.tsv" ,header=F)
	domains<-domains[domains$V1 %in% DF,] 

	#filter domains (as always)
	minDFsize=20
	maxDFsize=0.9
	minClustsize = minDFsize; 
	maxClustsize = (maxDFsize*nrow(A));
	D_m = loadAnnFixedProteins(domains, A);
	Ds = colSums(D_m);
	sel = which(Ds >= minClustsize & Ds <= maxClustsize);
	D_m = D_m[,sel];


	#restrict GO file to network.
	GOfile<-GOfile[GOfile$label %in% DF,] #The GO file after removing rows either in network or in GO file
	print("GOfile")
	print(dim(GOfile))
	GOfile<-unique(GOfile)
	print(dim(GOfile))
	L_m = loadAnnFixedProteins(GOfile, A)


	proteins_per_go = rowSums(L_m);	   #rowSums: give one value per row (sum columns). proteins_per_go:list counting GOs per gene		
	minus_ones = names(proteins_per_go)[proteins_per_go == 0] 
	#minus_ones: list of the names of the genes that will be "minus_ones" (0 GO terms associated with that gene)
	#i.e. for minGOsize:50, maxGOsize:0.9, only_EES_BP, PR=0, length(minus_ones):4416
	#i.e. for minGOsize:1, maxGOsize:1, 	only_EES_BP, PR=0, length(minus_ones):1554


	#The number of minus ones is very large (1554) with no filters, and many more if I apply filters and/or I remove rows. 
	#It may be thus a good idea to remove some of these genes. (i.e:0.8) RMO:remove_minus_ones
	genes_remove<-sample(minus_ones, length(minus_ones)*RMO)
	minus_ones<-minus_ones[!minus_ones %in% genes_remove]
	D_m<-D_m[!rownames(D_m) %in% genes_remove,]
	A<-A[!rownames(A) %in% genes_remove,]
	A<-A[,!colnames(A) %in% genes_remove]
	#GOfile does not need to be updated because "genes_remove" are not in GOfile
	L_m<-L_m[!rownames(L_m) %in% genes_remove,]

	Ret<-list("L_m" = L_m, "D_m" = D_m, "A" = A, "minus_ones" = minus_ones, "GOfile"=GOfile)

	return(Ret)
}

