AUCf<-function(i,k,all)
{
	library(Matrix);
	data_fold<-create_Lsingle2(k,all=all)
	glmnetpred = try(glmnetDalpha(Y = unlist(data_fold$dat[i]), X = loaded$D_m, MAXVAR = (ncol(loaded$D_m)-1)), silent=FALSE)
	if(class(glmnetpred) == 'try-error')
	{
	  return(NA);
	}
	posteriors = try(BMRFz(loaded$A, unlist(data_fold$dat[i]), glmnetpred, burnin = noit_GS, niter = noit_GS), silent=FALSE)
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
	return(AUCc)
}



create_Lsingle2<- function(k=10,all)
{
"given noGP, and k-cross-validation, it returns Lsingle2"
"iGP: the index of GP in the L_m GO file"
"Lsingle2: a list of labels for that GO term, with k folds. In each fold, k% of the proteins have been set to 0"	
	if(all){
		go_number<-as.character(colnames(loaded$L_m)[go])	
	}
	labels<-loaded$L_m[,colnames(loaded$L_m) == as.character(go_number)]
	valid_positive_GO<-loaded$valid_positive[loaded$valid_positive$V2 == as.character(go_number),]
	NONvalid_positive_GO<-loaded$NONvalid_positive[loaded$NONvalid_positive$V2 == as.character(go_number),]
	labels_nega<-labels[!names(labels) %in% valid_positive_GO$V1 & !names(labels) %in% NONvalid_positive_GO$V1]
	labels[loaded$U_after] = -1;
	labels_check<-labels[labels!=-1] 
	labels_check_posi<-labels_check[names(labels_check) %in% valid_positive_GO$V1]
	labels_check_nega<-labels_check[names(labels_check) %in% names(labels_nega)]
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
	labels<-loaded$L_m[,colnames(loaded$L_m) == as.character(go_number)]
	lista<-list("test"=prots_test_L,"dat"=Lsingle2_L,"labels"=labels)
	return(lista)
}


folds<-function(v,k=10){
"Given a vector of 1s or 0s, it returns k folds"
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



load_L_m<-function(P,subset=F,minGOsize=20,maxGOsize=0.1,only_EES_BP=F){
	
	#subset network
	if(subset){	
		path<-paste("/home/bueno002/tables_monday/data/fileS/subset",subset,sep="")
		filenames <- list.files(path=path, full.names=TRUE)
		if(length(filenames) == 1){
			data<-read.table(filenames)
		} else {
			import.list <- llply(filenames, read.table)
			library(reshape)
			data <- merge_recurse(import.list)
		}
		data[3]<-NULL
		dat.sort = t(apply(data, 1, sort))
		data<-data[!duplicated(dat.sort),]
		data$V3<-1
		network<-data
	}

	#reduce network size randomly
	if(subset==F){
		network<-read.table("/home/bueno002/tables_monday/data/big_topless.tsv" ,header=F)
		out_network<-network[sample(nrow(network), dim(network)[1]*P), ]
		network<-network[!rownames(network) %in% rownames(out_network),]
	}

	#restrict domains file to network
	A = loadSingleNet(network);
	prots<-as.vector(unique(network$V1))
	DF<-append(prots,as.vector(unique(network$V2)))
	DF<-unique(DF)
	domains<-read.table("/home/bueno002/tables_monday/data/big_ipr_labels.tsv" ,header=F)
	domains<-domains[domains$V1 %in% DF,] 
	D_m = loadAnnFixedProteins(domains, A);

	#filter domains (as always)
	minDFsize=10
	maxDFsize=0.9
	minClustsize = minDFsize; 
	maxClustsize = (maxDFsize*nrow(A));
	Ds = colSums(D_m);
	sel = which(Ds >= minClustsize & Ds <= maxClustsize);
	D_m = D_m[,sel];

	#restrict go file to network
	GOfile<-read.table("/home/bueno002/tables_monday/data/rbind_upvalid_upNONvalid" ,header=F)
	rb_upvalid_upNONvalid<-GOfile[GOfile$V1 %in% DF,]

	#Go file: apply filter only on valid
	valid_positive<-rb_upvalid_upNONvalid[rb_upvalid_upNONvalid$V3 == "valid",]
	NONvalid_positive<-rb_upvalid_upNONvalid[rb_upvalid_upNONvalid$V3 == "NONvalid",]
	valid_positive$V1<-NULL
	valid_positive$V3<-NULL
	table<-table(valid_positive$V2)
	minGOsize=minGOsize
	maxGOsize=maxGOsize*nrow(A)
	table<-table[table >= minGOsize & table <= maxGOsize]
	valid_positive<-rb_upvalid_upNONvalid[rb_upvalid_upNONvalid$V3 == "valid",]	
	valid_positive<-valid_positive[valid_positive$V2 %in% names(table),]

	#unless "only_EES_BP", add the info from NONvalid for those GOs that passed the filter 
	NONvalid_positive<-NONvalid_positive[NONvalid_positive$V2 %in% valid_positive$V2,]
	rb_upvalid_upNONvalid<-rbind(valid_positive,NONvalid_positive)
	if(only_EES_BP){	
		rb_upvalid_upNONvalid<-valid_positive
	}
	L_m = loadAnnFixedProteins(rb_upvalid_upNONvalid, A);

	#provide some info
	o<-sort(-table)				#This contain #valid proteins  
	u_after = rowSums(L_m);				
	U_after = which(u_after == 0);	#This is the # of "-1s"
	o_total<-table(rb_upvalid_upNONvalid$V2)		#This contain #total proteins 


	#count # connexions
	


	l<-list("L_m" = L_m, "valid_positive" = valid_positive, "NONvalid_positive" = NONvalid_positive, "#GOs" = dim(L_m)[2], "#-1s" = length(U_after), "table_EES_BP"=o_total, "table"=o, "U_after" = U_after, "D_m" = D_m, "A" = A, "network_size" = dim(network)[1], "network" = network)
	return(l)
}





