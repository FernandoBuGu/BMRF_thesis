#Script with the functions that are requred to compute Area under the curve (AUC) using ppi (protein-protein-interactions) data. 

    #MSc thesis bioinfomratics WUR. Protein function prediction for poorly annotated species.
    #Author: Fernando Bueno Gutierrez
    #email1: fernando.buenogutierrez@wur.nl
    #email2: fernando.bueno.gutie@gmail.com


AUC_func = function(minGOsize=aa, minDFsize=bb, maxGOsize=cc, maxDFsize=dd, only_biological=ee, only_EES=ff, portion_test=gg, note=hh)
{
#function to return AUC for a given GO term. The function needs to be run for each GO term.

#In:
    #minGOsize: int, minimum nº of genes that the GO terms should be associated with
    #minDFsize: int, minimum nº of genes that the domains should be associated with
    #maxGOsize: int, maximum nº of genes that the GO terms should be associated with
    #maxDFsize: int, maximum nº of genes that the domains should be associated with
    #only_biological: logical, whether only the biologial process category should be considered
    #only_EES: logical, whether only Experimental eviden scores (EES) should be considered
    #portion_test: numeric form 0 to 1. Which protion of the genes should be in the test set
    #note: character, to help in visualizing the results.

#Out:
    #Results: a list with several results. "mean_AUCs" is the main resuts and refers to the AUC of the GO term considered. 
	
	setwd("/home/bueno002/From_gitHub/ppi_data")

	GO = read.table("GO_general.txt",header=T)
	GO = GO[GO$term == 'P',]
	GO = GO[GO$evid %in% c('EXP', 'IDA', 'IEP', 'IMP', 'IPI', 'IGI'),]
	Go_complete <- unique(GO)
	Go_complete$label<-NULL
	Go_complete<-Go_complete[!(duplicated(Go_complete[c("go")])),]

	ppiGO<-read.table("GO_file.tsv",header=F)
	colnames(ppiGO)[1]<-"go"
	Go_complete<-merge(Go_complete,ppiGO,by="go",all.y=T)
	colnames(Go_complete)[4]<-"label"
	Go_complete <- Go_complete[c("go", "label", "evid", "term")]
	#out what is not in network
	network<-read.table("collins.dat")
	labels<-c(as.character(network$V1),as.character(network$V2))
	Go_complete<-Go_complete[Go_complete$label %in% labels,]
	if(only_biological) {
		Go_complete = Go_complete[Go_complete$term == 'P',]
	}
	if(only_EES) {
		Go_complete = Go_complete[Go_complete$evid %in% c('EXP', 'IDA', 'IEP', 'IMP', 'IPI', 'IGI'),]
	}

	Go_complete_2c <- subset(Go_complete[,1:2])
	Go_complete_2c <- unique(Go_complete_2c)
	Go_complete_3c<- subset(Go_complete[,1:3])


	#Merge domain info 
	domains<-read.table("big_ipr_labels.tsv")
	colnames(domains)[1]<-"label"
	domains_ppi<-merge(Go_complete,domains,by="label",all=F)
	domains_ppi$go<-NULL
	domains_ppi$evid<-NULL
	domains_ppi$term<-NULL
	domains_ppi<-unique(domains_ppi)
	write.table(domains_ppi,"ppi_ipr_labels.tsv", sep="\t", col.names = F, row.names = F, quote = F)




	
	#Separate test and train. GO terms in test set need to be also in train set, otherwise no predictions. So, I split the GO_label combinations of some of the GO labels between test and train sets.
	df_Freq<-as.data.frame(table(Go_complete_2c$go)) 
	df_Freq$removable<-df_Freq$Freq-minGOsize	#For each GO term we can "remove" (put in the test set) "df_Freq$removable" labels
	candidates<-df_Freq[df_Freq$removable>0,]
	candidate_rows <- Go_complete_3c	
	candidate_rows	<- candidate_rows[candidate_rows$go %in% candidates$Var1,] #Candidate_rows are the GO_label comb. for which some labels can go to test set 

	amount<-round(length(unique(candidates$Var1))*portion_test,0)#Amount of GO terms in test set.
	sampled_go<-sample(candidates$Var1,amount)
	extract <- candidates[candidates$Var1 %in% sampled_go,]
	extractable_rows<- candidate_rows[candidate_rows$go %in% extract$Var1,]	#Candidate_rows are the GO_label comb. of those GO terms that were chosen of test set (sampling). (extractable_rows are some of the candidate_rows)

	rows_test = data.frame()
	for(i in 1:length(extract$Var1)) #For each GO term to be in test and in train set. 
	{
		extractable_rows_af <- extractable_rows[extractable_rows$evid %in% c('EXP', 'IDA', 'IEP', 'IMP', 'IPI', 'IGI'),]#all comb. in test need to be experim evid scores.

#The following are to place in the test set all combinations except the minimum that are required in the train set. "number_rows_out" are the no. of combs that go to test set.
		extractable_real<-length(which( extractable_rows_af$go == extract$Var1[i] ))
		extractable_total<-length(which( extractable_rows$go == extract$Var1[i] ))
		if(extractable_real==extractable_total) {	#In some cases, all comb. have evidence scores and number rows (# comb) is simply the difference
					number_rows_out<-extractable_real-minGOsize
		}	else if(minGOsize-(extractable_total-extractable_real)<0) {#Some cases, the combs. without evidence scores are sufficient to cover that minimum in the train set, so evert coomb with validation goes to test set.
			number_rows_out<-extractable_real
		}	else { #Otherwise, some comb. with evidence scores need to be in train.
			number_rows_out<-extractable_real-(minGOsize-(extractable_total-extractable_real))
		}
		t<-extractable_rows_af[ sample( which( extractable_rows_af$go == extract$Var1[i] ) , number_rows_out) , ] #Of those with evidence scores, we sample number_rows_out combs.
		df <- data.frame(row.names(t))
		rows_test <- rbind(rows_test,df)
	}
	GO_test <- extractable_rows[rownames(extractable_rows) %in% rows_test$row.names.t.,]
	GO_test <- subset(GO_test[,1:2])
	GO_test <- GO_test[c("label", "go")]
	dim(GO_test)
	x <- rbind(Go_complete_2c, GO_test) #In train set, all other combs. 
	GO_train<-x[! duplicated(x, fromLast=TRUE) & seq(nrow(x)) <= nrow(Go_complete_2c), ]
	GO_train <- GO_train[c("label", "go")]

	name_file<-paste(note, "GO_train", ".tsv",sep="_") #"note" to avoid overwriting
	name_file_out<-paste(note, "Out", ".out",sep="_")
	file.remove(name_file)
	setwd("/home/bueno002/From_gitHub/reference_implementation_yiannis_kourmpetis")
	write.table(GO_train,file=name_file,sep="\t", col.names = F, row.names = F, quote = F)
	file.remove(name_file_out)
	write.table("",file=name_file_out,sep="\t", col.names = F, row.names = F, quote = F)

	



	FILEnet="../ppi_data/collins.dat"
	FILEann=name_file
	FILEclust="../ppi_data/ppi_ipr_labels.tsv"
	out=name_file_out
	source("bmrf_modified2.R"); 
	bmrf(FILEnet, FILEann, FILEclust, out, minGOsize,minDFsize, maxGOsize, maxDFsize)


	test<-GO_test
	train<-GO_train
	pred<-read.table(file=name_file_out,fill=T)
	file.remove(name_file_out)
	colnames(pred)<-c("labels","go","p")
	colnames(test)<-c("labels","go")
	colnames(train)<-c("labels","go")
	pred_2c<-subset(pred[,1:2])
	x <- rbind(pred_2c, train)
	to_see<-x[! duplicated(x, fromLast=TRUE) & seq(nrow(x)) <= nrow(pred_2c), ]#we look at the combs. that were predicetd, are not in train and are in test.
	to_see_merge<-merge(pred,to_see,by=c("labels","go"))
	to_see_merge<-to_see_merge[to_see_merge$go %in% test$go,]
	library(AUC)
	vecs<-c()
	to_see_GOxx_all<-c()
	AUCs<-c()
	failures<-c()
	for(i in 1:length(unique(to_see_merge$go))) #For each comb that we are going to evaluate
	{
		to_see_GOxx<-to_see_merge$p[to_see_merge$go == as.character(unique(to_see_merge$go)[i])] #We look at the probab in output for that comb. 
		to_see_GOxx_all<-append(to_see_GOxx_all,to_see_GOxx) #and append them to a vector of probabilities

		pred_labels<-pred$labels[pred$go==as.character(unique(to_see_merge$go)[i])]
		train_labels<-train$labels[train$go==as.character(unique(to_see_merge$go)[i])]
		pred_labels<-pred_labels[!pred_labels %in% train_labels]
		test_label<-test$labels[test$go==as.character(unique(to_see_merge$go)[i])]
		v <- c()
		for(j in 1:length(pred_labels))#If we find that a label (not present in train for that GO) is assigned to the GO term in the test, we give value "1" for that comb.  
		{
			if(pred_labels[j] %in% test_label) {
				v<-append(v,1)
			} else {
				v<-append(v,0)
			}
		}
		vecs<-append(vecs,v)
		AUC=try(auc(roc(to_see_GOxx,as.factor(v)))); #Roc requires a vector of probabilities and a binary vector. 
		if(class(AUC) == 'try-error')
		{
		  failures<-append(failures,as.character(unique(to_see_merge$go)[i]))
		  next;
		}
		g<-AUC
		AUCs<-c(AUCs,g)	
		no.unnanotated.prots=length(unique((pred$labels)))-length(unique((train$labels)))
	}
	Results <- list("mean_AUCs" = mean(AUCs), "labels_with_AUC" = length(AUCs), "sd_AUCs" = sd(AUCs), "no.go.above.AUC80" = sum(AUCs>0.8), "no.unnanotated.prots" = no.unnanotated.prots, "total.prots" = length(unique((pred$labels))), "total.predicted.go" = length(unique((pred$go))))
	return(Results)
}
