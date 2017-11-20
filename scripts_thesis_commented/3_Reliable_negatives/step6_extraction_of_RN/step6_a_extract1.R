# Script to prepare the Reliable negative (RN) data for each GO term


    #MSc thesis bioinfomratics WUR. Protein function prediction for poorly annotated species.
    #Author: Fernando Bueno Gutierrez
    #email1: fernando.buenogutierrez@wur.nl
    #email2: fernando.bueno.gutie@gmail.com


args<-commandArgs(T)
library(plyr)
library(reshape)
library(AUC)	
goes_pass_filter<<-read.table("take_30goes")

source("/home/WUR/bueno002/chickens035/validation_functions_chickens.R")

create_RN_file<-function(init_tol=1,minAUC=1,step=0.1,no_rep=1,fol=1,go_number){

#Creates an ".R image" object that can be used in step 6b. The object contains the information regarding the RN in the train set for teh GO term considered 

#In:
    #init_tol: int, tolerance at start. value 1 corresponds to the default value in Figure 8 of thesis. 
    #minAUC: int from 0 to 1, minimum AUC
    #step: int, difference in tolerance from one attempt of convergence to the next
    #no_rep: int, nº replicate
    #fol: int, number of fold
    #go_number: int, index of the GO term from "take_30goes"

#Out:
    #nameDATA: ".R image" object containing the information regarding the RN in the train set for teh GO term considered. 

    myGO<<-as.character(goes_pass_filter$V1[go_number])
    print(myGO)
    nameXX<<-paste("/home/WUR/bueno002/part2_k10/k10_gofile_fold",fol,"/",no_rep,"/gofile_h_",myGO,sep="")
    go<<-read.table(nameXX)

    genes_myGO<<-as.character(go$V1[go$V2==myGO]) #note that these files are only "valid". Note also that some genes have been hidden
        
    Nspec<<-read.table("AAAnoGOspec_info",header=T)


    #imprort the features files from all genes, for the GO term of interest
    name.dir=paste("/home/WUR/bueno002/part2_k10/feat/",fol,"/",no_rep,sep="")
    setwd(name.dir)
    lf <- list.files(name.dir, pattern=as.character(myGO))
    stopifnot(length(lf)>0)
    import.list <- llply(lf, read.table, header=T)


    #split import.list in bits of 200 files, since teh .R function "merge_recurse" cannot handle a larger number of files.
    till<-round(length(import.list)/200,0)
    start=0
    end=0
    for(i in 1:till){

        if(end !=0){
            start=end+1
        } else {
            start=1}


        if(i==till){
            end=length(lf)     
        } else {
            end=start+200}
        X<-import.list[start:end]
        if(i ==1){ #For 4 of the features we assign always value "NULL", because the fetures were defined at a given point but they will not be used.
            data<-merge_recurse(X)
            data$cl_interNXNP_m_large<-NULL
            data$bt_interNXNP_m<-NULL
            data$cl_UNLNXP_s<-NULL
            data$cl_UNLNXP_m_NX<-NULL
        } else {
            bit <- merge_recurse(X)
            bit$cl_interNXNP_m_large<-NULL
            bit$bt_interNXNP_m<-NULL
            bit$cl_UNLNXP_s<-NULL
            bit$cl_UNLNXP_m_NX<-NULL
            data<-rbind(data,bit)}
    }
    colnames(data)[1]<-"gene"

    features_genes<<-merge(Nspec,data,by="gene") #merge the NONspec with the spec fetures

    #scale the features
    name_gene<<-features_genes$gene
    features_genes<<-features_genes[,2:dim(features_genes)[2]]
    scalar1 <<- function(x) {x / sqrt(sum(x^2,na.rm=T))}
    features_genes<<-sapply(features_genes, function(x) scalar1(x)) 
    features_genes<<-as.data.frame(features_genes)
    features_genes<<-sapply(features_genes, function(x) round(x,4)) 
    features_genes<<-as.data.frame(features_genes)
    features_genes<<-cbind(as.character(name_gene),features_genes)
    features_genes<<-as.data.frame(features_genes)
    colnames(features_genes)[1]<<-"gene"
        
    #P_train: Positives train
    P_train<<-as.character(features_genes$gene[features_genes$gene %in% genes_myGO]) #These are the positives of the train

    #P_test: Positives test
    go<<-read.table("/home/WUR/bueno002/part2_k10/go_valid.tsv")
    P_test<<-as.character(go$V1[go$V2 == myGO & !go$V1 %in% P_train])

        #U_train: Unlabelled train. If unlabelled train is defined, unlabelled test can be automaitcally defined. 
    U_all_and_hidden<<-as.character(features_genes$gene[!features_genes$gene %in% P_train])
    U_all<<-as.character(features_genes$gene[!features_genes$gene %in% P_train & !features_genes$gene %in% P_test])
    S<-rep(0,length(U_all))
    names(S)<-U_all
    source("/home/WUR/bueno002/chickens035/validation_functions_chickens.R")
    U_all_FOLDS<-folds(S) 
    FO<-seq(1,10)
    FO<-FO[!FO %in% fol]
    U_train<-U_all_FOLDS[c(FO)]
    nameDATA<<-paste("/home/WUR/bueno002/part2_k10/step5a/",fol,"/",no_rep,"/",myGO,".RData",sep="")
    save.image(file=nameDATA)
}

sol<<-create_RN_file(init_tol=1,minAUC=1,step=0.1,no_rep=as.numeric(args[1]),fol=as.numeric(args[2]),go_number=as.numeric(args[3]))

