#Script to create the GO_term_files in each fold (after hiding the GO-gene-associations in the test set)
#The script also creates files "positive" (P) and "neighbours_of_positives" (NP) for each GO-term and each fold.

    #MSc thesis bioinfomratics WUR. Protein function prediction for poorly annotated species.
    #Author: Fernando Bueno Gutierrez
    #email1: fernando.buenogutierrez@wur.nl
    #email2: fernando.bueno.gutie@gmail.com


#The script is hard-coded for 30 GO terms, 10 folds and 4 replicates.


library("parallel")
library("doParallel")
library("foreach")


start.time <- Sys.time() #measure computational time
args<-commandArgs(T) #add arguments to sript
goes_pass_filter<-read.table("take_30goes") #import a file with (only) the GO terms. In this case, only 30 GO terms.
goes_pass_filter<-as.character(goes_pass_filter$V1)

go<<-read.table("/mnt/nexenta/bueno002/part2_k10/go_valid.tsv") #import GO_file
net<<-read.table("big_topless035.tsv")#import network file

source("/mnt/nexenta/bueno002/part2_k10/validation_functions_chickens.R")


hide_genes_AND_write_P_NP<-function(go_number,no_rep){
    #0) Creates the GO_term_files in each fold (after hiding the GO-gene-associations in the test set).
    #1-10) Creates files "positive" (P) and "neighbours_of_positives" (NP) for each GO-term and each fold.

    #In: 
        #go_number: int, from 1 to length of "goes_pass_filter", index of the GO term in the file "goes_pass_filter"
        #no_rep: int, from 1 to the number of replicates of the analysis. The analysis is computed independently for each replicate.

    #Out:
        #GO-files. One file for each fold, for the GOterm-replicate combination considered in the input.


    go_id<-as.character(goes_pass_filter[go_number]) #identification code of the GO term considered

    #0) create the GO file for analysis, after hidding reliable negatives (RN). One file per GO term
    table_goes<-read.table("/mnt/nexenta/bueno002/part2_k10/additional/take_30goes_table")
    extract_from_go<-go[go$V2 == go_id,] #select only the rows concerning the GO term considered
    extract_from_go_vect<-extract_from_go$V1 #extract the genes from these rows.
    labels<-rep(0,length(extract_from_go_vect)) #add label "0" to all genes, since the function "fold" requires a label.
    names(labels)<-extract_from_go_vect 
    labels_folds<-folds(labels,10) #creates 10 folds 
    extract_from_go_fold1<-extract_from_go[extract_from_go$V1 %in% as.character(names(unlist(labels_folds[1]))),] #extract genes each fold in labels_folds. These genes are going to be the test set.
    extract_from_go_fold2<-extract_from_go[extract_from_go$V1 %in% as.character(names(unlist(labels_folds[2]))),]
    extract_from_go_fold3<-extract_from_go[extract_from_go$V1 %in% as.character(names(unlist(labels_folds[3]))),]
    extract_from_go_fold4<-extract_from_go[extract_from_go$V1 %in% as.character(names(unlist(labels_folds[4]))),]
    extract_from_go_fold5<-extract_from_go[extract_from_go$V1 %in% as.character(names(unlist(labels_folds[5]))),]
    extract_from_go_fold6<-extract_from_go[extract_from_go$V1 %in% as.character(names(unlist(labels_folds[6]))),]
    extract_from_go_fold7<-extract_from_go[extract_from_go$V1 %in% as.character(names(unlist(labels_folds[7]))),]
    extract_from_go_fold8<-extract_from_go[extract_from_go$V1 %in% as.character(names(unlist(labels_folds[8]))),]
    extract_from_go_fold9<-extract_from_go[extract_from_go$V1 %in% as.character(names(unlist(labels_folds[9]))),]
    extract_from_go_fold10<-extract_from_go[extract_from_go$V1 %in% as.character(names(unlist(labels_folds[10]))),]
    


    go_fold1<-go[!rownames(go) %in% rownames(extract_from_go_fold1),]# for each fold, exclude the rows in the GO-file whose genes are in the the test set.
    go_fold2<-go[!rownames(go) %in% rownames(extract_from_go_fold2),] 
    go_fold3<-go[!rownames(go) %in% rownames(extract_from_go_fold3),] 
    go_fold4<-go[!rownames(go) %in% rownames(extract_from_go_fold4),] 
    go_fold5<-go[!rownames(go) %in% rownames(extract_from_go_fold5),] 
    go_fold6<-go[!rownames(go) %in% rownames(extract_from_go_fold6),] 
    go_fold7<-go[!rownames(go) %in% rownames(extract_from_go_fold7),] 
    go_fold8<-go[!rownames(go) %in% rownames(extract_from_go_fold8),] 
    go_fold9<-go[!rownames(go) %in% rownames(extract_from_go_fold9),] 
    go_fold10<-go[!rownames(go) %in% rownames(extract_from_go_fold10),] 
      
    name_fold1<-paste("/mnt/nexenta/bueno002/part2_k10/k10_gofile_fold1/",no_rep,"/gofile_h_",go_id,sep="") #name the GO-file after hidding the test set
    name_fold2<-paste("/mnt/nexenta/bueno002/part2_k10/k10_gofile_fold2/",no_rep,"/gofile_h_",go_id,sep="")
    name_fold3<-paste("/mnt/nexenta/bueno002/part2_k10/k10_gofile_fold3/",no_rep,"/gofile_h_",go_id,sep="")
    name_fold4<-paste("/mnt/nexenta/bueno002/part2_k10/k10_gofile_fold4/",no_rep,"/gofile_h_",go_id,sep="")
    name_fold5<-paste("/mnt/nexenta/bueno002/part2_k10/k10_gofile_fold5/",no_rep,"/gofile_h_",go_id,sep="")
    name_fold6<-paste("/mnt/nexenta/bueno002/part2_k10/k10_gofile_fold6/",no_rep,"/gofile_h_",go_id,sep="")
    name_fold7<-paste("/mnt/nexenta/bueno002/part2_k10/k10_gofile_fold7/",no_rep,"/gofile_h_",go_id,sep="")
    name_fold8<-paste("/mnt/nexenta/bueno002/part2_k10/k10_gofile_fold8/",no_rep,"/gofile_h_",go_id,sep="")
    name_fold9<-paste("/mnt/nexenta/bueno002/part2_k10/k10_gofile_fold9/",no_rep,"/gofile_h_",go_id,sep="")
    name_fold10<-paste("/mnt/nexenta/bueno002/part2_k10/k10_gofile_fold10/",no_rep,"/gofile_h_",go_id,sep="")
        
    write.table(go_fold1, file=name_fold1, col.names = F, row.names = F, quote = F, sep="\t") #save the GO-file after hidding the test set
    write.table(go_fold2, file=name_fold2, col.names = F, row.names = F, quote = F, sep="\t")
    write.table(go_fold3, file=name_fold3, col.names = F, row.names = F, quote = F, sep="\t")
    write.table(go_fold4, file=name_fold4, col.names = F, row.names = F, quote = F, sep="\t")
    write.table(go_fold5, file=name_fold5, col.names = F, row.names = F, quote = F, sep="\t")
    write.table(go_fold6, file=name_fold6, col.names = F, row.names = F, quote = F, sep="\t")
    write.table(go_fold7, file=name_fold7, col.names = F, row.names = F, quote = F, sep="\t")
    write.table(go_fold8, file=name_fold8, col.names = F, row.names = F, quote = F, sep="\t")
    write.table(go_fold9, file=name_fold9, col.names = F, row.names = F, quote = F, sep="\t")
    write.table(go_fold10, file=name_fold10, col.names = F, row.names = F, quote = F, sep="\t")



    #1) Creates files "positive" (P) and "neighbours_of_positives" (NP) for each GO-term and each fold.
    P<-as.character(go_fold1$V1[go_fold1$V2 == go_id]) #P (positives). The genes that are associated with the GO term
    NP<-c() #NP (neighbours of positives): The genes that are co-expressed with the genes in P.
    for(i in 1:length(P)){
        neig<-as.character(net$V2[net$V1==as.character(P[i])]) #co-expressed when the gene in P in in the 1st column of the network file
        neig2<-as.character(net$V1[net$V2==as.character(P[i])]) #co-expressed when the gene in P in in the 2nd column of the network file
        neig<-append(neig,neig2)
        NP<-append(NP,neig)
        NP<-unique(NP)
    }


    name1<-paste("/mnt/nexenta/bueno002/part2_k10/k10_P_NP_fold1/",no_rep,"/P_",go_id,sep="")
    name2<-paste("/mnt/nexenta/bueno002/part2_k10/k10_P_NP_fold1/",no_rep,"/NP_",go_id,sep="")
   

    write.table(P, file=name1, col.names = F, row.names = F, quote = F, sep="\t")
    write.table(NP, file=name2, col.names = F, row.names = F, quote = F, sep="\t")


    #2) Same as 1) but for another fold
    P<-as.character(go_fold2$V1[go_fold2$V2 == go_id]) 
    NP<-c()
    for(i in 1:length(P)){
        neig<-as.character(net$V2[net$V1==as.character(P[i])])
        neig2<-as.character(net$V1[net$V2==as.character(P[i])])
        neig<-append(neig,neig2)
        NP<-append(NP,neig)
        NP<-unique(NP)
    }


    name1<-paste("/mnt/nexenta/bueno002/part2_k10/k10_P_NP_fold2/",no_rep,"/P_",go_id,sep="")
    name2<-paste("/mnt/nexenta/bueno002/part2_k10/k10_P_NP_fold2/",no_rep,"/NP_",go_id,sep="")

    write.table(P, file=name1, col.names = F, row.names = F, quote = F, sep="\t")
    write.table(NP, file=name2, col.names = F, row.names = F, quote = F, sep="\t")



    #3) Same as 1) but for another fold
    P<-as.character(go_fold3$V1[go_fold3$V2 == go_id]) 
    NP<-c()
    for(i in 1:length(P)){
        neig<-as.character(net$V2[net$V1==as.character(P[i])])
        neig2<-as.character(net$V1[net$V2==as.character(P[i])])
        neig<-append(neig,neig2)
        NP<-append(NP,neig)
        NP<-unique(NP)
    }


    name1<-paste("/mnt/nexenta/bueno002/part2_k10/k10_P_NP_fold3/",no_rep,"/P_",go_id,sep="")
    name2<-paste("/mnt/nexenta/bueno002/part2_k10/k10_P_NP_fold3/",no_rep,"/NP_",go_id,sep="")

    write.table(P, file=name1, col.names = F, row.names = F, quote = F, sep="\t")
    write.table(NP, file=name2, col.names = F, row.names = F, quote = F, sep="\t")



    #4) Same as 1) but for another fold
    P<-as.character(go_fold4$V1[go_fold4$V2 == go_id]) 
    NP<-c()
    for(i in 1:length(P)){
        neig<-as.character(net$V2[net$V1==as.character(P[i])])
        neig2<-as.character(net$V1[net$V2==as.character(P[i])])
        neig<-append(neig,neig2)
        NP<-append(NP,neig)
        NP<-unique(NP)
    }


    name1<-paste("/mnt/nexenta/bueno002/part2_k10/k10_P_NP_fold4/",no_rep,"/P_",go_id,sep="")
    name2<-paste("/mnt/nexenta/bueno002/part2_k10/k10_P_NP_fold4/",no_rep,"/NP_",go_id,sep="")

    write.table(P, file=name1, col.names = F, row.names = F, quote = F, sep="\t")
    write.table(NP, file=name2, col.names = F, row.names = F, quote = F, sep="\t")



    #5) Same as 1) but for another fold
    P<-as.character(go_fold5$V1[go_fold5$V2 == go_id]) 
    NP<-c()
    for(i in 1:length(P)){
        neig<-as.character(net$V2[net$V1==as.character(P[i])])
        neig2<-as.character(net$V1[net$V2==as.character(P[i])])
        neig<-append(neig,neig2)
        NP<-append(NP,neig)
        NP<-unique(NP)
    }


    name1<-paste("/mnt/nexenta/bueno002/part2_k10/k10_P_NP_fold5/",no_rep,"/P_",go_id,sep="")
    name2<-paste("/mnt/nexenta/bueno002/part2_k10/k10_P_NP_fold5/",no_rep,"/NP_",go_id,sep="")

    write.table(P, file=name1, col.names = F, row.names = F, quote = F, sep="\t")
    write.table(NP, file=name2, col.names = F, row.names = F, quote = F, sep="\t")



    #6) Same as 1) but for another fold
    P<-as.character(go_fold6$V1[go_fold6$V2 == go_id]) 
    NP<-c()
    for(i in 1:length(P)){
        neig<-as.character(net$V2[net$V1==as.character(P[i])])
        neig2<-as.character(net$V1[net$V2==as.character(P[i])])
        neig<-append(neig,neig2)
        NP<-append(NP,neig)
        NP<-unique(NP)
    }


    name1<-paste("/mnt/nexenta/bueno002/part2_k10/k10_P_NP_fold6/",no_rep,"/P_",go_id,sep="")
    name2<-paste("/mnt/nexenta/bueno002/part2_k10/k10_P_NP_fold6/",no_rep,"/NP_",go_id,sep="")

    write.table(P, file=name1, col.names = F, row.names = F, quote = F, sep="\t")
    write.table(NP, file=name2, col.names = F, row.names = F, quote = F, sep="\t")


    #7) Same as 1) but for another fold
    P<-as.character(go_fold7$V1[go_fold7$V2 == go_id]) 
    NP<-c()
    for(i in 1:length(P)){
        neig<-as.character(net$V2[net$V1==as.character(P[i])])
        neig2<-as.character(net$V1[net$V2==as.character(P[i])])
        neig<-append(neig,neig2)
        NP<-append(NP,neig)
        NP<-unique(NP)
    }


    name1<-paste("/mnt/nexenta/bueno002/part2_k10/k10_P_NP_fold7/",no_rep,"/P_",go_id,sep="")
    name2<-paste("/mnt/nexenta/bueno002/part2_k10/k10_P_NP_fold7/",no_rep,"/NP_",go_id,sep="")

    write.table(P, file=name1, col.names = F, row.names = F, quote = F, sep="\t")
    write.table(NP, file=name2, col.names = F, row.names = F, quote = F, sep="\t")

    #8) Same as 1) but for another fold
    P<-as.character(go_fold8$V1[go_fold8$V2 == go_id]) 
    NP<-c()
    for(i in 1:length(P)){
        neig<-as.character(net$V2[net$V1==as.character(P[i])])
        neig2<-as.character(net$V1[net$V2==as.character(P[i])])
        neig<-append(neig,neig2)
        NP<-append(NP,neig)
        NP<-unique(NP)
    }

    name1<-paste("/mnt/nexenta/bueno002/part2_k10/k10_P_NP_fold8/",no_rep,"/P_",go_id,sep="")
    name2<-paste("/mnt/nexenta/bueno002/part2_k10/k10_P_NP_fold8/",no_rep,"/NP_",go_id,sep="")

    write.table(P, file=name1, col.names = F, row.names = F, quote = F, sep="\t")
    write.table(NP, file=name2, col.names = F, row.names = F, quote = F, sep="\t")


    #9) Same as 1) but for another fold
    P<-as.character(go_fold9$V1[go_fold9$V2 == go_id]) 
    NP<-c()
    for(i in 1:length(P)){
        neig<-as.character(net$V2[net$V1==as.character(P[i])])
        neig2<-as.character(net$V1[net$V2==as.character(P[i])])
        neig<-append(neig,neig2)
        NP<-append(NP,neig)
        NP<-unique(NP)
    }


    name1<-paste("/mnt/nexenta/bueno002/part2_k10/k10_P_NP_fold9/",no_rep,"/P_",go_id,sep="")
    name2<-paste("/mnt/nexenta/bueno002/part2_k10/k10_P_NP_fold9/",no_rep,"/NP_",go_id,sep="")

    write.table(P, file=name1, col.names = F, row.names = F, quote = F, sep="\t")
    write.table(NP, file=name2, col.names = F, row.names = F, quote = F, sep="\t")

    #10) get NP,P fold10
    P<-as.character(go_fold10$V1[go_fold10$V2 == go_id]) 
    NP<-c()
    for(i in 1:length(P)){
        neig<-as.character(net$V2[net$V1==as.character(P[i])])
        neig2<-as.character(net$V1[net$V2==as.character(P[i])])
        neig<-append(neig,neig2)
        NP<-append(NP,neig)
        NP<-unique(NP)
    }


    name1<-paste("/mnt/nexenta/bueno002/part2_k10/k10_P_NP_fold10/",no_rep,"/P_",go_id,sep="")
    name2<-paste("/mnt/nexenta/bueno002/part2_k10/k10_P_NP_fold10/",no_rep,"/NP_",go_id,sep="")

    write.table(P, file=name1, col.names = F, row.names = F, quote = F, sep="\t")
    write.table(NP, file=name2, col.names = F, row.names = F, quote = F, sep="\t")


    end.time <- Sys.time()
    time.taken <- end.time - start.time
    print("time:")
    print(time.taken)


}



#run the 4 replicates in parallel for 30 GO terms.
cores=detectCores()
cl <- makeCluster(cores[1]-4)
registerDoParallel(cl)
clusterExport(cl,varlist=ls(),envir=environment())
for(kk in 1:30){
    w<-foreach(i = seq(1,4), .packages='Matrix', .combine="c") %dopar% { hide_genes_AND_write_P_NP(kk,i) }
}
stopCluster(cl)





