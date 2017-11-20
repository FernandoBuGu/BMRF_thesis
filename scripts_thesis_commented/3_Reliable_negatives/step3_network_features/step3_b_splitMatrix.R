#Script to extract a small subset of th esimilarity matrix. This subset corresponds to the rows and columns of a spicific GO term.

    #MSc thesis bioinfomratics WUR. Protein function prediction for poorly annotated species.
    #Author: Fernando Bueno Gutierrez
    #email1: fernando.buenogutierrez@wur.nl
    #email2: fernando.bueno.gutie@gmail.com


start.time <- Sys.time()

library(igraph)

library("parallel")
library("doParallel")
library("foreach")


dom<-read.table("big_ipr_labels_as_humans.tsv")
dom<-unique(dom)

net<-read.table("big_topless035.tsv")





goes_pass_filter<-read.table("take_30goes")

Matrix_GOSim<-read.table("AAAMatrix_sim_Goes.txt") #a similarity matrix of dimensions nºGOterms x nºGOterms

split_big_matrix<-function(i){
#retuns a similarity matrix where all the rows and columns refer to the GO term of interest. The matrix is a fragment of "Matrix_GOSim"

#In:
    #i: int, index of the GO term of interest in "goes_pass_filter"

#Out: 
    #thisGO: file, similarity matrix where all the rows and columns refer to the GO term of interest

    go_number<-as.numeric(i)	
    go_id<-as.character(goes_pass_filter$V1[go_number])

    MYGO_withdot<-as.character(go_id) 
    MYGO_withdot<-gsub(":", ".", MYGO_withdot)   

    thisGO<-Matrix_GOSim[,colnames(Matrix_GOSim)==MYGO_withdot]
    allgoes<-colnames(Matrix_GOSim)
    names(thisGO)<-allgoes
    thisGO<-as.data.frame(thisGO)
    name<-paste("/mnt/nexenta/bueno002/part2_k10/splitMatrix/_",go_id,sep="_")
    write.table(thisGO, file=name, col.names = F, row.names = T, quote = F, sep="\t")
}



#Run in parallel for 30 GO terms
cores=detectCores()
cl <- makeCluster(cores[1]-4)
registerDoParallel(cl)
clusterExport(cl,varlist=ls(),envir=environment())
w<-foreach(i = seq(1,30), .packages='Matrix', .combine="c") %dopar% { split_big_matrix(i) } #seq(1,30)
stopCluster(cl)
#










