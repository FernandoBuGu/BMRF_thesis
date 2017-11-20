#Script to create networks used in "step 3A - Network features" of the thesis

    #MSc thesis bioinfomratics WUR. Protein function prediction for poorly annotated species.
    #Author: Fernando Bueno Gutierrez
    #email1: fernando.buenogutierrez@wur.nl
    #email2: fernando.bueno.gutie@gmail.com


#It takes around 38 minutes per fold


start.time <- Sys.time()
args<-commandArgs(T)

library("parallel")
library("doParallel")
library("foreach")

goes_pass_filter<-read.table("take_30goes")

net<-read.table("big_topless035.tsv")

getNet<-function(go_number,no_rep){
# Saves RDS objects corresponding to the networks described in "step 3A - Network features"

#In:
    #go_number: int, index of the GO term in "goes_pass_filter"
    #no_rep: int, nº replicate. Analysis is independent for each replicate. 

#Out: 
    #graph_betweenness,graph_closeness,AAA: betweenessm closeness and transitivity as described in "step 3A - Network features".
    

    go_id<-as.character(goes_pass_filter[go_number,])



    name_fold1<-paste("/mnt/nexenta/bueno002/part2_k10/k10_gofile_fold1/",no_rep,"/gofile_h_",go_id,sep="")
    go<-read.table(name_fold1)

    #generate betweeness, closeness and transitivty vector for the GO term

    nameP<-paste("/mnt/nexenta/bueno002/part2_k10/k10_P_NP_fold1/",no_rep,"/P_",go_id,sep="")
    P<-read.table(nameP)

    nameNP<-paste("/mnt/nexenta/bueno002/part2_k10/k10_P_NP_fold1/",no_rep,"/NP_",go_id,sep="")
    NP<-read.table(nameNP)



    library(igraph)
    DA_net_P<-net[net$V1 %in% P$V1 | net$V2 %in% P$V1,] #a network for all edges that have at least one node in P. Here P and U are possible, thus it helps in differinating P and U
    D<-get.adjacency(graph.edgelist(as.matrix(DA_net_P), directed=FALSE))
    aracne.net_as.graph<-graph_from_adjacency_matrix(D)
    graph_betweenness<-betweenness(aracne.net_as.graph,v=V(aracne.net_as.graph))
    graph_closeness<-closeness(aracne.net_as.graph,v=V(aracne.net_as.graph))
    #the transitivity here was 0 in almost all cases



    saveRDS(graph_betweenness, file=paste("/mnt/nexenta/bueno002/part2_k10/A/f1/",no_rep,"/bt_small",go_id,".Rda",sep=""))
    saveRDS(graph_closeness, file=paste("/mnt/nexenta/bueno002/part2_k10/A/f1/",no_rep,"/cl_small",go_id,".Rda",sep=""))




    DA_net_NP<-net[net$V1 %in% NP$V1 & net$V2 %in% NP$V1,] 
    DA_net_P<-net[net$V1 %in% P$V1 | net$V2 %in% P$V1,] 
    DA_net_P_large<-rbind(DA_net_NP,DA_net_P)  #a network for all edges whose both nodes are in NP, and all edges that have at least one node in P
    D<-get.adjacency(graph.edgelist(as.matrix(DA_net_P_large), directed=FALSE))
    aracne.net_as.graph<-graph_from_adjacency_matrix(D)
    graph_betweenness<-betweenness(aracne.net_as.graph,v=V(aracne.net_as.graph))
    graph_closeness<-closeness(aracne.net_as.graph,v=V(aracne.net_as.graph))
    AAA<-transitivity(aracne.net_as.graph,type="local")
    names(AAA)<-names(graph_closeness)


    saveRDS(graph_betweenness, file=paste("/mnt/nexenta/bueno002/part2_k10/A/f1/",no_rep,"/bt_large",go_id,".Rda",sep=""))
    saveRDS(graph_closeness, file=paste("/mnt/nexenta/bueno002/part2_k10/A/f1/",no_rep,"/cl_large",go_id,".Rda",sep=""))
    saveRDS(AAA, file=paste("/mnt/nexenta/bueno002/part2_k10/A/f1/",no_rep,"/tr_large",go_id,".Rda",sep=""))



    #huge netork. Perhaps the positives are more interconnected. Perhaps, in a network of mostly unlabelled, the RN have a very high cenctricty and th epositives a very low,
    start.time.huge <- Sys.time()
    D<-get.adjacency(graph.edgelist(as.matrix(net), directed=FALSE))
    aracne.net_as.graph<-graph_from_adjacency_matrix(D)
    graph_betweenness<-betweenness(aracne.net_as.graph,v=V(aracne.net_as.graph))
    graph_closeness<-closeness(aracne.net_as.graph,v=V(aracne.net_as.graph))
    AAA<-transitivity(aracne.net_as.graph,type="local")
    names(AAA)<-names(graph_closeness)


    saveRDS(graph_betweenness, file=paste("/mnt/nexenta/bueno002/part2_k10/U/f1/",no_rep,"/bt_large",go_id,".Rda",sep=""))
    saveRDS(graph_closeness, file=paste("/mnt/nexenta/bueno002/part2_k10/U/f1/",no_rep,"/cl_large",go_id,".Rda",sep=""))
    saveRDS(AAA, file=paste("/mnt/nexenta/bueno002/part2_k10/U/f1/",no_rep,"/tr_large",go_id,".Rda",sep=""))
    end.time <- Sys.time()
    time.taken <- end.time - start.time.huge
    print("time both")
    print(time.taken)



}



getNet(as.numeric(args[1]),as.numeric(args[2]))






end.time <- Sys.time()
time.taken <- end.time - start.time
print("time both")
print(time.taken)
