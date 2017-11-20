#Script to annotate GO-speific features for each gene. These features are used to extract Reliable negatives (RN) in PU-BMRF.

    #MSc thesis bioinfomratics WUR. Protein function prediction for poorly annotated species.
    #Author: Fernando Bueno Gutierrez
    #email1: fernando.buenogutierrez@wur.nl
    #email2: fernando.bueno.gutie@gmail.com




bt_inter<-function(inp1,inp2){

    #returns a list of numerics, corresponding to the annotation of GO-speific features for each gene

    #In:
        #inp1: char of char, genes. Either X or NX
        #inp2: char of char, genes. Either P or NP

    #Out:
        # bt_inter_sum: sum of the betweeness of the genes that are in both sets (inp1 and inp2)
        # bt_inter_mean: mean of the betweeness of the genes that are in both sets (inp1 and inp2)
        # cl_inter_sum: sum of the closeness of the genes that are in both sets (inp1 and inp2)
        # cl_inter_mean: mean of the closeness of the genes that are in both sets (inp1 and inp2)
        # tr_inter_sum: sum of the transitivity of the genes that are in both sets (inp1 and inp2)
        # tr_inter_mean: mean of the transitivity of the genes that are in both sets (inp1 and inp2)

    
    names_inter<-as.character(inp1[inp1 %in% inp2])


    if(length(names_inter) != 0){ #As long as there is a at least one gene common in inp1 adn inp2.

        #small
        bt_inter_sum<-round(sum(as.numeric(graph_betweenness[names(graph_betweenness) %in% names_inter]),na.rm=T),3)
        bt_inter_mean<-round(mean(as.numeric(graph_betweenness[names(graph_betweenness) %in% names_inter]),na.rm=T),3)
        cl_inter_sum<-round(sum(as.numeric(graph_closeness[names(graph_closeness) %in% names_inter]),na.rm=T),3)
        cl_inter_mean<-round(mean(as.numeric(graph_closeness[names(graph_closeness) %in% names_inter]),na.rm=T),3)

        #large
        bt_inter_large_sum<-round(sum(as.numeric(graph_betweenness_large[names(graph_betweenness_large) %in% names_inter]),na.rm=T),3)
        bt_inter_large_mean<-round(mean(as.numeric(graph_betweenness_large[names(graph_betweenness_large) %in% names_inter]),na.rm=T),3)
        cl_inter_large_sum<-round(sum(as.numeric(graph_closeness_large[names(graph_closeness_large) %in% names_inter]),na.rm=T),3)
        cl_inter_large_mean<-round(mean(as.numeric(graph_closeness_large[names(graph_closeness_large) %in% names_inter]),na.rm=T),3)
        tr_inter_large_sum<-round(sum(as.numeric(graph_transitivity_large[names(graph_transitivity_large) %in% names_inter]),na.rm=T),3)
        tr_inter_large_mean<-round(mean(as.numeric(graph_transitivity_large[names(graph_transitivity_large) %in% names_inter]),na.rm=T),3)
    } else {
        bt_inter_sum=NA
        bt_inter_mean=NA
        cl_inter_sum=NA
        cl_inter_mean=NA

        bt_inter_large_sum=NA
        bt_inter_large_mean=NA
        cl_inter_large_sum=NA
        cl_inter_large_mean=NA
        tr_inter_large_sum=NA
        tr_inter_large_mean=NA
    }

    L<-list("bt_inter_sum"=bt_inter_sum,"bt_inter_mean"=bt_inter_mean,"cl_inter_sum"=cl_inter_sum,"cl_inter_mean"=cl_inter_mean,"bt_inter_large_sum"=bt_inter_large_sum,"bt_inter_large_mean"=bt_inter_large_mean,"cl_inter_large_sum"=cl_inter_large_sum,"cl_inter_large_mean"=cl_inter_large_mean,"tr_inter_large_sum"=tr_inter_large_sum,"tr_inter_large_mean"=tr_inter_large_mean)
    return(L)
}



bt_UNL<-function(inp1){

#    inp1: char of char, genes. Either X or NX
#    inp2: char of char, genes. Either P or NP

    bt_UNL_sum<-round(sum(as.numeric(graph_betweenness_UNL[names(graph_betweenness_UNL) %in% inp1]),na.rm=T),3)
    bt_UNL_mean<-round(mean(as.numeric(graph_betweenness_UNL[names(graph_betweenness_UNL) %in% inp1]),na.rm=T),3)
    cl_UNL_sum<-round(sum(as.numeric(graph_closeness_UNL[names(graph_closeness_UNL) %in% inp1]),na.rm=T),3)
    cl_UNL_mean<-round(mean(as.numeric(graph_closeness_UNL[names(graph_closeness_UNL) %in% inp1]),na.rm=T),3)
    tr_UNL_sum<-round(sum(as.numeric(graph_transitivity_UNL[names(graph_transitivity_UNL) %in% inp1]),na.rm=T),3)
    tr_UNL_mean<-round(mean(as.numeric(graph_transitivity_UNL[names(graph_transitivity_UNL) %in% inp1]),na.rm=T),3)



    L<-list("bt_UNL_sum"=bt_UNL_sum,"bt_UNL_mean"=bt_UNL_mean,"cl_UNL_sum"=cl_UNL_sum,"cl_UNL_mean"=cl_UNL_mean,"tr_UNL_sum"=tr_UNL_sum,"tr_UNL_mean"=tr_UNL_mean)
    return(L)
}






common<-function(inp1,inp2){ #common: genes that are both in inp1 and in inp2.
    #given two sets of genes, find genes that are in both sets, and the portion of genes in set1 that are in set2
    
    #inp1: char (X or NX)
    #inp2: char (P or NP)
    
    #L: list with 2 num
        #no_genes_common: num, #genes in both sets
        #p_genes_inp1_common: num, portion of genes in set1 that are in set2

    no_genes_common<-length(intersect(inp1,inp2))
    p_genes_inp1_common<-round(no_genes_common*100/length(inp1),2)
    p_genes_inp2_common<-round(no_genes_common*100/length(inp2),2)
    L=list("no_genes_common"=no_genes_common,"p_genes_inp1_common"=p_genes_inp1_common,"p_genes_inp2_common"=p_genes_inp2_common)
    return(L)
}




links<-function(inp1,inp2){
    #given two ses of genes, find edges that link genes from both sets, and the portion of edges from genes in set1 that link to genes in set2

    #inp1: char (X or NX)
    #inp2: char (P or NP)

    #L: list with 2 num
        #edges_inp1_inp2: num, #edges that link genes from both sets
        #p_genes_inp1_lintTp_inp2: num, portion of edges from genes in set1 that link to genes in set2

    total_edges_from_inp1_1<-net[net$V1 %in% inp1,] 
    total_edges_from_inp1_2<-net[net$V2 %in% inp1,] 
    total_edges_from_inp1<-unique(rbind( total_edges_from_inp1_1,total_edges_from_inp1_2))   
    total_edges_from_inp1<-dim(total_edges_from_inp1)[1]

    edges_inp1_inp2_1<-net[net$V1 %in% inp1 & net$V2 %in% inp2,]
    edges_inp1_inp2_2<-net[net$V1 %in% inp2 & net$V2 %in% inp1,]
    edges_inp1_inp2<-unique(rbind(edges_inp1_inp2_1,edges_inp1_inp2_2))
    edges_inp1_inp2=dim(edges_inp1_inp2)[1]

    p_genes_inp1_lintTp_inp2=round(edges_inp1_inp2*100/total_edges_from_inp1,2)
    L=list("edges_inp1_inp2"=edges_inp1_inp2,"p_genes_inp1_lintTp_inp2"=p_genes_inp1_lintTp_inp2)
    return(L)
}




domains<-function(inp1,inp2){
    #given two ses of genes, find domains that are present in genes from both sets. 

    #inp1: char (X or NX)
    #inp2: char (P or NP)

    #L: list with 3 num
        #no_shared_assoc: num, #matches gene-domain. Some of these matches may be repeated
        #no_genes_share_domains: num, #genes that have domains that are in both sets
        #no_shared_domains: nom, #domains that are in both sets

    #note: If a gene is in both sets and has a domain, that domain will also be considerd as match

    dom_inp1<-dom[dom$V1 %in% inp1,]
    dom_inp2<-dom[dom$V1 %in% inp2,]
    shared_domain_asso<-c() 
    genes_share_domains<-c()
    if(dim(dom_inp1)[1]>0 & dim(dom_inp2)[1]>0){
        for(i in 1:dim(dom_inp1)[1]){
            if(dom_inp1$V2[i] %in% dom_inp2$V2){
                shared_domain_asso<-append(shared_domain_asso,as.character(dom_inp1$V2[i]))
                genes_share_domains<-append(genes_share_domains,as.character(dom_inp2$V1[i]))}
        }
        no_shared_assoc=length(shared_domain_asso)                      
        no_genes_share_domains=length(unique(genes_share_domains))      
        no_shared_domains=length(unique(shared_domain_asso))       
    } else {
        no_shared_assoc=0                    
        no_genes_share_domains=0
        no_shared_domains=0          
    }
    L<-list("no_shared_assoc"=no_shared_assoc,"no_genes_share_domains"=no_genes_share_domains,"no_shared_domains"=no_shared_domains)
    return(L)
}


write_gene<-function(GENE_ID,go_id){
    #given two ses of genes, return a vector of num, with info of that gene in differnet GO-sepcific features  

    #GENE_ID: char, name of gene. GENE_ID exist in "net"
    #go_id: char, od of GO term. go_id shoul exists in "go"

    #info: df with 1 row, each column is a numeric, with info of the gene (GENE_ID) in differnet GO-sepcific features    


    name<-paste("/home/WUR/bueno002/part2_k10/NX/NX_",GENE_ID,sep="_")
    NX<-read.table(name)
    NX<-as.character(NX$V1)
    X=GENE_ID
    name1<-paste("/home/WUR/bueno002/part2_k10/k10_P_NP_fold",fol,"/",as.numeric(no_rep),"/P_",go_id,sep="")
    P=read.table(name1)
    P=as.character(P$V1)
    name2<-paste("/home/WUR/bueno002/part2_k10/k10_P_NP_fold",fol,"/",as.numeric(no_rep),"/NP_",go_id,sep="")
    NP=read.table(name2)
    NP=as.character(NP$V1)

        #common (x4)
#    source("k10_features_funct_improved_foldX.R")
#    common<-common(NX,P)
#    no_genes_common_NXP<-common$no_genes_common
#    p_genes_inp1_common_NXP<-common$p_genes_inp1_common

#    source("k10_features_funct_improved_foldX.R")
#    common<-common(X,NP)
#    no_genes_common_XNP<-common$no_genes_common

    source("k10_features_funct_improved_foldX.R")
    common<-common(NX,NP)
    no_genes_common_NXNP<-common$no_genes_common
    p_genes_inp1_common_NXNP<-common$p_genes_inp1_common
    p_genes_inp2_common_NXNP<-common$p_genes_inp2_common


       #domains (x4)
    source("k10_features_funct_improved_foldX.R")
    domains<-domains(X,P)
    no_shared_assoc_XP=domains$no_shared_assoc 

    source("k10_features_funct_improved_foldX.R")
    domains<-domains(NX,P)
    no_shared_assoc_NXP=domains$no_shared_assoc
    no_genes_share_domains_NXP=domains$no_genes_share_domains

    source("k10_features_funct_improved_foldX.R")
    domains<-domains(X,NP)
    no_shared_assoc_XNP=domains$no_shared_assoc
    no_genes_share_domains_XNP=domains$no_genes_share_domains

    source("k10_features_funct_improved_foldX.R")
    domains<-domains(NX,NP)
    no_shared_assoc_NXNP=domains$no_shared_assoc
    no_genes_share_domains_NXNP=domains$no_genes_share_domains
    no_shared_domains_NXNP=domains$no_shared_domains


        #R_inter
    source("k10_features_funct_improved_foldX.R")
    R_interNXNP<-bt_inter(NX,NP)
    bt_interNXNP_s<-R_interNXNP$bt_inter_sum
    #bt_interNXNP_m<-R_interNXNP$bt_inter_mean. This gives computational problems
    cl_interNXNP_s<-R_interNXNP$cl_inter_sum
    cl_interNXNP_m<-R_interNXNP$cl_inter_mean

        #large
    bt_interNXNP_s_large<-R_interNXNP$bt_inter_large_sum
    bt_interNXNP_m_large<-R_interNXNP$bt_inter_large_mean
    cl_interNXNP_s_large<-R_interNXNP$cl_inter_large_sum
    #cl_interNXNP_m_large<-R_interNXNP$cl_inter_large_mean. This gives computational problems
    tr_interNXNP_s_large<-R_interNXNP$tr_inter_large_sum
    tr_interNXNP_m_large<-R_interNXNP$tr_inter_large_mean


#        #R_UNL
    source("k10_features_funct_improved_foldX.R")
    R_UNLNXP<-bt_UNL(X)
    bt_UNLNXP_s<-R_UNLNXP$bt_UNL_sum
    bt_UNLNXP_m<-R_UNLNXP$bt_UNL_mean
    #cl_UNLNXP_s<-R_UNLNXP$cl_UNL_sum
    cl_UNLNXP_m<-R_UNLNXP$cl_UNL_mean
    tr_UNLNXP_s<-R_UNLNXP$tr_UNL_sum
    tr_UNLNXP_m<-R_UNLNXP$tr_UNL_mean

    R_UNLNXP<-bt_UNL(NX)
    bt_UNLNXP_s_NX<-R_UNLNXP$bt_UNL_sum
    bt_UNLNXP_m_NX<-R_UNLNXP$bt_UNL_mean
    cl_UNLNXP_s_NX<-R_UNLNXP$cl_UNL_sum
    #cl_UNLNXP_m_NX<-R_UNLNXP$cl_UNL_mean
    tr_UNLNXP_s_NX<-R_UNLNXP$tr_UNL_sum
    tr_UNLNXP_m_NX<-R_UNLNXP$tr_UNL_mean



        #weighted neigh      
    goes_neigh<-as.character(go$V2[go$V1 %in% NX])
    name<-paste("/home/WUR/bueno002/part2_k10/splitMatrix/_",go_id,sep="_")
    Matrix_GOSim<-read.table(name)
    goes_neigh<-unique(gsub(":", ".", goes_neigh)) 
    goes_neigh_wei<-as.numeric(Matrix_GOSim$V1[as.character(Matrix_GOSim$V1) %in% goes_neigh])
    sum_we_GO_neigh<-round(sum(goes_neigh_wei),3)
    mean_we_GO_neigh<-round(mean(goes_neigh_wei,na.rm=T),3)

        #weighted thisGene     
    goes_this<-as.character(go$V2[go$V1 %in% X])
    goes_this<-unique(gsub(":", ".", goes_this)) 
    goes_this_wei<-as.numeric(Matrix_GOSim$V1[as.character(Matrix_GOSim$V1) %in% goes_this])
    sum_we_GO_this<-round(sum(goes_this_wei),3)






#   no_GO_neig and no_GO_v_neigS   
    start.time.neigh <- Sys.time()
    geness<-as.character(go$V1)
    G<-table(geness)
    no_GO_neig<-G[names(G) %in% NX]
    no_GO_neig<-sum(no_GO_neig,na.rm=T)
    end.time <- Sys.time()
     start.time.neigh <- end.time -  start.time.neigh
    print("time neigh")
    print(start.time.neigh)



        #to vector
    info_common<-c(no_genes_common_NXNP,p_genes_inp1_common_NXNP,p_genes_inp2_common_NXNP)
    names_common<-c("no_genes_common_NXNP","p_genes_inp1_common_NXNP","p_genes_inp2_common_NXNP")


    info_R_inter<-c(bt_interNXNP_s,cl_interNXNP_s,cl_interNXNP_m)
    names_R_inter<-c("bt_interNXNP_s","cl_interNXNP_s","cl_interNXNP_m")

    info_R_inter_large<-c(bt_interNXNP_s_large,bt_interNXNP_m_large,cl_interNXNP_s_large,tr_interNXNP_s_large,tr_interNXNP_m_large)
    names_R_inter_large<-c("bt_interNXNP_s_large","bt_interNXNP_m_large","cl_interNXNP_s_large","tr_interNXNP_s_large","tr_interNXNP_m_large")


    info_R_UNL<-c(bt_UNLNXP_s,bt_UNLNXP_m,tr_UNLNXP_s,tr_UNLNXP_m)
    names_R_UNL<-c("bt_UNLNXP_s","bt_UNLNXP_m","tr_UNLNXP_s","tr_UNLNXP_m")


    info_R_UNL_NX<-c(bt_UNLNXP_s_NX,bt_UNLNXP_m_NX,cl_UNLNXP_s_NX,tr_UNLNXP_s_NX,tr_UNLNXP_m_NX)
    names_R_UNL_NX<-c("bt_UNLNXP_s_NX","bt_UNLNXP_m_NX","cl_UNLNXP_s_NX","tr_UNLNXP_s_NX","tr_UNLNXP_m_NX")

#    info_edge<-c(edges_inp1_inp2_XP,p_genes_inp1_lintTp_inp2_XP,edges_inp1_inp2_NXP,p_genes_inp1_lintTp_inp2_NXP,edges_inp1_inp2_XNP,p_genes_inp1_lintTp_inp2_XNP,edges_inp1_inp2_NXNP,p_genes_inp1_lintTp_inp2_NXNP)
#    names_edge<-c("edges_inp1_inp2_XP","p_genes_inp1_lintTp_inp2_XP","edges_inp1_inp2_NXP","p_genes_inp1_lintTp_inp2_NXP","edges_inp1_inp2_XNP","p_genes_inp1_lintTp_inp2_XNP","edges_inp1_inp2_NXNP","p_genes_inp1_lintTp_inp2_NXNP")

    info_domains<-c(no_shared_assoc_XP,no_shared_assoc_NXP,no_genes_share_domains_NXP,no_shared_assoc_XNP,no_genes_share_domains_XNP,no_shared_assoc_NXNP,no_genes_share_domains_NXNP,no_shared_domains_NXNP)
    names_domains<-c("no_shared_assoc_XP","no_shared_assoc_NXP","no_genes_share_domains_NXP","no_shared_assoc_XNP","no_genes_share_domains_XNP","no_shared_assoc_NXNP","no_genes_share_domains_NXNP","no_shared_domains_NXNP")


    info_weigh_NEig_GO<-c(sum_we_GO_neigh,mean_we_GO_neigh,sum_we_GO_this)
    names_weigh_NEig_GO<-c("sum_we_GO_neigh","mean_we_GO_neigh","sum_we_GO_this")

    info_noGO_neigh<-c(no_GO_neig)
    names_noGO_neigh<-c("no_GO_neig")




    info<-(GENE_ID)
    info<-append(info,info_common)
#    info<-append(info,info_edge)
    info<-append(info,info_R_inter)
    info<-append(info,info_R_inter_large)
    info<-append(info,info_R_UNL)
    info<-append(info,info_R_UNL_NX)
    info<-append(info,info_domains)
    info<-append(info,info_weigh_NEig_GO)
    info<-append(info,info_noGO_neigh)

    names<-"GENE_ID"
    names<-append(names,names_common)
#    names<-append(names,names_edge)
    names<-append(names,names_R_inter_large)
    names<-append(names,names_R_inter)
    names<-append(names,names_R_UNL)
    names<-append(names,names_R_UNL_NX)
    names<-append(names,names_domains)
    names<-append(names,names_weigh_NEig_GO)
    names<-append(names,names_noGO_neigh)

    info<-as.data.frame(info)
    info<-t(info)
    colnames(info)<-as.character(names)

    return(info)
}


write_database_go<-function(GENE_ID,GOID,no_rep){
    #given a GOID, writes a df. Each row is a gene in the network and each column is a value corresponding to one GO-sepcific features   

    #GOID: char, od of GO term. GOID shoul exists in "go"

    #database: df with X row, where X is the number of genes in network each column is a numeric, with info of the gene (GENE_ID) in differnet GO-sepcific features


    #small network 
    name_a<-paste("/home/WUR/bueno002/part2_k10/A/","f",fol,"/",as.numeric(no_rep),"/bt_small",GOID,".Rda",sep="")
    name_b<-paste("/home/WUR/bueno002/part2_k10/A/","f",fol,"/",as.numeric(no_rep),"/cl_small",GOID,".Rda",sep="")
    graph_betweenness<<-readRDS(name_a)
    graph_closeness<<-readRDS(name_b)

    #large network 
    name_a<-paste("/home/WUR/bueno002/part2_k10/A/","f",fol,"/",as.numeric(no_rep),"/bt_large",GOID,".Rda",sep="")
    name_b<-paste("/home/WUR/bueno002/part2_k10/A/","f",fol,"/",as.numeric(no_rep),"/cl_large",GOID,".Rda",sep="")
    name_c<-paste("/home/WUR/bueno002/part2_k10/U/","f",fol,"/",as.numeric(no_rep),"/tr_large",GOID,".Rda",sep="")
    graph_betweenness_large<<-readRDS(name_a)
    graph_closeness_large<<-readRDS(name_b)
    graph_transitivity_large<<-readRDS(name_c)

    #huge network 
    name_aa<-paste("/home/WUR/bueno002/part2_k10/U/","f",fol,"/",as.numeric(no_rep),"/bt_large",GOID,".Rda",sep="")
    name_bb<-paste("/home/WUR/bueno002/part2_k10/U/","f",fol,"/",as.numeric(no_rep),"/cl_large",GOID,".Rda",sep="")
    name_cc<-paste("/home/WUR/bueno002/part2_k10/U/","f",fol,"/",as.numeric(no_rep),"/tr_large",GOID,".Rda",sep="")
    graph_betweenness_UNL<<-readRDS(name_aa)
    graph_closeness_UNL<<-readRDS(name_bb)
    graph_transitivity_UNL<<-readRDS(name_cc)
   

    row=write_gene(GENE_ID,GOID)
    database<-data.frame(row)
    name_file<-paste("/home/WUR/bueno002/part2_k10/feat/",fol,"/",as.numeric(no_rep),"/",GOID,GENE_ID,sep="")      
    write.table(database, file=name_file, col.names = T, row.names = F, quote = F, sep="\t")
}


    










