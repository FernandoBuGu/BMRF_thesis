#Script to create a file with the genes that are neighbours of the gene of interest. These files are used in Step 4.

    #MSc thesis bioinfomratics WUR. Protein function prediction for poorly annotated species.
    #Author: Fernando Bueno Gutierrez
    #email1: fernando.buenogutierrez@wur.nl
    #email2: fernando.bueno.gutie@gmail.com


#The memory memory required in the HPC is 10

args<-commandArgs(T)


net<-read.table("big_topless035.tsv")

only_two_sets<-function(GENE_ID){
    #given a GENE_ID and a go_id, define X,NX,P,NP
    
    #GENE_ID: char, name of gene. GENE_ID exist in "net"
    #go_id: char, od of GO term. go_id shoul exists in "go"

    #L: list of char, each char is a set of chars (GENE_IDs)
        #X,NX,P,NP: gene, neighbours of gene, positive genes, neighbours of positive genes 

    X=GENE_ID
    neig<-as.character(net$V2[as.character(net$V1)==X])
    neig2<-as.character(net$V1[as.character(net$V2)==X])
    NX<-append(neig,neig2) 
    
    name<-paste("/mnt/nexenta/bueno002/part2_k10/NX/NX_",GENE_ID,sep="_")
    write.table(NX, file=name, col.names = F, row.names = F, quote = F, sep="\t")
}



#execute the function "only_two_sets" for each  the genes in the network file.
genes<-read.table("genes12424.tsv")

for(i in 1:length(genes$V1)){
    gene<-as.numeric(i)
    GENE_ID<-as.character(genes$V1[gene])
    sets=only_two_sets(GENE_ID)
}




#




