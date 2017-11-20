#Script to invetsigate whether the genes are more coexpressed in the specific functions than in the more general functions

a<-read.table("info_h.txt",header=T,fill=T)
a<-a[complete.cases(a[ , 13]),]
go_file<-read.table("rbind_upvalid_upNONvalid_H")
net<-read.table("big_topless_H.tsv")
a<-a[a$spec>5 & a$genesV<3000,]

m<-median(a$spec,na.rm=T) #21

spec_goes<-a$go[a$spec<7]
gen_goes<-a$go[a$spec>2700]
print(length(gen_goes))



genes_only_spec<-unique(go_file$V1[!go_file$V2 %in% gen_goes & go_file$V2 %in% spec_goes])
genes_only_gen<-unique(go_file$V1[!go_file$V2 %in% spec_goes & go_file$V2 %in% gen_goes])

    #remove genes that are in both groups
genes_only_spec<-genes_only_spec[!genes_only_spec %in% genes_only_gen]
length(genes_only_spec)
genes_only_gen<-genes_only_gen[!genes_only_gen %in% genes_only_spec]


#are they more interconnected?
    #equal size
difs<-c()
for(i in 1:10000){
    start.time <- Sys.time()
    genes_only_spec_S<-sample(genes_only_spec,200)
    genes_only_gen_S<-sample(genes_only_gen,200)
    N_spe<-net[net$V1 %in% genes_only_spec_S & net$V1 %in% genes_only_spec_S,]
    N_gen<-net[net$V1 %in% genes_only_gen_S & net$V1 %in% genes_only_gen_S,]
    dim(N_gen)
    dim(N_spe)
    dif<-(dim(N_spe)[1]-dim(N_gen)[1])*100/dim(N_gen)[1]
    difs<-append(difs,dif)
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    time.taken
}




#are they more interconnected when they are positives with other positives
