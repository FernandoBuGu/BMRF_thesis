#This script creates the 3 .tsv files that are required to run https://github.com/jwbargsten/bmrf
########################################At one point it requires to go out from R and call Python

#Input data for yeast comes from:
	#Coexpression network: http://www.inetbio.org/yeastnet/downloadnetwork.php -> .txt files
	#GO: http://www.yeastgenome.org/download-data/curation -> gene_association.sgd
	#uniprot codes: http://www.uniprot.org/docs/yeast: (saveas)-> yeast.txt (gedit,calc)-> Uniprot_calc.txt


install.packages<-c("plyr","reshape","gsubfn","AUC","parallel","pbapply")

#Combine coexpression networks
filenames <- list.files(path="path_to_txt_files", full.names=TRUE)
library(plyr)
import.list <- llply(filenames, read.table)
u<-sapply(import.list,dim)
library(reshape)
data <- merge_recurse(import.list)
data[3]<-NULL
dat.sort = t(apply(data, 1, sort))
data<-data[!duplicated(dat.sort),]
data$V3<-1
write.table(data,file="coexpression_all.txt",sep="\t", col.names = F, row.names = F, quote = F)

#Obtain Uniprot codes to extract domains (only for labels that are in the network)
expre = read.table("coexpression_all.txt",sep="\t")
all_labels <- unique(c(as.character(expre[,1]),as.character(expre[,2])))	#domain and GO files should not contain any label that is not in this vector
DF<-as.data.frame(rep(1,length(all_labels)))
DF$label<-all_labels
DF[1]<-NULL
Uniprot <- read.csv("Uniprot_calc.txt", sep="", header=F)
colnames(Uniprot)[1]<-"label"
colnames(Uniprot)[2]<-"Unip"
Codes_labels_inData<-merge(DF,Uniprot,by="label")	#This ensures that we do not get domain info for proteins that are not in the network
Codes_labels_inData[1]<-NULL
write.table(DF,file="DF",sep="\t", col.names = F, row.names = F, quote = F)
write.table(Codes_labels_inData,file="Codes_labels_inData",sep="\t", col.names = F, row.names = F, quote = F)
write.table(expre,file="big_topless.tsv",sep="\t", col.names = F, row.names = F, quote = F)

#~ python create_data_frame.py Codes_labels_inData			This takes around 30 minutes

#python output -> domain file
domains = read.table("domains",sep=",")
domains=tail(domains, -1)
library("reshape")
colnames(domains)[1]<-"Unip"
DD<-melt.data.frame(domains, "Unip")
DD<-DD[with(DD, order(Unip)), ]
DD<-DD[!(is.na(DD$value) | DD$value==""), ]
DD$variable<-NULL

#UniprotCode/Domain -> label/Domain
Uniprot <- read.csv("Uniprot_calc.txt", sep="", header=F)
colnames(Uniprot)[1]<-"label"
colnames(Uniprot)[2]<-"Unip"
prot_domain<-merge(DD,Uniprot,by="Unip")
prot_domain[1]<-NULL
prot_domain$domain<-prot_domain$value
prot_domain[1]<-NULL
colnames(prot_domain)[1]<-"label"
write.table(prot_domain,file="big_ipr_labels.tsv",sep="\t", col.names = F, row.names = F, quote = F)

#Generate GO file
setwd("/home/WUR/bueno002/From_gitHub/large_coex")
library(gsubfn)
GO = read.csv("gene_association.sgd",skip = 8,sep="\t", header=F, quote="")
GO <- data.frame(GO$V5,GO$V7,GO$V9,GO$V11)
colnames(GO) <- c("go", "evid", "term","label")   
GO$label <- sub('(?<=\\|).*$', '', GO$label, perl=TRUE)
GO$label<-gsubfn(".", list("|" = ""), GO$label)
write.table(GO,file="GO_general.txt",sep="\t", col.names = T, row.names = F, quote = T)





