source("https://bioconductor.org/biocLite.R")
biocLite("GO.db")
require("igraph")

#generate files uppropagate small
go_file<-read.table("GO_general.txt",header=T)
go_file <- go_file[c("label", "go", "term","evid")]

#small. Although to_uppropagated_small.uppropagated == to_uppropagated.uppropagated, cause uppropagate.r discards the nonEES
go_file_small = go_file[go_file$evid %in% c('EXP', 'IDA', 'IEP', 'IMP', 'IPI', 'IGI'),]
go_file_small = go_file[go_file$term == "P",]  #it does not make any difference
go_file_small$evid<-NULL
go_file_small<-unique(go_file_small)
write.table(go_file_small,file="to_uppropagated_small",sep="\t", col.names = F, row.names = F, quote = F)
source("uppropagate.r")
uppropagate(FILE="to_uppropagated_small")

#large
go_file$evid<-NULL
go_file<-unique(go_file)
write.table(go_file,file="to_uppropagated",sep="\t", col.names = F, row.names = F, quote = F)
uppropagate(FILE="to_uppropagated")











#small 
small_up <- read.table(file="to_uppropagated_small.uppropagated", na.strings=c("", "NA"), sep="\t")
small_up<-small_up[complete.cases(small_up),]
colnames(small_up)<-c("label","go","term")
DF<-read.table("DF") #required cause uppropagate creates 643 new proteins and some other things in the column. 
colnames(DF)[1]<-"label"
small_up<-small_up[small_up$label %in% DF$label,]
small_up$term<-NULL												#this can be removed
write.table(small_up,file="go_file_small.uppropagated_ready",sep="\t", col.names = F, row.names = F, quote = F)


	#info small
	tabl<-table(small_up)
	j<-as.matrix(tabl)
	h<-colSums(j)
	z<-sort(as.numeric(h))
	sum(z)
	mean(z)
	median(z)
	range(z)
	length(unique(small_up$go))
	length(unique(small_up$label))
	d <- density(z)
	jpeg('rplot_BP_EES_UPoooouuu.jpg')
	plot(d, main="Density #proteins/GO BP&EES UPPR22",xlim=range(0, 60, 1),xlab="#proteins for 5272 GO terms")
	axis(1, at = seq(0, 60, by = 3), las=2)
	dev.off()


#large
up <- read.table(file="to_uppropagated.uppropagated", na.strings=c("", "NA"), sep="\t")
up<-up[complete.cases(up),]
colnames(up)<-c("label","go","term")
DF<-read.table("DF") 					#required cause uppropagate creates some new proteins
colnames(DF)[1]<-"label"
up<-up[up$label %in% DF$label,]

##set EES based on small
colnames(up)<-c("label","go","term")
small_up <- read.table(file="to_uppropagated_small.uppropagated", na.strings=c("", "NA"), sep="\t")
small_up<-small_up[complete.cases(small_up),]
small_up$conc<-paste(small_up$V1,small_up$V2,sep="_")
up$conc<-paste(up$label,up$go,sep="_")
with_EES<-up[up$conc %in% small_up$conc,]
with_EES$conc<-"EES"
colnames(with_EES)[4]<-"evid"
without_EES<-up[!up$conc %in% small_up$conc,]
without_EES$conc<-"noEES"
colnames(without_EES)[4]<-"evid"
with_EES<-rbind(with_EES,without_EES)
up_evid<-with_EES

#cbind those that were removed in uppropagate.r
general<-read.table("GO_general.txt",header=T)
general <- general[c("label", "go", "term","evid")]
general_EES<-general[general$evid %in% c('EXP', 'IDA', 'IEP', 'IMP', 'IPI', 'IGI'),]
general_EES$evid<-"EES"
general_noEES<-general[!general$evid %in% c('EXP', 'IDA', 'IEP', 'IMP', 'IPI', 'IGI'),]
general_noEES$evid<-"noEES"
general_EES<-rbind(general_EES,general_noEES)
general_evid<-general_EES
up_evid<-rbind(up_evid,general_evid) 
up_evid<-unique(up_evid)

#This is required when some are going to enter the train and some others the test
#If a comb has EES, remove the one without EES. Also if a comb has P, remove the other one
up_evid$term <- sub("P", "0P", up_evid$term )
up_evid<-up_evid[ order(up_evid[,2], as.character(up_evid[,1]), up_evid[,4], up_evid[,3]), ]
up_evid<-up_evid[!(duplicated(up_evid[c("go","label","term")])),]
up_evid<-unique(up_evid)

#remove proteins that are not in network
up_evid <- up_evid[sample(nrow(up_evid)),]
up_evid <- up_evid[c("label", "go", "term","evid")]
write.table(up_evid,file="go_file.uppropagated_extended_ready",sep="\t", col.names = F, row.names = F, quote = F) 	#This have 4 columns but should not matter
up_evid$term<-NULL
up_evid$evid<-NULL

	#info large
	tabl<-table(up_evid)	
	j<-as.matrix(tabl)
	h<-colSums(j)
	z<-sort(as.numeric(h))
	sum(z)
	mean(z)
	median(z)
	range(z)
	length(unique(up_evid$go))
	length(unique(up_evid$label))
	d <- density(z)
	jpeg('rplot_UP.jpg')
	plot(d, main="Density #proteins/GO UPPR",xlim=range(0, 60, 1),xlab="#proteins for 8035 GO terms")
	axis(1, at = seq(0, 60, by = 3), las=2)
	dev.off()

