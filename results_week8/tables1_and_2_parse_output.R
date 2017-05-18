#the output from pbreplicate needs some parsing, and therefore this code. It is ugly but it works (I have checked). 

#for 50 go and 15 rep
w_unlist<-unlist(w)
nm<-seq(1,1515,1)
names(w_unlist)<-as.character(nm)
rm<-seq(1,1515,101)
valid<-w_unlist[!names(w_unlist) %in% rm]
names<-c()
poss<-c()
all_15r<-c()
all_sd_15r<-c()
poss_sd<-c()
for(i in 1:50) {
	#start=start[-i]+2
	name<-c()
	name<-paste("go",i,sep="_")
	names<-append(names,name)
	pos<-c()
	start<-seq(2,1515,2)
	pos<-seq(start[i],1515,101)
	poss<-append(poss,pos)
	s<-15*(i-1)+1
	e<-s+14
	u<-poss[s:e]
	for_this<-mean(valid[names(valid) %in% u])
	all_15r<-append(all_15r,for_this)
	pos_sd<-seq((start[i]+1),1515,101)
	poss_sd<-append(poss_sd,pos_sd)
	s<-15*(i-1)+1
	e<-s+14
	u<-poss_sd[s:e]
	for_this<-mean(valid[names(valid) %in% u])
	all_sd_15r<-append(all_sd_15r,for_this)
}



#for 50go and 15 rep taking only 10 replicate
w_unlist<-unlist(w)
nm<-seq(1,1515,1)
names(w_unlist)<-as.character(nm)
rm<-seq(1,1515,101) 
valid<-w_unlist[!names(w_unlist) %in% rm]
names<-c()
poss<-c()
all_10r<-c()
all_sd_10r<-c()
poss_sd<-c()
for(i in 1:50) {
	#start=start[-i]+2
	name<-c()
	name<-paste("go",i,sep="_")
	names<-append(names,name)
	pos<-c()
	start<-seq(2,1515,2)
	pos<-seq(start[i],1010,101)		
	poss<-append(poss,pos)
	s<-10*(i-1)+1
	e<-s+9
	u<-poss[s:e]
	for_this<-mean(valid[names(valid) %in% u])
	all_10r<-append(all_10r,for_this)
	pos_sd<-seq((start[i]+1),1010,101)
	poss_sd<-append(poss_sd,pos_sd)
	s<-10*(i-1)+1
	e<-s+9
	u<-poss_sd[s:e]
	for_this<-mean(valid[names(valid) %in% u])
	all_sd_10r<-append(all_sd_10r,for_this)
}	


#for 50go and 15 rep taking only 5 replicate
w_unlist<-unlist(w)
nm<-seq(1,1515,1)
names(w_unlist)<-as.character(nm)
rm<-seq(1,1515,101) 
valid<-w_unlist[!names(w_unlist) %in% rm]
names<-c()
poss<-c()
all_5r<-c()
all_sd_5r<-c()
poss_sd<-c()
for(i in 1:50) {
	#start=start[-i]+2
	name<-c()
	name<-paste("go",i,sep="_")
	names<-append(names,name)
	start<-seq(2,1515,2)
	pos<-seq(start[i],505,101)		
	poss<-append(poss,pos)
	s<-5*(i-1)+1
	e<-s+4
	u<-poss[s:e]
	for_this<-mean(valid[names(valid) %in% u])
	all_5r<-append(all_5r,for_this)
	pos_sd<-seq((start[i]+1),505,101)
	poss_sd<-append(poss_sd,pos_sd)
	s<-5*(i-1)+1
	e<-s+4
	u<-poss_sd[s:e]
	for_this<-mean(valid[names(valid) %in% u])
	all_sd_5r<-append(all_sd_5r,for_this)
}	


table1<-as.data.frame(paste("go",seq(1:50),sep="_"))
table1$mean_5_rep<-all_5r
table1$sd_5_rep<-all_sd_5r
table1$mean_10_rep<-all_10r
table1$sd_10_rep<-all_sd_10r
table1$mean_15_rep<-all_15r
table1$sd_15_rep<-all_sd_15r
colnames(table1)[1]<-"go"
table1_2<-rbind(table1,last)
table1["means" ,-1] <- colMeans(table1[,-1])
table1[,-1] <-round(table1[,-1],4)
colnames(table1)<-c("go","mean_AUC_5r","sd_AUC_5r","mean_AUC_10r","sd_AUC_10r","mean_AUC_15r","sd_AUC_15r")
write.table(table1,file="table1.txt",sep="\t", col.names = T, row.names = F, quote = F)

table2 <- table1[-nrow(table1),]
table2<-table2[,1:3]

dim_goes<-readRDS(file="dim_goes.Rda")
mean_conex<-readRDS(file="mean_conex.Rda")
mean_conex<-round(mean_conex,4)
table2<-cbind(table2,dim_goes$no_prots)
table2<-cbind(table2,dim_goes$no_conexion)
table2<-cbind(table2,mean_conex)
colnames(table2)<-c("go","mean_AUC_5r","sd_AUC_5r","#proteins","#conexions","mean_#conexions")
table2<-table2[with(table2, order(-mean_AUC_5r)), ]
write.table(table2,file="table2.txt",sep="\t", col.names = T, row.names = F, quote = F)
