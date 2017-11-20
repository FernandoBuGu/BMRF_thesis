setwd("/home/bueno002/table1/r_tissues")
ti<-read.table("r_tissues_file.txt",header=T)
names_all<-colnames(ti)
names<-colnames(ti)[2:36]

ti<-ti[with(ti, order(go_number)), ]

ti_goes<-as.character(ti$go_number)
t<-table(ti_goes)
goes_no4<-names(which(t!=4))#change to 5


ti <- ti[!ti$go_number %in% goes_no4,]


d<-seq(1,length(ti$go_number)+1,4)

df_means<-data.frame()
df_sd<-data.frame()

for(i in 1:(length(d)-1)){

	
	start=d[i]
	end=d[i+1]-1
	ncols<-dim(ti)[2]	
	

	ti_go<-ti[start:end,2:ncols]
	colmeans<-colMeans(ti_go)
	colmeans<-round(colmeans,3)
	df_means<-rbind(df_means,colmeans)
	colnames(df_means)<-names
	SDs<-apply(ti_go, 2, sd)
	SDs<-round(SDs,3)
	df_sd<-rbind(df_sd,SDs)
	colnames(df_sd)<-names
}

df_means<-cbind(unique(ti$go_number),df_means)
df_sd<-cbind(unique(ti$go_number),df_sd)
colnames(df_means)<-names_all
colnames(df_sd)<-names_all
write.table(df_means, file="AUCtissues_perGO.txt", col.names = T, row.names = F, quote = F, sep="\t")
M<-read.table("/home/bueno002/table1/reduce_corrs_perGO.txt",header=T)

Mer<-merge(M,df_means,by="go_number",all.x=T,all.y=T)


####################################################################################################33Add info
info<-read.table("/home/bueno002/table1/info_perGO.txt",header=T) 
info$epn<-NULL
info2<-read.table("/home/bueno002/table1/info_perGO2.txt",header=T,fill=T) 
info2 <- subset(info2, !duplicated(info2[,1]))
info2<-info2[,c(1,6)]
info<-merge(info,info2,by="go_number",all.x=T,all.y=F)

Mer<-Mer[,1:7]   #To remove the AUC.info of each individual tissue
Mer2<-merge(Mer,info,by="go_number",all.x=T,all.y=T)


#add here the others
info_from2<-read.table("/home/bueno002/table2/info_from_table2_to1.txt",header=T) 
info_from2$Terms<-NULL

Mer3<-merge(Mer2,info_from2,by="go_number",all.x=T,all.y=F)


Mer3<-Mer3[,-c(1,19)]

#https://stackoverflow.com/questions/22282531/how-to-compute-correlations-between-all-columns-in-r-and-detect-highly-correlate
#https://stackoverflow.com/questions/10680658/how-can-i-create-a-correlation-matrix-in-r

library(dplyr)
library(reshape2)
Mer3<-Mer3[complete.cases(Mer3), ]

Mer4<-sapply(Mer3,as.numeric)
d_cor <- as.matrix(cor(Mer4))
d_cor_melt <- arrange(melt(d_cor), -abs(value))

d_cor_melt<-d_cor_melt[-c(1:dim(Mer4)[2]),]
d_cor_melt_AUC<-d_cor_melt[d_cor_melt$Var1=="AUC",]
write.table(d_cor_melt_AUC, file="/home/bueno002/table1/corr_variables_AUC.txt", col.names = T, row.names = F, quote = F, sep="\t")

S<-seq(1,dim(d_cor_melt)[1],2)#remove duplicates
d_cor_melt2<-d_cor_melt[-S,]

write.table(d_cor_melt2, file="/home/bueno002/table1/corr_bt_variables.txt", col.names = T, row.names = F, quote = F, sep="\t")













