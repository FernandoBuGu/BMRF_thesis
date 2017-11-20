#450, normal and rest have sometimes goven resut of same GOterm: remove the second or more lines of same GO with same reduce code 
##############################ONE SINGLE PLOT WITH MEAN(AUC) IN THE DIFFERENT SITUATIONS OF REDUCE...NOISE...

#Parse info so that I get mean AUC over 5, and sd
info<-read.table("/home/bueno002/table1/r_info/inform.txt",header=T,fill=T) #some random numbers appear between GO and GO (i.e line 6428)


	#some cleaning
info <- info[!is.na(info$no_labels),] #remove incomplete rows

info <- info[!is.na(info$epp_v.epn_v),]
goes<-as.character(info$go_number)
nchar_goes<-nchar(goes)
goes<-data.frame(goes)
goes2<-cbind(goes,nchar_goes)
remove<-which(goes2$nchar_goes>10)
info <- info[ -remove, ]


info <- info[!is.na(info$AUC),] #remove incomplete rows


	#mean for each GO
goes<-as.character(info$go_number)
t<-table(goes)
goes_no5<-names(which(t!=5))




info_complete <- info[!info$go_number %in% goes_no5,]

info_complete<-info_complete[with(info_complete, order(go_number)), ]


d<-seq(1,length(info_complete$AUC)+1,5)
means<-c()
sds<-c()
for(i in 1:(length(d)-1)){
	start=d[i]
	end=d[i+1]-1
	u<-mean(as.vector(info_complete$AUC[start:end]))
	sd<-sd(as.vector(info_complete$AUC[start:end]))
	means<-append(means,u)
	sds<-append(sds,sd)}

d<-d[-length(d)]



info_complete<-info_complete[ d, ] 
info_complete$AUC<-NULL
info_complete$AUC<-round(means,3)
info_complete$sd<-round(sds,3)


#add AUC reduce 0 to the reduce files


info<-info_complete
write.table(info, file="/home/bueno002/table1/info_perGO.txt", col.names = T, row.names = F, quote = F, sep="\t")
toMerge<-data.frame(info$go_number)
toMerge<-cbind(toMerge,info$AUC)
toMerge<-cbind(toMerge,info$sd)
colnames(toMerge)<-c("go_number","0","sd_0")


################################################################################################################################################
################################################################################################################################################



setwd("/home/bueno002/table1/OPEN_oa_epp_epn_enn")
filenames <- list.files(path="/home/bueno002/table1/OPEN_oa_epp_epn_enn", full.names=TRUE) #Correct manually, the first line
library(plyr)
library(reshape)
import.list <- llply(filenames, read.table)
data <- merge_recurse(import.list)

colnames(data)<-c("go_number","reduce","0.05","sd_0.05","0.1","sd_0.1","0.2","sd_0.2","0.4","sd_0.4","0.6","sd_0.6","0.8","sd_0.8","0.99","sd_0.99")

data_previous<-data[!(duplicated(data[c("go_number","reduce")]) | duplicated(data[c("go_number","reduce")], fromLast = TRUE)), ]


DA<-merge(data_previous,toMerge,by="go_number",all.x=T,all.y=T)
DA<-DA[with(DA, order(go_number,reduce)), ]
DA <- DA[c("go_number","reduce", "0", "sd_0","0.05","sd_0.05","0.1","sd_0.1","0.2","sd_0.2","0.4","sd_0.4","0.6","sd_0.6","0.8","sd_0.8","0.99","sd_0.99")]  #DA: the complete data

DA <- DA[!is.na(DA$reduce),] #remove incomplete rows
data<-DA

#add a column that is the correlation
	#keep only the means
data[,4]<-NULL
data[,5]<-NULL
data[,6]<-NULL
data[,7]<-NULL
data[,8]<-NULL
data[,9]<-NULL
data[,10]<-NULL
data[,11]<-NULL

return_means<-function(redu){
	#redu="oa"
	data<-data[data$reduce==redu,]  

	corrs<-c()
	for(i in 1:dim(data)[1]){
		co<-data[i,]
		co<-co[,3:length(co)]
		corrs<-append(corrs,cor(as.numeric(co),as.numeric(c(0,0.05,0.1,0.2,0.4,0.6,0.8,0.99))))}

	data <- DA
	data<-data[data$reduce==redu,]   
	data$corrs<-corrs
	colnames(data)<-c("go_number","reduce","0", "sd_0","0.05","sd_0.05","0.1","sd_0.1","0.2","sd_0.2","0.4","sd_0.4","0.6","sd_0.6","0.8","sd_0.8","0.99","sd_0.99","corr")

		#file
	data_redu_r<-data.frame(as.character(data$go_number))
	data_redu_r<-cbind(data_redu_r,round(data$corr,3))
	colnames(data_redu_r)<-c("go_number","corr_redu")


		#means
	x<-dim(data)[2]
	D<-data[,3:x]
	no_columns<-dim(D)[2]

	for(i in 1:no_columns){
		D[,i] <- as.numeric(D[,i])}
	corr<-colMeans(D,na.rm = T)[dim(D)[2]]
	lista<-list("corr"=corr,"all_means"=colMeans(D,na.rm = T),"file"=data_redu_r)
	return(lista)
}


out_oa<-return_means("oa")
out_epp<-return_means("epp")
out_epn<-return_means("epn")
out_enn<-return_means("enn")
out_oa$corr  #mean of correlations AUC_oa
out_epp$corr  #mean of correlations AUC_epp
out_epn$corr  #mean of correlations AUC_epn
out_enn$corr  #mean of correlations AUC_enn
#
corr_after_means<-function(out_redu){
	out_redu$all_means_ready<-out_redu$all_means[seq(1,length(out_redu$all_means)-1,2)]
	corr<-cor(as.numeric(out_redu$all_means_ready),as.numeric(c(0,0.05,0.1,0.2,0.4,0.6,0.8,0.99)))   
	return(corr)}

corr_after_means(out_oa)    #mean of correlations AUC_oa   (doing first the average)
corr_after_means(out_epp)   #mean of correlations AUC_epp  (doing first the average)
corr_after_means(out_epn)  #mean of correlations AUC_epn   (doing first the average)
corr_after_means(out_enn)   #mean of correlations AUC_enn  (doing first the average)

write.table(out_oa$file, file="corr_oa.txt", col.names = T, row.names = F, quote = F, sep="\t")
write.table(out_epp$file, file="corr_epp.txt", col.names = T, row.names = F, quote = F, sep="\t")
write.table(out_epn$file, file="corr_epn.txt", col.names = T, row.names = F, quote = F, sep="\t")
write.table(out_enn$file, file="corr_enn.txt", col.names = T, row.names = F, quote = F, sep="\t")



################################################################################################################################################
################################################################################################################################################
setwd("/home/bueno002/table1/OPEN_amg")
filenames <- list.files(path="/home/bueno002/table1/OPEN_amg", full.names=TRUE)
library(plyr)
library(reshape)
import.list <- llply(filenames, read.table)
data <- merge_recurse(import.list)
colnames(data)<-c("go_number","reduce","0.05","sd_0.05","0.1","sd_0.1","0.2","sd_0.2","0.4","sd_0.4","0.6","sd_0.6","0.8","sd_0.8")

data_previous<-data[!(duplicated(data[c("go_number","reduce")]) | duplicated(data[c("go_number","reduce")], fromLast = TRUE)), ]

DA<-merge(data_previous,toMerge,by="go_number",all.x=T,all.y=T)
DA<-DA[with(DA, order(go_number,reduce)), ]
DA <- DA[c("go_number","reduce", "0", "sd_0","0.05","sd_0.05","0.1","sd_0.1","0.2","sd_0.2","0.4","sd_0.4","0.6","sd_0.6","0.8","sd_0.8")]  #DA: the complete data

DA <- DA[!is.na(DA$reduce),] #remove incomplete rows
DA <- DA[!DA$reduce!="amg",] #remove incomplete rows
data<-DA


#add a column that is the correlation
	#keep only the means
data[,4]<-NULL
data[,5]<-NULL
data[,6]<-NULL
data[,7]<-NULL
data[,8]<-NULL
data[,9]<-NULL
data[,10]<-NULL


return_means<-function(redu){
	redu="amg"
	data<-data[data$reduce==redu,]  
	
	data <- data[!is.na(data[5]),] #remove incomplete rows. 0.1 is always there

	corrs<-c()
	for(i in dim(data)[1]){		
		co<-data[i,]
		co<-co[,3:length(co)]
		corrs<-append(corrs,cor(as.numeric(co),as.numeric(c(0,0.05,0.1,0.2,0.4,0.6,0.8)),use = "complete.obs",method = c("pearson")))}





	data <- DA
	data <- data[!is.na(data[5]),] 
	colnames(data)<-c("go_number","reduce","0","sd_0","0.05","sd_0.05","0.1","sd_0.1","0.2","sd_0.2","0.4","sd_0.4","0.6","sd_0.6","0.8","sd_0.8")
	data<-data[data$reduce==redu,]  
	data$corrs<-corrs
	colnames(data)<-c("go_number","reduce","0","sd_0","0.05","sd_0.05","0.1","sd_0.1","0.2","sd_0.2","0.4","sd_0.4","0.6","sd_0.6","0.8","sd_0.8","corr")

	#file
	data_r<-data.frame(as.character(data$go_number))
	data_r<-cbind(data_r,round(data$corr,3))
	colnames(data_r)<-c("go_number","corr_amg")

	#means
	x<-dim(data)[2]
	D<-data[3:x]
	no_columns<-dim(D)[2]

	for(i in 1:no_columns){
		D[,i] <- as.numeric(D[,i])}
	corr<-colMeans(D,na.rm = T)[dim(D)[2]]
	lista<-list("corr"=corr,"all_means"=colMeans(D,na.rm = T),"file"=data_r)
	return(lista)
}

out_amg=return_means("amg")
out_amg$corr   #mean of correlations AUC_amg of each GO term
out_amg$all_means

out_amg$all_means_ready<-out_amg$all_means[seq(1,length(out_amg$all_means)-1,2)]
cor(as.numeric(out_amg$all_means_ready),as.numeric(c(0,0.05,0.1,0.2,0.4,0.6,0.8)))   #mean of correlations AUC_amg (doing first the average)
write.table(out_amg$file, file="corr_amg.txt", col.names = T, row.names = F, quote = F, sep="\t")


################################################################################################################################################
################################################################################################################################################
setwd("/home/bueno002/table1/OPEN_noise")
filenames <- list.files(path="/home/bueno002/table1/OPEN_noise", full.names=TRUE) #Correct manually, the first line
library(plyr)
library(reshape)
import.list <- llply(filenames, read.table)
data <- merge_recurse(import.list)

colnames(data)<-c("go_number","reduce","0.1","sd_0.1","0.5","sd_0.5","1","sd_1","2","sd_2","4","sd_4","8","sd_8","16","sd_16")

data_previous<-data[!(duplicated(data[c("go_number","reduce")]) | duplicated(data[c("go_number","reduce")], fromLast = TRUE)), ]

DA<-merge(data_previous,toMerge,by="go_number",all.x=T,all.y=T)
DA<-DA[with(DA, order(go_number,reduce)), ]
DA <- DA[c("go_number","reduce", "0", "sd_0","0.1","sd_0.1","0.5","sd_0.5","1","sd_1","2","sd_2","4","sd_4","8","sd_8","16","sd_16")]  #DA: the complete data

DA <- DA[!is.na(DA$reduce),] #remove incomplete rows
DA <- DA[!DA$reduce!="noise",] #remove incomplete rows
data<-DA

data[,4]<-NULL
data[,5]<-NULL
data[,6]<-NULL
data[,7]<-NULL
data[,8]<-NULL
data[,9]<-NULL
data[,10]<-NULL
data[,11]<-NULL


return_means<-function(redu){
	redu="noise"
	data<-data[data$reduce==redu,]  
	data <- data[!is.na(data[5]),] 

	corrs<-c()
	for(i in 1:dim(data)[1]){
		co<-data[i,]
		co<-co[,3:length(co)]
		corrs<-append(corrs,cor(as.numeric(co),as.numeric(c(0,0.1,0.5,1,2,4,8,16)),use = "complete.obs",method = c("pearson")))}

	data <- DA
	data <- data[!is.na(data[5]),] 
	colnames(data)<-c("go_number","reduce","0","sd_0","0.1","sd_0.1","0.5","sd_0.5","1","sd_1","2","sd_2","4","sd_4","8","sd_8","16","sd_16")
	data<-data[data$reduce==redu,]  
	data$corrs<-corrs
	colnames(data)<-c("go_number","reduce","0","sd_0","0.1","sd_0.1","0.5","sd_0.5","1","sd_1","2","sd_2","4","sd_4","8","sd_8","16","sd_16","corr")



	#file
	data_r<-data.frame(as.character(data$go_number))
	data_r<-cbind(data_r,round(data$corr,3))
	colnames(data_r)<-c("go_number","corr_noise")

	#means
	x<-dim(data)[2]
	D<-data[3:x]
	no_columns<-dim(D)[2]
	for(i in 1:no_columns){
		D[,i] <- as.numeric(D[,i])}
	corr<-colMeans(D,na.rm = T)[dim(D)[2]]
	lista<-list("corr"=corr,"all_means"=colMeans(D,na.rm = T),"file"=data_r)

	return(lista)
}

out_noise=return_means("noise")
out_noise$corr   #mean of correlations AUC_noise of each GO term
out_noise$all_means

out_noise$all_means_ready<-out_noise$all_means[seq(1,length(out_noise$all_means)-1,2)]
cor(as.numeric(out_noise$all_means_ready),as.numeric(c(0,0.1,0.5,1,2,4,8,16)))   #mean of correlations AUC_noise (doing first the average)
write.table(out_noise$file, file="corr_noise.txt", col.names = T, row.names = F, quote = F, sep="\t")

################################################################################################################################################
################################################################################################################################################

noise<-read.table("/home/bueno002/table1/OPEN_noise/corr_noise.txt",header=T)
amg<-read.table("/home/bueno002/table1/OPEN_amg/corr_amg.txt",header=T)
epp<-read.table("/home/bueno002/table1/OPEN_oa_epp_epn_enn/corr_epp.txt",header=T)
oa<-read.table("/home/bueno002/table1/OPEN_oa_epp_epn_enn/corr_oa.txt",header=T)
epn<-read.table("/home/bueno002/table1/OPEN_oa_epp_epn_enn/corr_epn.txt",header=T)
enn<-read.table("/home/bueno002/table1/OPEN_oa_epp_epn_enn/corr_enn.txt",header=T)
M<-merge(noise,amg,by="go_number",all.x=T,all.y=T)
M<-merge(M,epp,by="go_number",all.x=T,all.y=T)
M<-merge(M,epn,by="go_number",all.x=T,all.y=T)
M<-merge(M,enn,by="go_number",all.x=T,all.y=T)
M<-merge(M,oa,by="go_number",all.x=T,all.y=T)
colnames(M)<-c("go_number","corr_noise","corr_amg","corr_epp","corr_epn","corr_enn","corr_oa")
write.table(M, file="reduce_corrs_perGO.txt", col.names = T, row.names = F, quote = F, sep="\t")
Mn<-M[,2:7]
colMeans(Mn,na.rm = T)













