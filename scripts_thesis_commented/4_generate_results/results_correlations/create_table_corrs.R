x<-read.table("/home/bueno002/table2/tissue_names") #2x2 (i.e. 1  Stomach; 2	Esophagus-Muscularis...)



x2<-read.table("/home/bueno002/table1/r_tissues/AUCtissues_perGO.txt",header=T) #1103x36 (go_number::AUC_tissue35)
x3<-x2[,2:dim(x2)[2]] #1103x35 (AUC_tissue1::AUC_tissue35)
means<-colMeans(x3,na.rm = T)
df<-data.frame(x$V2,round(means,4))
df<-df[with(df, order(-means)), ] #35x2 (Stomach	0.7138860...) 			mean of each tissue, sorted by mean
write.table(df, file="/home/bueno002/table2/all_tissues_same.txt", col.names = T, row.names = F, quote = F, sep="\t")

#1) Generate "table2.txt": best and worst tissue for each GO term

	#1a) find the worst tissue (and AUC) for each GO 

min_f<-function(ddff){
	u<-min(ddff,na.rm = T)
	return(u)}




cc<-dim(x3)[2]
x3after<-x3[rowSums(is.na(x3))!=cc, ]
out<-as.numeric(rownames(x3)[!rownames(x3) %in% rownames(x3after)]) #*A vector containing the rowname of x3 (and x2) of those GOes for which NA is in all tissues

minR<-apply(x3after, 1, min_f)
x4<-cbind(x3after,minR) #1102x36 (AUC_tissue1::minR)



this_tissues_MIN<-c()
mins<-c()
for(j in 1:dim(x4)[1]){
	nextt=F
	x5<-x4[j,]
	x6<-x5[colSums(!is.na(x5)) > 0]
	till=dim(x6)[2]-1
	for(i in 1:till){
		if(x6[,i]==x6$minR){  #This is require in case more than one have euql to min
			this_tissues_MIN<-append(this_tissues_MIN,colnames(x6)[i])
			mins<-append(mins,x6$minR)
			break}
	}
}


length(this_tissues_MIN)# 1102. vector with the tissue with lowes AUC for each GO term

	#1b) find the best tissue (and AUC) for each GO  

max_f<-function(ddff){
	u<-max(ddff,na.rm = T)
	return(u)}

cc<-dim(x3)[2]
x3after<-x3[rowSums(is.na(x3))!=cc, ]
out2<-as.numeric(rownames(x3)[!rownames(x3) %in% rownames(x3after)]) 

maxR<-apply(x3after, 1, max_f)
x4<-cbind(x3after,maxR)



this_tissues_MAX<-c()
maxes<-c()
for(j in 1:dim(x4)[1]){
	nextt=F
	x5<-x4[j,]
	x6<-x5[colSums(!is.na(x5)) > 0]
	till=dim(x6)[2]-1
	for(i in 1:till){
		if(x6[,i]==x6$maxR){
			this_tissues_MAX<-append(this_tissues_MAX,colnames(x6)[i])
			maxes<-append(maxes,x6$maxR)
			break}
	}
}


length(this_tissues_MAX)# 1102. vector with the tissue with highest AUC for each GO term


	#1c) extract the name(s) so that dimG and  length(this_tissues_MAX) match (with in max and min)

G<-x2$go_number #1103 GO_IDs in same order as the vectors "this_tissues_MIN" and "this_tissues_MAX". 1103>1102, so exclude (above*)
out_goes<-c()
for(i in 1:length(out)){	
	out_go<-as.character(x2$go_number[rownames(x2)==as.character(out)[i]])
	out_goes<-append(out_goes,out_go)}
for(i in 1:length(out2)){	
	out_go<-as.character(x2$go_number[rownames(x2)==as.character(out2)[i]])
	out_goes<-append(out_goes,out_go)}

G<-G[!G %in% out_goes]
G<-data.frame(G)


	#1d) give proper name to tissue
add_all<-c()
for(i in 1:35){
	add<-paste("AUC_tissue",i,sep="")
	add_all<-append(add_all,add)}
x$V1<-add_all #35x2 (i.e. AUC_tissue1  Stomach...)


hs<-c()
for(i in 1:length(this_tissues_MIN)){	
	for(j in 1:length(x$V1)){
		if(this_tissues_MIN[i]==x$V1[j]){
			h=as.character(x$V2[j])
			hs<-append(hs,h)}
	}
}
hs<-data.frame(hs)
colnames(hs)<-c("tissue_with_min_AUC")
G<-cbind(G,hs)


hs<-c()
for(i in 1:length(this_tissues_MAX)){	
	for(j in 1:length(x$V1)){
		if(this_tissues_MAX[i]==x$V1[j]){
			h=as.character(x$V2[j])
			hs<-append(hs,h)}
	}
}
hs<-data.frame(hs)
colnames(hs)<-c("tissue_with_max_AUC")
G<-cbind(G,hs)

maxes<-data.frame(maxes)
mins<-data.frame(mins)

G<-cbind(G,maxes)
G<-cbind(G,mins)

G<-G[,c(1,3,4,2,5)]# 1102x5 (i.e. GO:0000003  Pancreas	0.645	Skin-Not_Sun_Exposed(Suprapubic)	0.617)



	#1e) write table2.txt
		#after 2)



#2) generate "info_from_table2_to1.txt": mean, min, max (and sd) of each GO term across tissues and GO term description

	#2a) mean (and sd) of each GO term across tissues
x2<-x2[!x2$go_number %in% out_goes,]
means<-c()
sds<-c()
mean_all<-c()
for(i in 1:dim(x2)[1]){
	x2_go<-x2[i,1]
	x2_go_aucs<-x2[i,2:36]
	x2_go_aucs <- x2_go_aucs[!is.na(x2_go_aucs)]
	sds<-append(sds,sd(x2_go_aucs))
	mean_all<-append(mean_all,mean(x2_go_aucs))
	med<-median(as.numeric(x2_go_aucs))
	above<-x2_go_aucs[x2_go_aucs>=med]
	below<-x2_go_aucs[x2_go_aucs<med]
	dif<-mean(above)-mean(below)
	means<-append(means,dif)
}

	#2b)get name of GO
library(GO.db)    #google: GO.db GOTERM
Fer_GOTERM<-as.character(x2$go_number) #1102. (i.e. "GO:0000003"...)

Terms<-c()  #1102. (i.e. "reproduction"...)
for(i in 1:length(Fer_GOTERM)){	
	f<-Term(Fer_GOTERM[i])
	f_n<-as.vector(unlist(f))
	if(is.null(f_n)){print(i)}
	Terms<-append(Terms,f_n[1])
}

	#2c)write file
G<-cbind(G,Terms)
G$dif_max_min<-G$maxes-G$mins
difference<-G$maxes-G$mins

G<-G[,c(1,6,3,2,5,4,7)]
#the info form mean and sd across tissues is to be added. Remove tissue_with_max_AUC...
info<-cbind(G,mean_all)
info<-cbind(info,sds)
info<-cbind(info,means)
colnames(info)[1]<-"go_number"
colnames(info)[2]<-"term"
colnames(info)[3]<-"AUC_best_prdictd_tis"
info$tissue_with_max_AUC<-NULL
info$tissue_with_min_AUC<-NULL
colnames(info)[4]<-"AUC_worst_prdictd_tis"
colnames(info)[6]<-"mean_across_tissues"
colnames(info)[7]<-"sd_across_tissues"
colnames(info)[8]<-"difAUC_50best50worst_tis"
info<-info[with(info, order(-difference_AUC_50Pbest_50Pworst_tissues)), ]#sorted from more to les difference
info<-info[,c(1,3,2,4,5,6,7,8)]
info<-info[,-c(3,4,5)]#comment this if table should contain info on best and worst predicted tissue-AUC
info<-cbind(info,Terms)
info$Terms<-gsub(" ","_",info$Terms)

write.table(info, file="/home/bueno002/table2/info_from_table2_to1.txt", col.names = T, row.names = F, quote = F, sep="\t")



#1) Generate "table2.txt" [continue]

	#1e) write table2.txt
mas<-c()
mis<-c()
for(i in 1:dim(G)[1]){
	ma<-paste(G$tissue_with_max_AUC[i],G$maxes[i])
	mas<-append(mas,ma)
	mi<-paste(G$tissue_with_min_AUC[i],G$mins[i])
	mis<-append(mis,mi)}

G<-x2$go_number
G<-G[!G %in% out_goes]
G<-data.frame(G)
G<-cbind(G,mas)
G<-cbind(G,mis)
G<-cbind(G,Terms)
G$G<-NULL
G<-G[,c(3,1,2)]
colnames(G)<-c("go_term","max_AUC","min_AUC")
difference<-sort(difference,decreasing = T)
G$dif_max_min<-difference
write.table(G, file="/home/bueno002/table2/table2.txt", col.names = T, row.names = F, quote = F, sep="\t")


#3) generate "table2_35rows.txt": Top 3 GOES for each tissue

x2<-read.table("/home/bueno002/table1/r_tissues/AUCtissues_perGO.txt",header=T)
x2<-x2[!x2$go_number %in% out_goes,]
x2$term<-Terms
x2$go_number<-NULL

tops1<-c()
tops2<-c()
tops3<-c()
for(K in 1:35){

	x2<-read.table("/home/bueno002/table1/r_tissues/AUCtissues_perGO.txt",header=T)#These 5 lines are required because x2 is redefined inside the loop
	x2<-x2[!x2$go_number %in% out_goes,]
	x2$term<-Terms
	x2$go_number<-NULL
	x3<-x2[,1:35]

	xx<-x3[,-K]
	MEAN<-rowMeans(xx,na.rm = T)
	FF<-x2[,K]
	diff<-(MEAN-FF)/MEAN	#GOs that differ a lot bt tissue 1 and the mean of the rest, correcting by the mean
	diff_df<-data.frame(diff)
	diff_df<-cbind(diff_df,x2$term)
	colnames(diff_df)[2]<-"term"
	x2or<-diff_df[with(diff_df, order(-diff)), ]
	top1<-x2or$term[1]
	tops1<-append(tops1,top1)
	top2<-x2or$term[2]
	tops2<-append(tops2,top2)
	top3<-x2or$term[3]
	tops3<-append(tops3,top3)
}

#
Tr<-data.frame(x$V2)
Tr<-cbind(Tr,tops1)
Tr<-cbind(Tr,tops2)
Tr<-cbind(Tr,tops3)
colnames(Tr)<-c("tissue","top1_GOterm","top2_GOterm","top3_GOterm")
write.table(Tr, file="/home/bueno002/table2/table2_35rows.txt", col.names = T, row.names = F, quote = F, sep="\t")








#4) generate "difference_less_corr_tissues"
x2<-read.table("/home/bueno002/table1/r_tissues/AUCtissues_perGO.txt",header=T) #1103x36 (go_number::AUC_tissue35)
x2<-x2[!x2$go_number %in% out_goes,]
x2<-x2[,2:36]
H<-colMeans(x2,na.rm = T)

x<-read.table("/home/bueno002/table2/tissue_names") #2x2 (i.e. 1  Stomach; 2	Esophagus-Muscularis...)

x2<-x2[complete.cases(x2), ]

#colnames(x2)<-x$V2


library(dplyr)
library(reshape2)
d_cor <- as.matrix(cor(x2))
d_cor_melt <- arrange(melt(d_cor), -abs(value))
d_cor_melt<-d_cor_melt[-c(1:35),]
S<-seq(1,dim(d_cor_melt)[1],2)#remove duplicates
d_cor_melt2<-d_cor_melt[-S,]







write.table(d_cor_melt2, file="/home/bueno002/table2/corr_bt_tissues.txt", col.names = T, row.names = F, quote = F, sep="\t")


t<-(x2$AUC_tissue33-x2$AUC_tissue35) 
jpeg('difference_less_corr_tissues.jpg')
plot(density(t), main="difference_less_corr_tissues", xlab="differenceAUC", ylab="density GO terms", xlim=c(0, 0.08),ylim=c(0, 42), col="red")
dev.off()





