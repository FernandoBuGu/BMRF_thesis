library(plyr)
library(parallel)
library(pbapply)

source("AUC_function.R");

resu<-function(aa=aa,bb=bb, cc=cc, dd=dd, ee=ee, ff=ff, gg=gg, hh=hh)
{
	source("AUC_function.R");
	library(pbapply)
	w<-pbreplicate(10, AUC_func(aa=aa,bb=bb, cc=cc, dd=dd, ee=ee, ff=ff, gg=gg, hh=hh))
}
j<-as.data.frame(c(20,20,0.9,0.9,T,T,0.2,34343411))
colnames(j)[1]<-"zzero"
j$one<-c(15,20,0.9,0.9,T,T,0.2,34343412)
j$two<-c(25,20,0.9,0.9,T,T,0.2,34343413)
j$thr<-c(5,20,0.9,0.9,T,T,0.2,34343414)
j$fou<-c(20,20,0.9,0.9,F,F,0.2,34343421)
j$fif<-c(15,20,0.9,0.9,F,F,0.2,34343422)
j$six<-c(25,20,0.9,0.9,F,F,0.2,34343423)
j$sev<-c(5,20,0.9,0.9,F,F,0.2,34343424)
j$sv2<-c(20,20,0.9,0.9,T,T,0.5,34343431)
j$eig<-c(15,20,0.9,0.9,T,T,0.5,34343432)
j$nin<-c(25,20,0.9,0.9,T,T,0.5,34343433)
j$ten<-c(5,20,0.9,0.9,T,T,0.5,34343434)
j$ele<-c(30,20,0.9,0.9,T,T,0.2,34343441)
j$twe<-c(20,20,0.9,0.9,F,T,0.2,34343451)
j$the<-c(20,20,0.9,0.9,T,F,0.2,34343461)

cluster <- makeCluster(detectCores())
clusterEvalQ(cluster, library(xts))
result <- clusterMap(cluster, resu, aa=j[1,],bb=j[2,], cc=j[3,], dd=j[4,], ee=j[5,], ff=j[6,], gg=j[7,], hh=j[8,])
save(j,file="j.Rdata")
save(result,file="re.Rdata")


#parse Output
df <- do.call('rbind', result)
u<-rownames(df)
rownames(df)<-NULL
dfg<-apply(df, 2, unlist)
DF<-as.data.frame(dfg)
DF$means<-rowMeans(DF)
DF<-round(DF[,11],2)
DF<-matrix(unlist(t(DF)), byrow=F, 8, 15)
rownames(DF)<-u[1:8]
DF<-as.data.frame(t(DF))
v<-c()
for(i in 1:dim(DF)[1])
{
	v<-append(v,paste(j[,i][1],j[,i][5],j[,i][6],j[,i][7],"subset_x", sep = '-'))
}
rownames(DF)<-v
DF<-DF[ order(-DF[,1]), ]
write.table(DF,"file",sep="\t", col.names = T, row.names = T)

















