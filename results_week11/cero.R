#qsub -q all.q -t 1-7***** one.sh		******To choose array size
#data: /home/fern/connect/tables_monday/data/
args<-commandArgs(T)
library(Matrix);
library(methods) #This is required when running Rscript
library(plyr)
source("/home/bueno002/tables_monday/bmrf_functions.R");
source("/home/bueno002/tables_monday/validation_functions.R")	


k=10
no_R=20
que=rep(1:k,no_R)
noit_GS<-30

minGOsize=20
maxGOsize=0.1
only_EES_BP=F
subset=F							#Choose [1,2,3,4,5,F], depending on folders: "From_gitHub/large_coex/additional_inputs_and_plots/fileS/" 
									#To see size and description of subsets: http://www.inetbio.org/yeastnet/downloadnetwork.php	
Ps<-c(0.1,0.05,0.03,0.01,0.007)		#Portions to substarct. If subsets != F, Ps will take value 0.	
Ps<-c(0)							#For conventional run: Ps<-c(0) and subset=F

all=T
			


if(all){
	go<-as.numeric(args[1])	
} else {
	goes<-c("GO:0006417","GO:0031670","GO:0006414","GO:0051054","GO:0045931","GO:0007533","GO:0000209")
	go<-as.numeric(args[1])
	go_number<-as.character(goes[go]) 
}			
		


#Call functions
options(warn=2)
reps_all<-c()
for(P in 1:length(Ps)){
	loaded<-load_L_m(P=Ps[P],subset=subset,minGOsize=minGOsize,maxGOsize=maxGOsize,only_EES_BP=only_EES_BP)
	print(dim(loaded$L_m)[2])	#******To choose array size
	ws<-c()
	for(i in 1:length(que)){
		w<-try(AUCf(que[i],k,all=all), silent=FALSE)
		if(class(w) == 'try-error')
		{
		  w<-NA
		}
		ws<-append(ws,w)
	}
	reps<-c()
	for(i in 0:(no_R-1)){
		start<-(i*k)+1
		end<-start+k
		reps<-append(reps,mean(ws[start:end],na.rm = TRUE))			#If any fold from any replicate gives error OR WARNING, it will average over the other folds of the replicate
	}
	reps_all<-append(reps_all,mean(reps,na.rm = TRUE))
}



#print details
if(all){
	go_number<-as.character(colnames(loaded$L_m)[go])	
}
total_labels<-loaded$table_EES_BP[names(loaded$table_EES_BP) == as.character(go_number)]
EES_BP_labels<-loaded$table[names(loaded$table) == as.character(go_number)]
	# # connexions
labels<-loaded$L_m[,colnames(loaded$L_m) == as.character(go_number)]
assoc<-labels[labels == 1]
assoc1<-loaded$network$V1[loaded$network$V1 %in% names(assoc)]
assoc2<-loaded$network$V2[loaded$network$V2 %in% names(assoc)]
no_conn<-length(assoc1)+length(assoc2)


if(length(Ps) == 1){
	column_name<-as.character(c("go_number","total_labels","EES_BP_labels","no_connexions","AUC"))
} else {
	column_name<-as.character(c("go_number","total_labels","EES_BP_labels","no_connexions",paste("AUC_substracting_",Ps[1], "%_of_the_network", sep=""),Ps[2:length(Ps)]))
}
column<-as.character(c(go_number,total_labels,abs(EES_BP_labels),no_conn,round(reps_all,3)))


out.file <- paste("/home/bueno002/tables_monday/r/", all, subset, minGOsize, maxGOsize,  only_EES_BP, ".txt", sep="_")
if(file.exists(out.file)){
	DF<-data.frame(column)
	DF<-t(DF)
	write(DF, file=out.file,ncolumns=length(column),append=T)
} else {
	DF<-data.frame(column_name)
	DF<-t(DF)
	DF<-rbind(DF,column)
	write.table(DF, file=out.file, col.names = F, row.names = F, quote = F, sep="\t")
}



