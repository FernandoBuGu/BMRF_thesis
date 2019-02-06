#script to run PU-BMRF analysis with tissues data. This R script runs in parallel
#qsub -q all.q -t 1-1887 tissues.sh
args<-commandArgs(T)
old_start <- Sys.time()
library(Matrix);
library(methods) #This is required when running Rscript
library(plyr)
library(foreach)
library(doParallel)
source("/home/WUR/bueno002/tissues/data_and_functions/bmrf_functions.R");
source("/home/WUR/bueno002/tissues/data_and_functions/validation_functions_with_tissue.R")	


#fixed parameters. 
k=7; no_R=1; noit_GS=20


#Parameters to choose
minGOsize=20
maxGOsize=0.05
only_EES_BP=F
subset=F
all=T							#Choose [1,2,3,4,5,F], depending on folders: "From_gitHub/large_coex/additional_inputs_and_plots/fileS/" 
reduce="F"		#amg,oa,epp,epn,enn   associationsof my go; other associations; edges of proteins of my go; other edges	
spec_tissue=T



if(reduce != F){PRs<-c(0.05,0.1,0.2,0.4,0.6,0.8,0.99)}		#Portions to substarct	
							
			
if(spec_tissue){Tissues<-seq(1:35)}


loaded<-return_L_m(subset=subset,minGOsize=minGOsize,maxGOsize=maxGOsize,only_EES_BP=only_EES_BP)
options(warn=0)
print("max_array_size="); print(dim(loaded$L_m)[2])

G<-as.numeric(args[1])
go<-as.character(colnames(loaded$L_m)[G])
print(go)

loaded_and_AUCf<-function(i,go,k,all=all,noit_GS){
	loaded<-return_L_m(subset=subset,minGOsize=minGOsize,maxGOsize=maxGOsize,only_EES_BP=only_EES_BP,tissue=i)
	aucf_result<-AUCf(go,k,all=all,noit_GS)
	return(aucf_result)}


#Calculate AUC for each Tissue "no_R" times
cores=detectCores()
cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)
clusterExport(cl,varlist=ls(),envir=environment())
w<-foreach(i = sort(rep(Tissues,no_R)), .packages='Matrix', .combine="c") %dopar% { loaded_and_AUCf(i,go,k,all=all,noit_GS) }
stopCluster(cl)
#


options(warn=0)
print("max_array_size="); print(dim(loaded$L_m)[2])



#Organize output from cluster
len=length(Tissues)
till=len-1
reps<-c()
for(i in 1:len){
	if(i==1){start=1}
	else{start<-end+1}
	end<-start+no_R-1
	reps<-append(reps,mean(w[start:end],na.rm = TRUE))
}
reps_sd<-c()
for(i in 1:len){
	if(i==1){start=1}
	else{start<-end+1}
	end<-start+no_R-1
	reps_sd<-append(reps_sd,sd(w[start:end],na.rm = TRUE))
}
ps<-c()
for(i in 1:length(reps)){
	p<-paste(round(reps[i],3), round(reps_sd[i],3), sep="\t")
	ps<-append(ps,p)}
header<-c()
for(i in 1:length(Tissues)){
	hm<-Tissues[i]
	hsd<-paste("sd",Tissues[i],sep="_")
	header<-append(header,hm)
	header<-append(header,hsd)}


#write in file
column_name<-as.character(c("go_number",header))
column<-as.character(c(go,ps))

out.file <- paste("/home/WUR/bueno002/tissues/r_tissues/r_tissues_file", ".txt", sep="")
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



new_end <- Sys.time() - old_start
print(new_end)
print("FINSHED")





