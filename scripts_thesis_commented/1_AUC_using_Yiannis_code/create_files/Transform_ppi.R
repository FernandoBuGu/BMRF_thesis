#Script to transform the files GO_labellings.dat, GO_terms.dat and "ORF_names.dat" in Yiannis code (@meyers /mnt/scratch/dijk097/Fernando/BMRF-R/) to the format used in version 2 of the code (https://github.com/jwbargsten/bmrf).   

    #MSc thesis bioinfomratics WUR. Protein function prediction for poorly annotated species.
    #Author: Fernando Bueno Gutierrez
    #email1: fernando.buenogutierrez@wur.nl
    #email2: fernando.bueno.gutie@gmail.com


mat<-read.table("GO_labellings.dat")
go<-read.table("GO_terms.dat")[,1]
label<-read.table("ORF_names.dat")[,1]
h<-c()
for(i in 1:dim(mat)[2])
{
	for(j in 1:dim(mat)[1])
	{
		if(mat[j,i]==1) {
			o<-paste(go[i],label[j],sep=" ")
			h<-rbind(h,o)
		} else {
			next;
		}
	}
}
write.table(h,"GO_file.tsv", sep= "\t", col.names = F, row.names = F, quote = F)
r<-read.table("../ppi_data/collins.dat")
write.table(r,file="../ppi_data/collins.dat", sep="\t", col.names = F, row.names = F, quote = F)

