# Takes a tab delimited files with annotations and up-propagates
uppropagate = function(FILE)
{


	library(GO.db);
	library(igraph);	
	BP = toTable(GOBPPARENTS);
	MF = toTable(GOMFPARENTS);
	CC = toTable(GOCCPARENTS);
	
	ancBP = as.list(GOBPANCESTOR);
	ancMF = as.list(GOMFANCESTOR);
	ancCC = as.list(GOCCANCESTOR);
	
	outann = matrix(0,nrow=0, ncol = 3);
	
	inann = as.matrix(read.table(FILE, sep="\t"));


	proteins = unique(as.vector(inann[,1]));

	plist = list();	
	
	for(r in 1:nrow(inann))
	{
		ptemp = inann[r,1];
		gtemp = inann[r,2];
		
		if (inann[r,3] == ' F' || inann[r,3] == 'F')
		{
			plist[[ptemp]] = c(plist[[ptemp]], gtemp, ancMF[[gtemp]])
		}

		if (inann[r,3] == ' P' || inann[r,3] == 'P')
		{
			plist[[ptemp]] = c(plist[[ptemp]], gtemp, ancBP[[gtemp]])
		}
		
		if (inann[r,3] == ' C' || inann[r,3] == 'C')
		{
			plist[[ptemp]] = c(plist[[ptemp]], gtemp, ancCC[[gtemp]])
		}
		
		plist[[ptemp]] = unique(plist[[ptemp]]);
	}

	outanntemp = matrix(0,nrow=0, ncol=2);
	for (p in 1:length(proteins))
	{

		ptemp = proteins[p];

		ptempv = rep(ptemp, length(plist[[ptemp]]));
		ctemp = cbind(ptempv, plist[[ptemp]]);

		outanntemp = rbind(outanntemp,ctemp);
	}
		
	

	ontology = vector(mode="character", length=nrow(outanntemp));
	
        print(nrow(outanntemp));
	for (r in 1:nrow(outanntemp))
	{

		gtemp = as.vector(outanntemp[r,2]);

		if(gtemp == "all" || gtemp == " all") { next;}

		if(is.element(gtemp, MF$go_id))
		{
			ontology[r] = 'F';
		}
		if(is.element(gtemp, BP$go_id))
		{
			ontology[r] = 'P';
		}
		if(is.element(gtemp, CC$go_id))
		{
			ontology[r] = 'C';
		}
		
		outann = cbind(outanntemp, ontology);
		
	}	
	outfile = paste(FILE,".uppropagated", sep="");
	write.table(outann, outfile, quote=F, row.names=F, col.names=F, sep = "\t")

		
		
		


}


