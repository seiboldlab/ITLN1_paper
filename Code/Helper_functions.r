


#####Read in excel tabs
readExcelTabs<-function(xlsxFile,returnSingleTable=F){
sheetNames<-openxlsx::getSheetNames(xlsxFile)
sheetList<-lapply(sheetNames,openxlsx::read.xlsx,xlsxFile=xlsxFile)
names(sheetList)<-sheetNames
  if(returnSingleTable){
  	sheetDF<-data.frame()
  	for(sn in sheetNames){
  		if(nrow(sheetList[[sn]]) > 0){
  			sheetDF<-rbind(sheetDF,cbind(sheetList[[sn]],"comparison"= sn))
  		}
  	}
  	return(sheetDF)
  }else{
  	return(sheetList)
  }
}







####### Read in excel tabs more generically
readExcelTabsGeneric<-function(xlsxFile) {
  sheet_names = openxlsx::getSheetNames(xlsxFile)
  sheet_list = lapply(sheet_names, function(x)openxlsx::read.xlsx(xlsxFile, sheet=x))
  names(sheet_list)<-sheet_names
 return(sheet_list)
}






#### Perform hypergeometric/Fisher test of enrichment, where the gene2module table contains
#a set of reference gene lists and where the "genes" vector contains the sample gene list to 
#compare against the reference genes. 
DoEnrichment_alt<-function(genes,gene2module,background=20441){

	moduleColors<-gene2module$Module
	
	#Create vectors to store results in for each module
	enrichment_hyp_p<-data.frame("pval"=rep(NA, length(unique(moduleColors))),
		"module"=rep(1, length(unique(moduleColors))),"genes"=rep(NA,length(unique(moduleColors))),
		"N_sample"=rep(NA,length(unique(moduleColors))),"N_reference"=rep(NA,length(unique(moduleColors))),
		"N_overlap"=rep(NA,length(unique(moduleColors)))) #store hypergeometric pvalues

	##Run loop through each module
	for(i in 1:length(unique(moduleColors))){
		
		#Relevant variables
		A=background #the total population size of genes.
		B=length(genes)# number of genes sampled
		C=length(gene2module$Gene[gene2module$Module==unique(moduleColors)[i]]) #total number of genes in the target or reference list
		D=length(intersect(gene2module$Gene[gene2module$Module==unique(moduleColors)[i]],genes)) #number of sampled genes that overlap 
		                 #with the reference list
		
		####Hypergeometric Test
		targetsInSample = D #The number of target genes in the sample
		targetsInBackground = C #The number of target genes overall
		backgroundMinusTargets = A-targetsInBackground #The number of genes in the total population/background NOT in the target gene list
		sampleSize = B #The size of the sample
		
		#Do the hypergeometric test
		enrichment_hyp_p[i,1]<-phyper((targetsInSample - 1),targetsInBackground,backgroundMinusTargets,
			sampleSize,lower.tail = FALSE)
		
		#Add in module and overlap info
		enrichment_hyp_p[i,2]<-unique(moduleColors)[i]
		#Add in the overlapping genes (these are sorted by marker p-value)
		enrichment_hyp_p[i,3]<-paste(intersect(gene2module$Gene[gene2module$Module==unique(moduleColors)[i]],genes),collapse=", ")
		#Add in N sample
		enrichment_hyp_p[i,4]<-B
		#Add in N reference
		enrichment_hyp_p[i,5]<-C
		#Add in N overlap
		enrichment_hyp_p[i,6]<-D
	}	
	
	#Calculate FDR
	enrichment_hyp_p<-cbind(enrichment_hyp_p,"padj"=as.numeric(p.adjust(enrichment_hyp_p$pval,method="fdr")))
	enrichment_hyp_p<-enrichment_hyp_p[,c(2,4,5,6,1,7,3)]
		
	return(enrichment_hyp_p[order(enrichment_hyp_p$padj),])
}

