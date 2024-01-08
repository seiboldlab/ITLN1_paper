


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









####### Pathway enrichment function
doEnrichOneAtATime<-function(DEG_table,dataset,ngenes=NULL,Enrichr_dir="Enrichr",
	enrichr.libraries=c("GO_Cellular_Component_2023","GO_Biological_Process_2023","GO_Molecular_Function_2023",
	"KEGG_2021_Human","Reactome_2022","BioPlanet_2019"),minOverlap=3,minAdjPval=0.2){
	
	masterList<-list()
	for(i in 1:length(unique(DEG_table$comparison))){
		#Set current set of genes
		currData<-as.character(DEG_table$gene[which(DEG_table$comparison == unique(DEG_table$comparison)[i])])
		#Cull genes if desired
		if(!is.null(ngenes)){
			currData<-currData[1:ngenes]
		}
		#Run enrichr
		enriched<-enrichr(currData,enrichr.libraries)
		#Place results from different libraries into a single table
		enriched_tab<-data.frame()
		for(j in 1:length(enriched)){
			if(nrow(enriched[[j]] > 0)){
				enriched_tab<-rbind(enriched_tab,cbind("Gene.Set"=names(enriched)[j],enriched[[j]]))
			}
		}
		#Prepare table for export
		if(nrow(enriched_tab) > 0){
			enriched_tab<-enriched_tab[,-c(6:8)] #weed out unwanted columns
			enriched_tab$Overlap<-gsub("/","_",enriched_tab$Overlap) #Separate slashes with underscores
			enriched_tab$Genes<-gsub(";",", ",enriched_tab$Genes) #Separate genes with commas
			enriched_tab$Genes<-sapply(strsplit(enriched_tab$Genes,", "),function(x)paste(sort(x),collapse=", ")) #sort genes alphabetically
			enriched_tab<-enriched_tab[order(enriched_tab$Adjusted.P.value),] #sort by padj
			enriched_tab<-enriched_tab[which(enriched_tab$Adjusted.P.value < minAdjPval),] #set min padj
			enriched_tab<-enriched_tab[which(sapply(strsplit(enriched_tab$Genes,","), #set min overlap
				function(x)length(x)) >= minOverlap),]
		}
		#If there are zero columns, then add in empty columns
		if(ncol(enriched_tab) == 0){
			enriched_tab<-data.frame(matrix(ncol=7,nrow=0))
			colnames(enriched_tab)<-c("Gene.Set","Term","Overlap","P.value","Adjusted.P.value","Combined.Score","Genes")
		}
		#Add table to master list
		masterList[[i]]<-enriched_tab
		names(masterList)[i]<-as.character(unique(DEG_table$comparison)[i])
	
		#Add slight pause since I'm getting issues with the previous result copying over to the next
		Sys.sleep(3) 
			
	}
	
	#Export master table to file
	wb<-createWorkbook()
	for(i in 1:length(masterList)){
		addWorksheet(wb,sheetName=names(masterList)[i])
		writeData(wb,sheet=names(masterList)[i],masterList[[i]])
		widths<-c(25,80,8,8,12,15,100)
	 	setColWidths(wb, sheet = i, cols = 1:ncol(masterList[[i]]), widths = widths) #or use auto
	 	freezePane(wb,i,firstRow=T)
	 	saveWorkbook(wb, paste(Enrichr_dir,"/",dataset,"_Enrichr.output.xlsx",sep=""), overwrite = TRUE)  
	 }
}














