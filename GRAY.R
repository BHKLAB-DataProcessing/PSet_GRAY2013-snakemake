library(PharmacoGx)
library(readxl)
library(openxlsx)
library(tximport)
library(Biobase)
library(data.table)
library(reshape2)
library(CoreGx)
library(SummarizedExperiment)

options(stringsAsFactors = FALSE)
myDirPrefix <- "/pfs/"
args = commandArgs(trailingOnly=TRUE)
rnaseq_select <- args
print(rnaseq_select)
rnaseq_results <- list()
ORCESTRA_ID = tail(rnaseq_select, n=1)

	  
tools <- grep(pattern = 'Kallisto|Salmon', x = rnaseq_select)
tools <- rnaseq_select[tools]
tools <- gsub("-", "_", tools)
transcriptome <- grep(pattern = 'Gencode|Ensembl', x = rnaseq_select)
transcriptome <- rnaseq_select[transcriptome]
tool_path = expand.grid(a = tools,b = transcriptome)
tool_path = paste0(tool_path$a, "_",tool_path$b)
	  
print(tool_path)


    
	  
    #match to curations
    
    
    matchToIDTable <- function(ids,tbl, column, returnColumn="unique.cellid") {
      sapply(ids, function(x) {
        myx <- grep(paste0("((///)|^)",Hmisc::escapeRegex(x),"((///)|$)"), tbl[,column])
        if(length(myx) > 1){
          stop("Something went wrong in curating ids, we have multiple matches")
        }
        if(length(myx) == 0){return(NA_character_)}
        return(tbl[myx, returnColumn])
      })
    }
    
    badchars <- "[\xb5]|[]|[ ,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[\\]|[.]|[_]|[ ]"
    
    #get curations
    
    cell_all <- read.csv(file = "/pfs/downAnnotations/cell_annotation_all.csv", na.strings=c("", " ", "NA"))
    curationCell <- cell_all[which(!is.na(cell_all[ , "GRAY.cellid"])),]
    curationCell <- curationCell[ , c("unique.cellid", "GRAY.cellid")]
    rownames(curationCell) <- curationCell[ , "unique.cellid"]
    
    curationTissue <- cell_all[which(!is.na(cell_all[ , "GRAY.cellid"])),]
    curationTissue <- curationTissue[ , c("unique.tissueid", "GRAY.tissueid")]
    rownames(curationTissue) <- curationCell[ , "unique.cellid"]
    
    drug_all <- read.csv("/pfs/downAnnotations/drugs_with_ids.csv", na.strings=c("", " ", "NA"))
    curationDrug <- drug_all[which(!is.na(drug_all[ , "GRAY.drugid"])),]
    curationDrug <- curationDrug[,c("unique.drugid","GRAY.drugid")]
    rownames(curationDrug) <- curationDrug[ , "unique.drugid"]
    
    
    load("/pfs/GRAYRawSensitivity/drug_norm_post.RData")
    
    # cell information (cell slot)
    
    cellineinfo <- read.xlsx("/pfs/getGRAY2013/gb-2013-14-10-r110-s1.xlsx", sheet = 1)
    cellineinfo[!is.na(cellineinfo) & cellineinfo == ""] <- NA
    rn <- cellineinfo[-1, 1]
    cn <- t(cellineinfo[1, -1])
    cn <- gsub(badchars, ".", cn)
    cellineinfo <- cellineinfo[-1, -1]
    dimnames(cellineinfo) <- list(rn, cn)
    cellineinfo <- data.frame("cellid"=rn, "tissueid"="breast", cellineinfo[,1:10])
    cellineinfo <- cellineinfo[which(!is.na(cellineinfo$Transcriptional.subtype)), ]
    
    #duplicate MB157, need to merge with MDAMB157
    cellineinfo[31,"RNASeq.availability"] <- "1"
    cellineinfo[31,"Transcriptional.subtype"] <- "Claudin-low/Basal"
    cellineinfo[31,"ERBB2.status"] <- "Claudin-low/Basal"
    cellineinfo <- cellineinfo[which(!cellineinfo$cellid == "MB157"),]
    
    cellineinfo$cellid <- as.character(matchToIDTable(ids=cellineinfo$cellid, tbl=curationCell, column = "GRAY.cellid", returnColumn = "unique.cellid"))
    rownames(cellineinfo) <-  cellineinfo$cellid
    
    #curationCell <- curationCell[rownames(cellineinfo),]
    
    #published sensitivity
    
    profiles <- read.xlsx("/pfs/getGRAY2013/gb-2013-14-10-r110-s1.xlsx", sheet = 1)
    profiles[!is.na(profiles) & profiles == ""] <- NA
    rn <- profiles[-1, 1]
    cn <- t(profiles[1, -1])
    profiles <- profiles[-1, -1]
    dimnames(profiles) <- list(rn, cn)
    profiles <- profiles[which(!is.na(profiles[, "Transcriptional subtype"])), ]
    colnames(profiles)[which(colnames(profiles) == "L-779450")] <-  "L-779405"
    indices <- 11:ncol(profiles)
    GI50 <- as.numeric(array(apply(profiles, 1, function(x)(x[indices]))))
    
    drugs <- as.character(matchToIDTable(ids=colnames(profiles)[indices], tbl=curationDrug, column = "GRAY.drugid", returnColumn = "unique.drugid"))
    celllines <- as.character(matchToIDTable(ids=rownames(profiles), tbl=curationCell, column = "GRAY.cellid", returnColumn = "unique.cellid"))
    x <- expand.grid(drugs,celllines)
    names(GI50) <- paste("drugid", paste(x[,1],x[,2],sep = "_"), sep="_")
    GI50 <- GI50[which(!is.na(GI50))]
    sensitivity.profiles <- matrix(NA, dimnames = list(rownames(sensitivity.info), "GI50_published"), nrow=nrow(sensitivity.info))
    for(nn in names(GI50)) {
      sensitivity.profiles[grep(nn, rownames(sensitivity.profiles)), "GI50_published"] <- GI50[nn]
    }
    
    
    load("/pfs/gray2013ProfilesAssemble/profiles.RData")
    
    #compile sensitivity profiles
               
    sensitivity.profiles <-  data.frame("aac_recomputed" = as.numeric(res[,"AAC"]), "ic50_recomputed"=as.numeric(res[,"IC50"]), "HS"=as.numeric(res[,"HS"]), "E_inf"=as.numeric(res[,"E_inf"]), "EC50"=as.numeric(res[,"EC50"]), "GI50_published"= as.numeric(sensitivity.profiles[,"GI50_published"]))
    
    sensitivity.profiles$aac_recomputed <- sensitivity.profiles$aac_recomputed/100
    
    #compute slope and add to sensitivity profiles
                                   
    slope <- NULL
    for(exp in rownames(sensitivity.info)){
      slope <- c(slope, computeSlope(raw.sensitivity[exp, , "Dose"], raw.sensitivity[exp, , "Viability"])) #computeSlope (returns normalized slope of drug response curve)
    }
    
    names(slope) <- rownames(sensitivity.info)
    sensitivity.profiles <- cbind(sensitivity.profiles, "slope_recomputed"=slope)
    head(sensitivity.profiles)
                                   
    sensitivity.info <- as.data.frame(sensitivity.info)
    # drug info (drug slot)

    curationDrug <- curationDrug[as.character(unique(sensitivity.info[,"drugid"])),]
    druginfo <- data.frame("drugid"=curationDrug$unique.drugid)
    rownames(druginfo) <- druginfo$drugid
    
    
summarizeRnaSeq <- function (dir, 
                             features_annotation,
                             samples_annotation,
			      method) {
  library(Biobase)
  library(readr)
  library(tximport)
  
  load(features_annotation)
    
  tx2gene <- as.data.frame(cbind("transcript"=tx2gene$transcripts, "gene"=tx2gene$genes))
  
  files <- list.files(dir, recursive = TRUE, full.names = T)
  if(method=="kallisto"){
  resFiles <- grep("abundance.h5", files)
  }else{
  resFiles <- grep("quant.sf", files)
  }
  resFiles <- files[resFiles]
  length(resFiles)
  names(resFiles) <- basename(dirname(resFiles))
  
  if(features_annotation == "/pfs/downAnnotations/Ensembl.v99.annotation.RData"){
  txi <- tximport(resFiles, type=method, tx2gene=tx2gene, ignoreAfterBar = TRUE, ignoreTxVersion = TRUE)
  } else{
  txi <- tximport(resFiles, type=method, tx2gene=tx2gene, ignoreAfterBar = TRUE, ignoreTxVersion = FALSE)	  
  }
	  
  head(txi$counts[,1:5])
  dim(txi$counts)
	  
  xx <- txi$abundance
  gene.exp <- Biobase::ExpressionSet(log2(xx + 0.001))
  fData(gene.exp) <- features_gene[featureNames(gene.exp),]
  pData(gene.exp) <- samples_annotation[sampleNames(gene.exp),]
  annotation(gene.exp) <- "rnaseq"
  
  xx <- txi$counts
  gene.count <- Biobase::ExpressionSet(log2(xx + 1))
  fData(gene.count) <- features_gene[featureNames(gene.count),]
  pData(gene.count) <- samples_annotation[sampleNames(gene.count),]
  annotation(gene.count) <- "rnaseq"
  
  txii <- tximport(resFiles, type=method, txOut=T)
  
  if(features_annotation == "/pfs/downAnnotations/Ensembl.v99.annotation.RData"){
  #remove non-coding transcripts in ensembl 	  
  rownames(txii$abundance) <-  gsub("\\..*","",rownames(txii$abundance))
  txii$abundance[which(!rownames(txii$abundance)  %in% features_transcript$transcript_id)]
  missing_transcript <- rownames(txii$abundance)[which(!rownames(txii$abundance)  %in% features_transcript$transcript_id)]
  txii$abundance <- txii$abundance [-which(rownames(txii$abundance) %in% missing_transcript),]
  }
  	  
  xx <- txii$abundance
  transcript.exp <- Biobase::ExpressionSet(log2(xx[,1:length(resFiles)] + 0.001))
  if(features_annotation == "/pfs/downAnnotations/Gencode.v33.annotation.RData" || features_annotation == "/pfs/downAnnotations/Gencode.v33lift37.annotation.RData"){
  featureNames(transcript.exp) <- gsub("\\|.*","",featureNames(transcript.exp))
  fData(transcript.exp) <- features_transcript[featureNames(transcript.exp),]
  }else{
  fData(transcript.exp) <- features_transcript[featureNames(transcript.exp),]
  }
  pData(transcript.exp) <- samples_annotation[sampleNames(transcript.exp),]
  annotation(transcript.exp) <- "isoform"
  
	  
  if(features_annotation == "/pfs/downAnnotations/Ensembl.v99.annotation.RData"){
  #remove non-coding transcripts in ensembl
  rownames(txii$counts) <-  gsub("\\..*","",rownames(txii$counts))
  txii$counts <- txii$counts [-which(rownames(txii$counts) %in% missing_transcript),]	  
  }	  
  xx <- txii$counts
  transcript.count <- Biobase::ExpressionSet(log2(xx[,1:length(resFiles)] + 1))
  if(features_annotation == "/pfs/downAnnotations/Gencode.v33.annotation.RData" || features_annotation == "/pfs/downAnnotations/Gencode.v33lift37.annotation.RData"){
  featureNames(transcript.count) <- gsub("\\|.*","",featureNames(transcript.count))
  fData(transcript.count) <- features_transcript[featureNames(transcript.count),]
  }else{
  fData(transcript.count) <- features_transcript[featureNames(transcript.count),]
  }
  pData(transcript.count) <- samples_annotation[sampleNames(transcript.count),]
  annotation(transcript.count) <- "isoform"
  
  return(list("rnaseq"=gene.exp, 
              "rnaseq.counts"=gene.count, 
              "isoforms"=transcript.exp, 
              "isoforms.counts"=transcript.count))
}

    

    rnaseq.sampleinfo <- read.csv(file="/pfs/downAnnotations/JRGraySRRMapping.csv", stringsAsFactors=FALSE, row.names=1)
    
    rnaseq.sampleinfo[ , "cellid"] <- as.character(matchToIDTable(ids=rnaseq.sampleinfo[ , "cellid"], tbl=curationCell, column = "GRAY.cellid", returnColumn = "unique.cellid"))
    rnaseq.sampleinfo[ , "batchid"] <- NA
				   
    for (r in 1:length(tool_path)){
  print(tool_path[r])
  if (length(grep(pattern = 'Kallisto', x = tool_path[r])) > 0){
    tool <- sub("(_[^_]+)_.*", "\\1", tool_path[r])
    tdir = paste0("gray_rnaseq_",gsub(".","_",tolower(tool), fixed = T), "/",  tool, "/", tool, "/")  
    rnatool="kallisto"	  
  } else {
    tool <- sub("(_[^_]+)_.*", "\\1", tool_path[r])
    tdir = paste0("gray_rnaseq_",gsub(".","_",tolower(tool), fixed = T), "/",  tool, "/", tool, "/")
    rnatool="salmon"	  
  }
  
  
  if (length(grep(pattern = 'lift37', x = tool_path[r])) > 0){
    annot = "/pfs/downAnnotations/Gencode.v33lift37.annotation.RData"
  } else if (length(grep(pattern = 'v33', x = tool_path[r])) > 0){
    annot = "/pfs/downAnnotations/Gencode.v33.annotation.RData"
  } else {
    annot = "/pfs/downAnnotations/Ensembl.v99.annotation.RData"
  }
    print(annot)
  
  print(tdir)
  rnaseq <- summarizeRnaSeq(dir=file.path(paste0(myDirPrefix, tdir, tool_path[r])),
                            features_annotation=annot,
                            samples_annotation=rnaseq.sampleinfo,
			    method = rnatool)
  rnaseq_results <- c(rnaseq_results,c(
    rnaseq <- setNames(rnaseq,  paste0(tool,".", names(rnaseq)))
  )
  )
}

                                   
  rnaseq_cellid_all <- pData(rnaseq_results[[1]])[,"cellid"]  
  #add missing cells to cell info                                 
  cellnall <- CoreGx::.unionList(rownames(cellineinfo),rnaseq_cellid_all, sensitivity.info$cellid)
  newcells <- setdiff(cellnall, rownames(cellineinfo))
  newRows <- matrix(NA_character_, nrow=length(newcells), ncol=ncol(cellineinfo))
  # newRows <- cell.info[newcells,]

  rownames(newRows) <- newcells
  colnames(newRows) <- colnames(cellineinfo)
  newRows[,"cellid"] <- newcells

  cellineinfo <- rbind(cellineinfo, newRows)
    
    
  cellsPresent <- sort(CoreGx::.unionList(sensitivity.info$cellid,rnaseq_cellid_all))
  cellineinfo <- cellineinfo[cellsPresent,]                             

  cellineinfo$tissueid <- curationTissue[rownames(cellineinfo), "unique.tissueid"]
  cellineinfo$cellid <- rownames(cellineinfo)
  
  curationTissue <- curationTissue[rownames(cellineinfo),]
  curationCell <- curationCell[rownames(cellineinfo),]
				   
  drug_all <- read.csv("/pfs/downAnnotations/drugs_with_ids.csv", na.strings=c("", " ", "NA"))
  drug_all <- drug_all[which(!is.na(drug_all[ , "GRAY.drugid"])),]
  drug_all <- drug_all[ , c("unique.drugid", "GRAY.drugid","smiles","inchikey","cid","FDA")]
  rownames(drug_all) <- drug_all[ , "unique.drugid"]

  drug_all <- drug_all[rownames(druginfo),]
  druginfo[,c("smiles","inchikey","cid","FDA")] <- drug_all[,c("smiles","inchikey","cid","FDA")] 
 	 
z <- list()

z <- c(z,c(
  rnaseq_results
  )
)
				   
		 
	 
#add cellosaurus disease type to cell-info
disease <- cell_all$Cellosaurus.Disease.Type[match(cellineinfo$cellid, cell_all$unique.cellid)]
cellineinfo$Cellosaurus.Disease.Type <- disease
	 
#add cellosaurus assession to cell-info
assession <- cell_all$Cellosaurus.Accession.id[match(cellineinfo$cellid, cell_all$unique.cellid)]
cellineinfo$Cellosaurus.Accession.id <- assession
		 
#add pharmacodb id to cell-info
pdb <- cell_all$PharmacoDB.id[match(cellineinfo$cellid, cell_all$cellid)]
cellineinfo$PharmacoDB.id <- pdb

#add study tissue id to cell_info
study_tissue <- cell_all$unique.tissueid.fromstudies[match(cellineinfo$cellid, cell_all$unique.cellid)]
cellineinfo$unique.tissueid.fromstudies <- study_tissue
		 
#add study cell-line type to cell_info
cell_type <- cell_all$CellLine.Type[match(cellineinfo$cellid, cell_all$unique.cellid)]
cellineinfo$CellLine.Type <- cell_type
		 
#add metastatic info to cell_info		 
metastatic <- cell_all$Metastatic[match(cellineinfo$cellid, cell_all$unique.cellid)]
cellineinfo$Metastatic <- metastatic	
		 
.converteSetToSE <- function(eSets) {
  
  SEfinal <- lapply(eSets,
         function(eSet){
             # Change rownames from probes to EnsemblGeneId for rna data type
             if (grepl("^rna$", Biobase::annotation(eSet))) {
               rownames(eSet) <- Biobase::fData(eSet)$EnsemblGeneId
             }
             
             # Build summarized experiment from eSet
             SE <- SummarizedExperiment::SummarizedExperiment(
               ## TODO:: Do we want to pass an environment for better memory efficiency?
               assays=S4Vectors::SimpleList(as.list(Biobase::assayData(eSet))
               ),
               # Switch rearrange columns so that IDs are first, probes second
               rowData=S4Vectors::DataFrame(Biobase::fData(eSet),
                                            rownames=rownames(Biobase::fData(eSet)) 
               ),
               colData=S4Vectors::DataFrame(Biobase::pData(eSet),
                                            rownames=rownames(Biobase::pData(eSet))
               ),
               metadata=list("experimentData" = eSet@experimentData, 
                             "annotation" = Biobase::annotation(eSet), 
                             "protocolData" = Biobase::protocolData(eSet)
               )
             )
             ## TODO:: Determine if this can be done in the SE constructor?
             # Extract names from expression set
             SummarizedExperiment::assayNames(SE) <- Biobase::assayDataElementNames(eSet)
             mDataType <- Biobase::annotation(eSet)
             eSets[[mDataType]] <- SE
         })
  #setNames(pSet@molecularProfiles, names(eSets))
  return(SEfinal)
}
		 
z <- .converteSetToSE(z)		 
		 
	 
		 				   

    GRAY2013 <- PharmacoGx::PharmacoSet(molecularProfiles=z,
                            name="GRAY", 
                            cell=cellineinfo, 
                            drug=druginfo, 
                            sensitivityInfo=sensitivity.info, 
                            sensitivityRaw=raw.sensitivity, 
                            sensitivityProfiles=sensitivity.profiles, 
                            sensitivityN=NULL,
                            curationCell=curationCell, 
                            curationDrug=curationDrug, 
                            curationTissue=curationTissue, 
                            datasetType="sensitivity")
    

GRAY2013@annotation$version <- 2
saveRDS(GRAY2013,file="/pfs/out/GRAY2013.rds")				   
dataset <- "GRAY"
#output ORCESTRA_ID and Pachyderm commit id
write.table(dataset, file="/pfs/out/dataset.txt", row.names = F ,quote = F, sep = "\t", col.names = F)	
write.table(ORCESTRA_ID, file="/pfs/out/orcestra_id.txt", row.names = F ,quote = F, sep = "\t", col.names = F)				   
pach_commit_id <- Sys.getenv("PACH_OUTPUT_COMMIT_ID")
write.table(pach_commit_id, file="/pfs/out/commit_id.txt", row.names = F ,quote = F, sep = "\t", col.names = F) 				   
