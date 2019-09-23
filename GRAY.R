library(PharmacoGx)
library(readxl)
library(openxlsx)
library(tximport)
library(Biobase)

getGRAYP <-
  function (
    verbose=FALSE,
    nthread=1){
    
    options(stringsAsFactors=FALSE)
    
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
    
    curationCell <- curationCell[rownames(cellineinfo),]
    
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
                                   
    
    # drug info (drug slot)

    curationDrug <- curationDrug[as.character(unique(sensitivity.info[,"drugid"])),]
    druginfo <- data.frame("drugid"=curationDrug$unique.drugid)
    rownames(druginfo) <- druginfo$drugid
    
    
    #summarize rnaseq quantifications into expression sets (Kallisto)
                                   
    summarizeRnaSeq <- function (dir, 
                                 tool=c("kallisto", "stringtie", "cufflinks", "rsem", "salmon"), 
                                 features_annotation,
                                 samples_annotation) {
      library(Biobase)
      library(readr)
      library(tximport)
      
      load(features_annotation)
      tx2gene <- as.data.frame(cbind("transcript"=toil.transcripts$transcript_id, "gene"=toil.transcripts$gene_id))
      
      files <- list.files(dir, recursive = TRUE, full.names = T)
      resFiles <- grep("abundance.h5", files)
      resFiles <- files[resFiles]
      length(resFiles)
      names(resFiles) <- basename(dirname(resFiles))
      
      txi <- tximport(resFiles, type="kallisto", tx2gene=tx2gene)
      head(txi$counts[,1:5])
      dim(txi$counts)
      
      xx <- txi$abundance
      gene.exp <- Biobase::ExpressionSet(log2(xx + 0.001))
      fData(gene.exp) <- toil.genes[featureNames(gene.exp),]
      pData(gene.exp) <- samples_annotation[sampleNames(gene.exp),]
      annotation(gene.exp) <- "rnaseq"
      
      xx <- txi$counts
      gene.count <- Biobase::ExpressionSet(log2(xx + 1))
      fData(gene.count) <- toil.genes[featureNames(gene.count),]
      pData(gene.count) <- samples_annotation[sampleNames(gene.count),]
      annotation(gene.count) <- "rnaseq"
      
      txii <- tximport(resFiles, type="kallisto", txOut=T)
      
      xx <- txii$abundance
      transcript.exp <- Biobase::ExpressionSet(log2(xx[,1:length(resFiles)] + 0.001))
      fData(transcript.exp) <- toil.transcripts[featureNames(transcript.exp),]
      pData(transcript.exp) <- samples_annotation[sampleNames(transcript.exp),]
      annotation(transcript.exp) <- "isoforms"
      
      xx <- txii$counts
      transcript.count <- Biobase::ExpressionSet(log2(xx[,1:length(resFiles)] + 1))
      fData(transcript.count) <- toil.transcripts[featureNames(transcript.count),]
      pData(transcript.count) <- samples_annotation[sampleNames(transcript.count),]
      annotation(transcript.count) <- "isoforms"
      
      return(list("rnaseq"=gene.exp, 
                  "rnaseq.counts"=gene.count, 
                  "isoforms"=transcript.exp, 
                  "isoforms.counts"=transcript.count))
    }
    

    rnaseq.sampleinfo <- read.csv("/pfs/downloadrna/Kallisto_0.43.1_processed/Kallisto_0.43.1_processed/JRGraySRRMapping.csv", stringsAsFactors=FALSE, row.names=1)
    
    rnaseq.sampleinfo[ , "cellid"] <- as.character(matchToIDTable(ids=rnaseq.sampleinfo[ , "cellid"], tbl=curationCell, column = "GRAY.cellid", returnColumn = "unique.cellid"))
   
    rnaseq <- summarizeRnaSeq(dir="/pfs/downloadrna/Kallisto_0.43.1_processed/Kallisto_0.43.1_processed", 
                                tool="kallisto", 
                                features_annotation="/pfs/downloadrna/Kallisto_0.43.1_processed/Kallisto_0.43.1_processed/Gencode.v23.annotation.RData",
                                samples_annotation=rnaseq.sampleinfo)
    

    
    GRAY2013 <- PharmacoSet(molecularProfiles=rnaseq,
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
    
    saveRDS(GRAY2013,file="/pfs/out/GRAY_2013.rds")
    
    return (GRAY2013)
    
  }

getGRAYP(verbose=FALSE, nthread=1)
