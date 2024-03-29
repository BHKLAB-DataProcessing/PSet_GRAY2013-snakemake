library(PharmacoGx)
library(readxl)
library(openxlsx)
library(tximport)
library(Biobase)
library(data.table)
library(reshape2)
library(CoreGx)
library(SummarizedExperiment)
library(parallel)
# library(biocompute)

options(stringsAsFactors = FALSE)
args = commandArgs(trailingOnly=TRUE)

data_dir <- paste0(args[[1]], "data")
download_dir <- paste0(args[[1]], "download")
processed_dir <- paste0(args[[1]], "processed")
output_dir <- args[[1]]
filename <- args[[2]]
tools <- args[[3]]
transcriptome <- args[[4]]

# data_dir <- "/Users/minoru/Code/bhklab/DataProcessing/PSet/getGRAY2013/data"
# download_dir <- "/Users/minoru/Code/bhklab/DataProcessing/PSet/getGRAY2013/download"
# processed_dir <- "/Users/minoru/Code/bhklab/DataProcessing/PSet/getGRAY2013/processed"
# output_dir <- "/Users/minoru/Code/bhklab/DataProcessing/PSet/getGRAY2013/"
# tools <- "Kallisto-0.46.1"
# transcriptome <- "Gencode_v33"

tools <- gsub("-", "_", tools)
untar(file.path(download_dir, paste0(tools, ".tar.gz")), exdir = download_dir)
tool_path <- expand.grid(a = tools, b = transcriptome)
tool_path <- paste0(tool_path$a, "_", tool_path$b)

print(tool_path)

standardize <- args[grep("filtered", args)]

standardizeRawDataConcRange <- function(sens.info, sens.raw){
  unq.drugs <- unique(sens.info$drugid)
  
  conc.m <- data.table(melt(sens.raw[,,1], as.is=TRUE))
  conc.m[,drugid := sens.info$drugid[match(Var1, rownames(sens.info))]]
  conc.ranges <- conc.m[,.(l = min(value, na.rm=T), r = max(value, na.rm=T)), c("drugid", "Var1")]
  conc.ranges[,Var1 := NULL]
  conc.ranges <- conc.ranges[,unique(.SD), drugid]	
  # conc.ranges[,N := .N, drugid]
  conc.ranges.disj <- conc.ranges[, {sq <- sort(unique(c(l,r))); 
  l = sq[seq(1,length(sq)-1)];
  r = sq[seq(2,length(sq))];
  .(l=l,r=r)}, drugid]
  ## Function below returns all consecutive ranges of ints between 1 and N
  returnConsInts <- function(N) {
    stopifnot(N>0)
    unlist(sapply(seq(1,N), function(ii) return(sapply(seq(ii, N), function(jj) return(seq(ii,jj))))), recursive=FALSE)
  }
  rangeNoHoles <- function(indicies, lr.tbl){
    if(length(indicies) == 1) return(TRUE)
    sq <- seq(indicies[1], indicies[length(indicies)]-1)
    all(lr.tbl[["l"]][sq+1] <= lr.tbl[["r"]][sq])
  }
  per.drug.range.indicies <- sapply(conc.ranges.disj[,.N,drugid][,N], returnConsInts)
  
  names(per.drug.range.indicies) <- conc.ranges.disj[,unique(drugid)] ## checked this: conc.ranges.disj[,.N,drugid][,drugid] == conc.ranges.disj[,unique(drugid)]
  
  
  # Check if there are any holes in the chosen range combination
  per.drug.range.indicies <- sapply(names(per.drug.range.indicies), function(drug){
    
    lr.tbl <- conc.ranges.disj[drugid == drug]
    per.drug.range.indicies[[drug]][sapply(per.drug.range.indicies[[drug]], rangeNoHoles, lr.tbl = lr.tbl)]
    
  })
  per.drug.range.indicies.2 <- sapply(names(per.drug.range.indicies), function(drug){
    
    lr.tbl <- conc.ranges.disj[drugid == drug]
    res <- t(sapply(per.drug.range.indicies[[drug]], function(x) return(c(lr.tbl[x[1],l], lr.tbl[x[length(x)],r]))))
    colnames(res) <- c("l", "r")
    res <- data.frame(res)
    res <- cbind(drugid = drug, res)
  }, simplify=FALSE)
  per.drug.range.indicies.dt <- rbindlist(per.drug.range.indicies.2)
  
  conc.ranges <- conc.m[,.(l = min(value, na.rm=T), r = max(value, na.rm=T)), c("drugid", "Var1")]
  setkey(conc.m, Var1)
  conc.m <- na.omit(conc.m)
  setkey(conc.m, drugid, Var1, value)
  setkey(conc.ranges, drugid, l, r)
  # tic()
  ## NOTE:: Data.table used for maximum speed. Probably possible to do this more intelligently by 
  ## NOTE:: being aware of which conditions overlap, but its fast enough right now as it is.
  chosen.drug.ranges <- lapply(unq.drugs, function(drug){
    num.points.in.range <- apply(per.drug.range.indicies.dt[drugid==drug, .(l,r)], 1, function(rng){
      conc.m[drugid==drug][conc.ranges[drugid==drug][l<=rng["l"]][r>=rng["r"],Var1], on="Var1"][value >= rng["l"]][value <= rng["r"],.N]
      # conc.m[drugid==drug][, Var1]
    })
    max.ranges <- per.drug.range.indicies.dt[drugid==drug][which(num.points.in.range==max(num.points.in.range))]
    max.ranges[which.max(log10(r) - log10(l)), ]
  })
  # toc()
  names(chosen.drug.ranges) <- sapply(chosen.drug.ranges, `[[`, "drugid")
  removed.experiments <- unlist(lapply(unq.drugs, function(drug){
    rng <- unlist(chosen.drug.ranges[[drug]][,.(l,r)])
    exp.out.range <- conc.ranges[drugid==drug][l>rng["l"] | r<rng["r"],Var1]
    return(exp.out.range)
  }))
  
  sens.raw[removed.experiments,,] <- NA_real_
  conc.ranges.kept <- conc.ranges[!Var1 %in% removed.experiments]
  
  for(drug in unq.drugs){
    rng <- unlist(chosen.drug.ranges[[drug]][,.(l,r)])
    myx <- conc.ranges.kept[drugid==drug,Var1]
    doses <- sens.raw[myx, ,"Dose"]
    which.remove <- (doses < rng["l"] | doses > rng["r"])
    sens.raw[myx, ,"Dose"][which(which.remove,arr.ind=TRUE)] <- NA_real_
    sens.raw[myx, ,"Viability"][which(which.remove,arr.ind=TRUE)] <- NA_real_
    
    ## Annotate sens info with chosen range
    sens.info[sens.info$drugid==drug,"chosen.min.range"] <- rng["l"]
    sens.info[sens.info$drugid==drug,"chosen.max.range"] <- rng["r"]
  }
  sens.info$rm.by.conc.range <- FALSE
  sens.info[removed.experiments,"rm.by.conc.range"] <- TRUE
  
  return(list("sens.info" = sens.info, sens.raw = sens.raw))
}


#filter noisy curves from PSet (modified function to take into account standardized conc range)
filterNoisyCurves2 <- function(pSet, epsilon=25 , positive.cutoff.percent=.80, mean.viablity=200, nthread=1) {
  acceptable <- mclapply(rownames(sensitivityInfo(pSet)), function(xp) {
    #for(xp in rownames(sensitivityInfo(pSet))){
    drug.responses <- as.data.frame(apply(pSet@sensitivity$raw[xp , ,], 2, as.numeric), stringsAsFactors=FALSE)
    if (!all(is.na(drug.responses))){
      
      
      drug.responses <- drug.responses[complete.cases(drug.responses), ]
      doses.no <- nrow(drug.responses)
      drug.responses[,"delta"] <- .computeDelta(drug.responses$Viability)
      
      delta.sum <- sum(drug.responses$delta, na.rm = TRUE)
      
      max.cum.sum <- .computeCumSumDelta(drug.responses$Viability)
      
      if ((table(drug.responses$delta < epsilon)["TRUE"] >= (doses.no * positive.cutoff.percent)) &
          (delta.sum < epsilon) &
          (max.cum.sum < (2 * epsilon)) &
          (mean(drug.responses$Viability) < mean.viablity)) {
        return (xp)
      }
    }
    
  }, mc.cores=nthread)
  acceptable <- unlist(acceptable)
  noisy <- setdiff(rownames(sensitivityInfo(pSet)), acceptable)
  return(list("noisy"=noisy, "ok"=acceptable))
}

.computeDelta <- function(xx ,trunc = TRUE) {
  xx <- as.numeric(xx)
  if(trunc)
  {
    return(c(pmin(100, xx[2:length(xx)]) - pmin(100, xx[1:length(xx)-1]), 0))
  }else{
    return(c(xx[2:length(xx)] - xx[1:length(xx)-1]), 0)
  }
}

#' @importFrom utils combn
.computeCumSumDelta <- function(xx, trunc = TRUE) {
  xx <- as.numeric(xx)
  if(trunc) {
    xx <- pmin(xx, 100)
  }
  tt <- t(combn(1:length(xx), 2 , simplify = TRUE))
  tt <- tt[which(((tt[,2] - tt[,1]) >= 2) == TRUE),]
  if (is.null(nrow(tt))){
    tt <- matrix(tt, ncol = 2)
  }
  cum.sum <- unlist(lapply(1:nrow(tt), function(x){xx[tt[x,2]]-xx[tt[x,1]]}))
  return(max(cum.sum))
}    
    
	  
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

cell_all <- read.csv(file = file.path(download_dir, "cell_annotation_all.csv"), na.strings=c("", " ", "NA"))
curationCell <- cell_all[which(!is.na(cell_all[ , "GRAY.cellid"])),]
curationCell <- curationCell[ , c("unique.cellid", "GRAY.cellid")]
rownames(curationCell) <- curationCell[ , "unique.cellid"]

curationTissue <- cell_all[which(!is.na(cell_all[ , "GRAY.cellid"])),]
curationTissue <- curationTissue[ , c("unique.tissueid", "GRAY.tissueid")]
rownames(curationTissue) <- curationCell[ , "unique.cellid"]

drug_all <- read.csv(file.path(download_dir, "drugs_with_ids.csv"), na.strings=c("", " ", "NA"))
curationDrug <- drug_all[which(!is.na(drug_all[ , "GRAY.drugid"])),]
curationDrug <- curationDrug[,c("unique.drugid","GRAY.drugid")]
rownames(curationDrug) <- curationDrug[ , "unique.drugid"]


load(file.path(processed_dir, "drug_norm_post.RData"))

# cell information (cell slot)

cellineinfo <- read.xlsx(file.path(data_dir, "gb-2013-14-10-r110-s1.xlsx"), sheet = 1)
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

profiles <- read.xlsx(file.path(data_dir, "gb-2013-14-10-r110-s1.xlsx"), sheet = 1)
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


load(file.path(processed_dir, "profiles.RData"))

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
    
    
summarizeRnaSeq <- function (dir, features_annotation, samples_annotation, method) {
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
  
  if(features_annotation == file.path(download_dir, "Ensembl.v99.annotation.RData")){
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
  
  if(features_annotation == file.path(download_dir, "Ensembl.v99.annotation.RData")){
  #remove non-coding transcripts in ensembl 	  
  rownames(txii$abundance) <-  gsub("\\..*","",rownames(txii$abundance))
  txii$abundance[which(!rownames(txii$abundance)  %in% features_transcript$transcript_id)]
  missing_transcript <- rownames(txii$abundance)[which(!rownames(txii$abundance)  %in% features_transcript$transcript_id)]
  txii$abundance <- txii$abundance [-which(rownames(txii$abundance) %in% missing_transcript),]
  }
  	  
  xx <- txii$abundance
  transcript.exp <- Biobase::ExpressionSet(log2(xx[,1:length(resFiles)] + 0.001))
  if(features_annotation == file.path(download_dir, "Gencode.v33.annotation.RData") || features_annotation == file.path(download_dir, "Gencode.v33lift37.annotation.RData")){
  featureNames(transcript.exp) <- gsub("\\|.*","",featureNames(transcript.exp))
  fData(transcript.exp) <- features_transcript[featureNames(transcript.exp),]
  }else{
  fData(transcript.exp) <- features_transcript[featureNames(transcript.exp),]
  }
  pData(transcript.exp) <- samples_annotation[sampleNames(transcript.exp),]
  annotation(transcript.exp) <- "isoform"
  
	  
  if(features_annotation == file.path(download_dir, "Ensembl.v99.annotation.RData")){
  #remove non-coding transcripts in ensembl
  rownames(txii$counts) <-  gsub("\\..*","",rownames(txii$counts))
  txii$counts <- txii$counts [-which(rownames(txii$counts) %in% missing_transcript),]	  
  }	  
  xx <- txii$counts
  transcript.count <- Biobase::ExpressionSet(log2(xx[,1:length(resFiles)] + 1))
  if(features_annotation == file.path(download_dir, "Gencode.v33.annotation.RData") || features_annotation == file.path(download_dir, "Gencode.v33lift37.annotation.RData")){
  featureNames(transcript.count) <- gsub("\\|.*","",featureNames(transcript.count))
  fData(transcript.count) <- features_transcript[featureNames(transcript.count),]
  }else{
  fData(transcript.count) <- features_transcript[featureNames(transcript.count),]
  }
  pData(transcript.count) <- samples_annotation[sampleNames(transcript.count),]
  annotation(transcript.count) <- "isoform"
  
  return(
    list(
      "rnaseq"=gene.exp, 
      "rnaseq.counts"=gene.count, 
      "isoforms"=transcript.exp, 
      "isoforms.counts"=transcript.count
    )
  )
}

    

rnaseq.sampleinfo <- read.csv(file=file.path(download_dir, "JRGraySRRMapping.csv"), stringsAsFactors=FALSE, row.names=1)
    
rnaseq.sampleinfo[ , "cellid"] <- as.character(matchToIDTable(ids=rnaseq.sampleinfo[ , "cellid"], tbl=curationCell, column = "GRAY.cellid", returnColumn = "unique.cellid"))
rnaseq.sampleinfo[ , "batchid"] <- NA

rnaseq_results <- c()			   
for (r in 1:length(tool_path)){
  print(tool_path[r])
  if (length(grep(pattern = 'Kallisto', x = tool_path[r])) > 0){
    tool <- sub("(_[^_]+)_.*", "\\1", tool_path[r])
    # tdir = paste0("gray_rnaseq_",gsub(".","_",tolower(tool), fixed = T), "/",  tool, "/", tool, "/")  
    rnatool="kallisto"	  
  } else {
    tool <- sub("(_[^_]+)_.*", "\\1", tool_path[r])
    # tdir = paste0("gray_rnaseq_",gsub(".","_",tolower(tool), fixed = T), "/",  tool, "/", tool, "/")
    rnatool="salmon"	  
  }
  
  
  if (length(grep(pattern = 'lift37', x = tool_path[r])) > 0){
    annot = file.path(download_dir, "Gencode.v33lift37.annotation.RData")
  } else if (length(grep(pattern = 'v33', x = tool_path[r])) > 0){
    annot = file.path(download_dir, "Gencode.v33.annotation.RData")
  } else {
    annot = file.path(download_dir, "Ensembl.v99.annotation.RData")
  }
  print(annot)
  
  #print(tdir)
  rnaseq <- summarizeRnaSeq(
    dir=file.path(download_dir, tool, tool_path[r]),
    features_annotation=annot,
    samples_annotation=rnaseq.sampleinfo,
		method = rnatool
	)
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
			   
drug_all <- read.csv(file.path(download_dir, "drugs_with_ids.csv"), na.strings=c("", " ", "NA"))
drug_all <- drug_all[which(!is.na(drug_all[ , "GRAY.drugid"])),]
drug_all <- drug_all[ , c("unique.drugid", "GRAY.drugid","smiles","inchikey","cid","FDA")]
rownames(drug_all) <- drug_all[ , "unique.drugid"]

drug_all <- drug_all[rownames(druginfo),]
druginfo[,c("smiles","inchikey","cid","FDA")] <- drug_all[,c("smiles","inchikey","cid","FDA")] 
 	 
z <- list()
z <- c(z,c(rnaseq_results))
				   
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
		 
if (length(standardize) > 0){

  # standardize <- standardizeRawDataConcRange(sens.info = sensitivity.info, sens.raw = raw.sensitivity)
  # sensitivity.info <- standardize$sens.info
  # raw.sensitivity <- standardize$sens.raw

} else {
  print("unfiltered PSet")
}	 
		 				   

GRAY2013 <- PharmacoGx::PharmacoSet(
  molecularProfiles=z,
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
  datasetType="sensitivity"
)

if (length(standardize) > 0){
 noisy_out <- filterNoisyCurves2(GRAY2013)
 print("filter done")
 GRAY2013@sensitivity$profiles[noisy_out$noisy, ] <- NA

} else {
  print("unfiltered PSet")
}

GRAY2013@annotation$version <- 2
saveRDS(GRAY2013,file=paste0(output_dir, filename))				   
# dataset <- "GRAY"
# #output ORCESTRA_ID and Pachyderm commit id
# write.table(dataset, file="/pfs/out/dataset.txt", row.names = F ,quote = F, sep = "\t", col.names = F)	
# write.table(ORCESTRA_ID, file="/pfs/out/orcestra_id.txt", row.names = F ,quote = F, sep = "\t", col.names = F)				   
# pach_commit_id <- Sys.getenv("PACH_OUTPUT_COMMIT_ID")
# write.table(pach_commit_id, file="/pfs/out/commit_id.txt", row.names = F ,quote = F, sep = "\t", col.names = F) 
# 				   
# 				   
# 				   
# ###CREATE BIOCOMPUTE OBJECT###
# tools <- grep(pattern = 'Kallisto|Salmon', x = rnaseq_select)
# tools <- rnaseq_select[tools]
# tools <- gsub("-", "_", tools)
# 				   
# tool <- gsub("\\_.*","", tools)
# version <- gsub(".*_","", tools)
# 
# print(tool)
# print(version)				   
# 
# if (length(tools) > 1){
#   
#   if (isTRUE(any(duplicated(tool)))){
#     
#     tool_name <- tool[1]
#     version <- c(version[1],version[2])
#     
#     if (tool[1] == "Kallisto" || tool[2] == "Kallisto"){
#       
#       parameter <- "-t"
#       parameter_value <- "2"
#       uri <- rep("https://pachterlab.github.io/kallisto/",2)
#       
#     } else{
#       parameter <- c("-l", "--validateMappings")
#       parameter_value <- c("A", "")
#       
#       uri <- rep("https://combine-lab.github.io/salmon/",2)
#     }
#     
#   } else {
#     
#     tool <- c(tool[1], tool[2])
#     v <- version
#     #version <- c(version[1],version[2])
#     if (v[1] == "0.8.2" || v[2] == "0.8.2"){
#       parameter <- c("-t","-l") 
#       parameter_value <- c("2", "A")
#     } else{
#       parameter <- c("-t","-l", "--validateMappings")
#       parameter_value <- c("2", "A","")
#     }
#     uri <- c("https://pachterlab.github.io/kallisto/","https://combine-lab.github.io/salmon/")
#   }
#   
# } else {
#   
#   if (tool == "Kallisto"){
#     
#     parameter <- "-t"
#     parameter_value <- c("2")
#     uri <- "https://pachterlab.github.io/kallisto/"
#   } else{
#     if (version == "0.8.2"){
#       parameter <- paste0("-l")
#       parameter_value <- c("A")
#     } else{
#       parameter <- c("-l", "--validateMappings")
#       parameter_value <- c("A","")
#     }
#     uri <- "https://combine-lab.github.io/salmon/"
#   }
# }
# 
# ###Selection Transcriptome###
# 
# if (transcriptome == "Gencode_v33"){
#   
#   transcriptome_link <- "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.transcripts.fa.gz"
#   
# } else if (transcriptome == "Gencode_v33lift37"){
#   
#   transcriptome_link <- "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_23/GRCh37_mapping/gencode.v23lift37.annotation.gtf.gz"
# } else {
#   
#   transcriptome_link <- "ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz"
# }
# 
# print(tool)
# print(version)
# print(parameter)
# print(parameter_value)
# print(uri)
# print(transcriptome)
# print(transcriptome_link)
# 
# ###########################
# #####Provenance Domain#####
# ###########################
# 
# #Created and modified dates
# #Sys.setenv(TZ = "EST")
# created <- as.POSIXct(Sys.time(), format = "%Y-%m-%dT%H:%M:%S", tz = "EST")
# modified <- as.POSIXct(Sys.time(), format = "%Y-%m-%dT%H:%M:%S", tz = "EST")
# 
# #Contributions
# contributors <- data.frame(
#   "name" = c("Anthony Mammoliti", "Petr Smirnov", "Benjamin Haibe-Kains"),
#   "affiliation" = c(rep("University Health Network", 3)),
#   "email" = c("anthony.mammoliti@uhnresearch.ca", "petr.smirnov@utoronto.ca", "Benjamin.Haibe-Kains@uhnresearch.ca"),
#   "contribution" = c("createdBy","createdBy","authoredBy"),
#   "orcid" = c(NA,NA,"https://orcid.org/0000-0002-7684-0079"),
#   stringsAsFactors = FALSE
# )
# 
# #License
# license <- "https://opensource.org/licenses/Apache-2.0"
# 
# #Name of biocompute object
# name <- "GRAY_2013"
# 
# #Version of biocompute object
# bio_version <- "1.0.0"
# 
# #Embargo (none)
# embargo <- c()
# 
# #Derived from and obsolete after (none)
# derived_from <- c()
# obsolete_after <- c()
# 
# #reviewers (none)
# review <- c()
# 
# #compile domain
# provenance <- compose_provenance_v1.3.0(
#   name, bio_version, review, derived_from, obsolete_after,
#   embargo, created, modified, contributors, license
# )
# provenance %>% convert_json()
# 
# 
# ############################
# #####Description Domain#####
# ############################
# times_rnaseq <- as.POSIXct("2020-01-20T1:04:10", format = "%Y-%m-%dT%H:%M:%S", tz = "EST")
# #Keywords and platform info
# keywords <- c("Biomedical", "Pharmacogenomics", "Cellline", "Drug")
# platform <- c("Pachyderm", "ORCESTRA (orcestra.ca)", "Linux/Ubuntu")
# 
# #Metadata for each pipeline step
# pipeline_meta <- data.frame(
#   "step_number" = c("1","2","3","4"),
#   "name" = c("Expression processing",
#              "Curated Sample and treatment identifier compilation",
#              "Drug sensitivity processing",
#              "Build data object"),
#   "description" = c("Pseudoalignment of reads", 
#                     "Download of appropriate sample and treatment identifiers from GitHub (curations performed by BHK Lab - http://bhklab.ca)",
#                     "Process sensitivity data",
#                     "Building of ORCESTRA data object"),
#   "version" = c(1.0,1.0,1.0,1.0),
#   stringsAsFactors = FALSE
# )
# 
# #Inputs for each pipeline step
# pipeline_input <- data.frame(
#   "step_number" = c("1","1","2","2","3","4"),
#   "filename" = c("Raw RNA-seq data",
#                  "Transcriptome",
#                  "Sample annotation data",
#                  "Treatment annotations",
#                  "Raw sensitivity data",
#                  "Script for data object generation"),
#   "uri" = c(
#     "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48213",
#     transcriptome_link,
#     "https://github.com/BHKLAB-Pachyderm/Annotations/blob/master/cell_annotation_all.csv",
#     "https://github.com/BHKLAB-Pachyderm/Annotations/blob/master/drugs_with_ids.csv",
#     "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3937590/bin/gb-2013-14-10-r110-S9.txt",
#     "https://github.com/BHKLAB-Pachyderm/getGRAY2013/GRAY.R"
#   ),
#   "access_time" = c(times_rnaseq,times_rnaseq,created,created,created,created),
#   stringsAsFactors = FALSE
# )
# 
# 
# #Outputs for each pipeline step
# pipeline_output <- data.frame(
#   "step_number" = c("1","2","2","3","3","3","4"),
#   "filename" = c("Processed expression data", 
#                  "Downloaded sample annotations",
#                  "Downloaded treatment annotations",
#                  "Downloaded raw sensitivity data",
#                  "Processed sensitivity data in parallel",
#                  "Compiled sensitivity data",
#                  "Data object"),
#   "uri" = c(
#     "https://orcestradata.blob.core.windows.net/gray/RNA-seq/",
#     "/pfs/downAnnotations/cell_annotation_all.csv",
#     "/pfs/downAnnotations/drugs_with_ids.csv",
#     "/pfs/GRAYRawSensitivity/drug_norm_post.RData",
#     "/pfs/calculategray2013Raw/slice_$.rds",
#     "/pfs/gray2013ProfilesAssemble/profiles.RData",
#     "/pfs/out/GRAY.rds"
#   ),
#   "access_time" = c(times_rnaseq,created,created,created,created,created,created),
#   stringsAsFactors = FALSE
# )
# 
# #xref (none)
# xref <- c()
# 
# #pipeline prereq (none)
# pipeline_prerequisite <- c()
# 
# #compile domain
# description <- compose_description_v1.3.0(
#   keywords, xref, platform,
#   pipeline_meta, pipeline_prerequisite, pipeline_input, pipeline_output
# )
# description %>% convert_json()
# 
# 
# ############################
# ######Execution Domain######
# ############################
# 
# script <- c()
# script_driver <- c()
# 
# #software/tools and its versions used for data object creation
# software_prerequisites <- data.frame(
#   "name" = c(tool, "Pachyderm", "Docker Image"),
#   "version" = c(version, "1.9.3", "v3"),
#   "uri" = c(
#     uri,
#     "https://www.pachyderm.com", "https://hub.docker.com/r/bhklab/pharmacogx2.0"
#   ),
#   stringsAsFactors = FALSE
# )
# 
# software_prerequisites[,"access_time"] <- rep(NA, length(software_prerequisites$name))
# software_prerequisites[,"sha1_chksum"] <- rep(NA, length(software_prerequisites$name))
# 
# external_data_endpoints <- c()
# environment_variables <- c()
# 
# execution <- compose_execution_v1.3.0(
#   script, script_driver, software_prerequisites, external_data_endpoints, environment_variables
# )
# execution %>% convert_json()
# 
# 
# ############################
# ######Extension Domain######
# ############################
# 
# #repo of scripts/data used
# scm_repository <- data.frame("extension_schema"= c("https://github.com/BHKLAB-Pachyderm"))
# scm_type <- "git"
# scm_commit <- c()
# scm_path <- c()
# scm_preview <- c()
# 
# scm <- compose_scm(scm_repository, scm_type, scm_commit, scm_path, scm_preview)
# scm %>% convert_json()
# 
# extension <- compose_extension_v1.3.0(scm)
# extension %>% convert_json()
# 
# ############################
# ######Parametric Domain#####
# ############################
# 
# df_parametric <- data.frame(
#   "param" = parameter,
#   "value" = parameter_value,
#   "step" = c(1),
#   stringsAsFactors = FALSE
# )
# 
# parametric <- compose_parametric_v1.3.0(df_parametric)
# parametric %>% convert_json()
# 
# 
# 
# ############################
# ######Usability Domain######
# ############################
# 
# #usability of our data objects
# text <- c(
#   "Pipeline for creating GRAY 2013 data object through ORCESTRA (orcestra.ca), a platform for the reproducible and transparent processing, sharing, and analysis of biomedical data."
# )
# 
# usability <- compose_usability_v1.3.0(text)
# usability %>% convert_json()
# 
# 
# ######################
# ######I/O Domain######
# ######################
# 
# input_subdomain <- data.frame(
#   "filename" = c("Raw RNA-seq data",
#                  "Transcriptome",
#                  "Sample annotation data",
#                  "Treatment annotations",
#                  "Raw sensitivity data",
#                  "Script for data object generation"),
#   "uri" = c(
#     "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48213",
#     transcriptome_link,
#     "https://github.com/BHKLAB-Pachyderm/Annotations/blob/master/cell_annotation_all.csv",
#     "https://github.com/BHKLAB-Pachyderm/Annotations/blob/master/drugs_with_ids.csv",
#     "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3937590/bin/gb-2013-14-10-r110-S9.txt",
#     "https://github.com/BHKLAB-Pachyderm/getGRAY2013/GRAY.R"
#   ),
#   "access_time" = c(times_rnaseq,times_rnaseq,created,created,created,created),
#   stringsAsFactors = FALSE
# )
# 
# output_subdomain <- data.frame(
#   "mediatype" = c("tar.gz", "csv", "csv", "RData", "RDS","RData","RDS"),
#   "uri" = c(
#     "https://orcestradata.blob.core.windows.net/gray/RNA-seq/",
#     "/pfs/downAnnotations/cell_annotation_all.csv",
#     "/pfs/downAnnotations/drugs_with_ids.csv",
#     "/pfs/GRAYRawSensitivity/drug_norm_post.RData",
#     "/pfs/calculategray2013Raw/slice_$.rds",
#     "/pfs/gray2013ProfilesAssemble/profiles.RData",
#     "/pfs/out/GRAY.rds"
#   ),
#   "access_time" = c(times_rnaseq,created,created,created,created,created,created),
#   stringsAsFactors = FALSE
# )
# 
# io <- compose_io_v1.3.0(input_subdomain, output_subdomain)
# io %>% convert_json()
# 
# 
# ########################
# ######Error Domain######
# ########################
# 
# empirical <- c()
# algorithmic <- c()
# 
# error <- compose_error(empirical, algorithmic)
# error %>% convert_json()
# 
# 
# ####Retrieve Top Level Fields####
# tlf <- compose_tlf_v1.3.0(
#   provenance, usability, extension, description,
#   execution, parametric, io, error
# )
# tlf %>% convert_json()
# 
# 
# ####Complete BCO####
# 
# bco <- biocompute::compose_v1.3.0(
#   tlf, provenance, usability, extension, description,
#   execution, parametric, io, error
# )
# bco %>% convert_json() %>% export_json("/pfs/out/GRAY_2013_BCO.json") %>% validate_checksum()
# 				   
