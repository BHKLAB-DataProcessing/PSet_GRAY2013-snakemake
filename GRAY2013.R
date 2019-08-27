library(PharmacoGxPrivate)
library(readxl)
library(openxlsx)
library(tximport)
library(Biobase)

getGRAYP <-
  function (
    verbose=FALSE,
    nthread=1){
    
options(stringsAsFactors=FALSE)

badchars <- "[\xb5]|[]|[ ,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[\\]|[.]|[_]|[ ]"

#match to cell curations

matchToIDTableCELL <- function(ids,tbl, column) {
  sapply(ids, function(x) {
    myx <- grep(paste0("((///)|^)",x,"((///)|$)"), tbl[,column])
    if(length(myx) > 1){
      stop("Something went wrong in curating cell ids")
    }
    return(tbl[myx, "unique.cellid"])
  })
}


#match to drug curation 

matchToIDTableDRUG <- function(ids,tbl, column) {
  sapply(ids,function(x) {
    myx <- grep(paste0("((///)|^)",x,"((///)|$)"), tbl[,column])
    if(length(myx) > 1){
      stop("Something went wrong in curating drug ids")
    }
    return(tbl[myx, "unique.drugid"])
  })
}


cell_all <- read.csv(file = "/pfs/downAnnotations/cell_annotation_all.csv", na.strings=c("", " ", "NA"))
curationCell <- cell_all[which(!is.na(cell_all[ , "GRAY.cellid"])),]
curationCell <- curationCell[ , c("unique.cellid", "GRAY.cellid")]
rownames(curationCell) <- curationCell[ , "unique.cellid"]
removed_cells <- c("DU-4475","T47D_KBluc","HCC2157","MB157","HCC1500","SUM 190","HCC1007","184A1N4", "MDA-MB-435")
curationCell <- curationCell[which(!curationCell$unique.cellid %in% removed_cells),]
rownames(curationCell) <- curationCell$unique.cellid

curationTissue <- cell_all[which(!is.na(cell_all[ , "GRAY.cellid"])),]
curationTissue <- curationTissue[which(curationTissue$unique.cellid %in% curationCell$unique.cellid),]
curationTissue <- curationTissue[ , c("unique.tissueid", "GRAY.tissueid")]
rownames(curationTissue) <- curationCell[ , "unique.cellid"]

drug_all <- read.csv("/pfs/downAnnotations/drugs_with_ids.csv", na.strings=c("", " ", "NA"))
curationDrug <- drug_all[which(!is.na(drug_all[ , "GRAY.drugid"])),]
curationDrug <- curationDrug[,c("unique.drugid","GRAY.drugid")]
rownames(curationDrug) <- curationDrug[ , "unique.drugid"]

removed_drugs <- c("Nilotinib","CI-1040","PD 173074","GW843682X","Dichloroacetic acid","2-deoxyglucose","Celecoxib","LY-294002","Sulindac sulfide","Chloroquine","Ribavirin","TPCA-1","Cetuximab","Bromopyruvic acid","Methyl Glyoxol","Z-LL-Nva-CHO","TPT","Trastuzumab",        
                   "UO126","PS-1145","Tyrphostin AG 1024","GW5074","Ilomastat","QNZ",             
                   "TAPI-0","2C4","L779450","Chk2 Inhibitor II","Ro 32-0432","TCS JNK 5a",         
                   "PP242 hydrate")

curationDrug <- curationDrug[which(!curationDrug$unique.drugid %in% removed_drugs),]

load("/pfs/GRAYRawSensitivity/drug_norm_post.RData")

 ## cell information
    
    cellineinfo <- read.xlsx("/pfs/getGRAY2013/gb-2013-14-10-r110-s1.xlsx", sheet = 1)
    cellineinfo[!is.na(cellineinfo) & cellineinfo == ""] <- NA
    rn <- cellineinfo[-1, 1]
    cn <- t(cellineinfo[1, -1])
    cn <- gsub(badchars, ".", cn)
    cellineinfo <- cellineinfo[-1, -1]
    dimnames(cellineinfo) <- list(rn, cn)
    cellineinfo <- data.frame("cellid"=rn, "tissueid"="breast", cellineinfo[,1:10])
    cellineinfo <- cellineinfo[which(!is.na(cellineinfo$Transcriptional.subtype)), ]
    
    cellineinfo <- cellineinfo[which(!cellineinfo$cellid=="T47D_KBluc"),]
    cellineinfo <- cellineinfo[which(!cellineinfo$cellid=="MB157"),]
    c1 <- matchToIDTableCELL(cellineinfo$cellid, curationCell, "GRAY.cellid")
    c1[c1=="character(0)"] <- NA
    c1 <- c1[which(!is.na(c1))]
    cellineinfo$cellid <- c1
    cellineinfo$cellid[is.na(cellineinfo$cellid)]<-"NA"
    rownames(cellineinfo) <-  cellineinfo$cellid
    head(cellineinfo)
    
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
    profiles <- profiles[1:70, ]
    GI50 <- as.numeric(array(apply(profiles, 1, function(x)(x[indices]))))
    

    profilecol <- gsub("\\s*\\([^\\)]+\\)","",as.character(colnames(profiles)[indices]))
    drug <- matchToIDTableDRUG(profilecol, curationDrug, "GRAY.drugid")
    drug <- as.character(drug)
    drug <- rep(drug, times=70)
    #drug <- as.character(unique(unlist(y)))
    
    
    row_prof <- rownames(profiles)
    row_prof[grep("157",row_prof)] <- "MDAMB157"
    cell <- matchToIDTableCELL(row_prof, curationCell, "GRAY.cellid")
    cell <- as.character(cell)
    cell[cell=="character(0)"] <- "T47D_KBluc"
    cell <- cell[which(!is.na(cell))]
    cell <- rep(cell, each=90)
    
    GI50initial <- data.frame("cellid"=cell, "drugid"=drug, "GI50_published"=GI50)
    GI50match <- data.frame("cellid"=sensitivity.info[,"cellid"], "drugid"=sensitivity.info[,"drugid"], "GI50"=NA)
    
    GI50match$GI50 <- GI50initial[match(paste(GI50match$cellid,GI50match$drugid),paste(GI50initial$cellid,GI50initial$drugid)),"GI50_published"]
    
    sensitivity.profiles <- matrix(NA, dimnames = list(rownames(sensitivity.info), "GI50_published"), nrow=nrow(sensitivity.info))
    sensitivity.profiles[,"GI50_published"] <- GI50match$GI50
    
    load("/pfs/gray2013ProfilesAssemble/profiles.RData")
    
    sensitivity.profiles <- cbind(sensitivity.profiles, "AAC"=as.numeric(res[,"AAC"]), "IC50"=as.numeric(res[,"IC50"]), "HS"=as.numeric(res[,"HS"]), "E_inf"=as.numeric(res[,"E_inf"]), "EC50"=as.numeric(res[,"EC50"]))
    
    slope <- NULL
    for(exp in rownames(sensitivity.info)){
      slope <- c(slope, computeSlope(raw.sensitivity[exp, , "Dose"], raw.sensitivity[exp, , "Viability"])) #computeSlope (returns normalized slope of drug response curve)
    }
    
    names(slope) <- rownames(sensitivity.info)
    sensitivity.profiles <- cbind(sensitivity.profiles, "slope_recomputed"=slope)
    head(sensitivity.profiles)
    
    druginfo <- data.frame("drugid"=curationDrug$unique.drugid)
    rownames(druginfo) <- druginfo$drugid
    
    #RNA-seq Processed Data (Kallisto)
    
    load("/pfs/getGRAY2013/GRAY_Kallisto.RData")
    
    GRAY2013_updated <- PharmacoSet(molecularProfiles=GRAY_RNA,
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
    
    save(GRAY2013,file="/pfs/out/GRAY_2013.RData")
    
    return (GRAY2013)
    
    }
    
    getGRAYP(verbose=FALSE, nthread=1)
