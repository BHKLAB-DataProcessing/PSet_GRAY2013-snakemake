library(downloader)
library(PharmacoGx)
library(openxlsx)

options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
data_dir <- paste0(args[1], "data")
download_dir <- paste0(args[1], "download")
processed_dir <- paste0(args[1], "processed")

# data_dir <- '/Users/minoru/Code/bhklab/DataProcessing/PSet/getGRAY2013/data'
# download_dir <- '/Users/minoru/Code/bhklab/DataProcessing/PSet/getGRAY2013/download'
# processed_dir <- '/Users/minoru/Code/bhklab/DataProcessing/PSet/getGRAY2013/processed'

badchars <- "[\xb5]|[]|[ ,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"

matchToIDTable <- function(ids, tbl, column, returnColumn = "unique.cellid") {
  sapply(ids, function(x) {
    myx <- grep(paste0("((///)|^)", Hmisc::escapeRegex(x), "((///)|$)"), tbl[, column])
    if (length(myx) > 1) {
      stop("Something went wrong in curating ids, we have multiple matches")
    }
    if (length(myx) == 0) {
      return(NA_character_)
    }
    return(tbl[myx, returnColumn])
  })
}

cell_all <- read.csv(file.path(download_dir, "cell_annotation_all.csv"), na.strings = c("", " ", "NA"))

## read cell and tissue curations

curationCell <- cell_all[which(!is.na(cell_all[, "GRAY.cellid"])), ]
curationTissue <- cell_all[which(!is.na(cell_all[, "GRAY.cellid"])), ]
curationCell <- curationCell[, c("unique.cellid", "GRAY.cellid", "Cellosaurus.Accession.id")]


curationTissue <- curationTissue[, c("unique.tissueid", "GRAY.tissueid")]

rownames(curationTissue) <- curationCell[, "unique.cellid"]
rownames(curationCell) <- curationCell[, "unique.cellid"]


# read drug curations
drug_all <- read.csv(file.path(download_dir, "drugs_with_ids.csv"), na.strings = c("", " ", "NA"))

curationDrug <- drug_all[which(!is.na(drug_all[, "GRAY.drugid"])), ]
curationDrug <- curationDrug[, c("unique.drugid", "GRAY.drugid")]

rownames(curationDrug) <- curationDrug[, "unique.drugid"]


## create sensitivity slices

getGRAYrawData <-
  function(result.type = c("array", "list")) {
    gray.raw.drug.sensitivity <- read.csv(file = file.path(data_dir, "gb-2013-14-10-r110-s9.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    gray.raw.drug.sensitivity.list <- do.call(c, apply(gray.raw.drug.sensitivity, 1, list))

    concentrations.no <- length(grep("^c[1-9]", colnames(gray.raw.drug.sensitivity)))

    if (result.type == "array") {
      ## create the gray.drug.response object including information viablilities and concentrations for each cell/drug pair
      obj <- array(NA, dim = c(length(unique(gray.raw.drug.sensitivity[, "cellline"])), length(unique(gray.raw.drug.sensitivity[, "compound"])), 2, concentrations.no), dimnames = list(unique(gray.raw.drug.sensitivity[, "cellline"]), unique(gray.raw.drug.sensitivity[, "compound"]), c("concentration", "viability"), 1:concentrations.no))
    }
    fnexperiment <-
      function(values) {
        cellline <- values["cellline"]
        drug <- values["compound"]
        doses <- as.numeric(values[grep("^c[1-9]", names(values))]) * 10^6 # micro molar

        if (concentrations.no > length(doses)) {
          doses <- c(doses, rep(NA, concentrations.no - length(doses)))
        }

        # responses <- as.numeric(unlist(strsplit(input.matrix["Activity Data\n(raw median data)"], split=",")))  #nature paper raw data

        responses <- NULL
        background <- median(as.numeric(values[grep("^background_od", names(values))]))
        ctrl <- median(as.numeric(values[sprintf("od%s.%s", 0, 1:3)])) - background
        ctrl <- ifelse(ctrl < 0, 1, ctrl)
        for (i in 1:concentrations.no)
        {
          res <- median(as.numeric(values[sprintf("od%s.%s", i, 1:3)])) - background
          res <- ifelse(res < 0, 0, res)
          responses <- c(responses, res / ctrl)
        }
        responses <- responses * 100
        if (result.type == "array") {
          obj[cellline, drug, "concentration", 1:length(doses)] <<- doses
          obj[cellline, drug, "viability", 1:length(responses)] <<- responses
        } else {
          return(list(cell = cellline, drug = drug, doses = doses, responses = responses)) # paste(doses, collapse = ","), responses=paste(responses, collapse = ",")))
        }
      }

    gray.raw.drug.sensitivity.res <- mapply(fnexperiment, values = gray.raw.drug.sensitivity.list)
    if (result.type == "array") {
      return(list("data" = obj, "concentrations.no" = concentrations.no))
    } else {
      return(list("data" = gray.raw.drug.sensitivity.res, "concentrations.no" = concentrations.no))
    }
  }

raw.sensitivity <- getGRAYrawData(result.type = "list")
con_tested <- raw.sensitivity$concentrations.no
raw.sensitivity <- t(raw.sensitivity$data)
raw.sensitivity <- t(apply(raw.sensitivity, 1, function(x) {
  unlist(x)
}))
rownames(raw.sensitivity) <- sprintf("drugid_%s_%s", as.character(raw.sensitivity[, 2]), as.character(raw.sensitivity[, 1]))

## handle replicates
tt <- rownames(raw.sensitivity)
for (i in 1:length(tt)) {
  xx <- which(tt == tt[i])
  if (length(xx) > 1) {
    for (j in 1:length(xx)) {
      tt[xx[j]] <- paste(tt[xx[j]], j, sep = "_")
    }
  }
}
rownames(raw.sensitivity) <- tt

tt <- paste0("doses", con_tested)
sensitivity.info <- raw.sensitivity[, c(1, 2, 3, grep(tt, colnames(raw.sensitivity)))]
colnames(sensitivity.info) <- c("cellid", "drugid", "min.Dose.uM", "max.Dose.uM")
sensitivity.info <- cbind(sensitivity.info, "nbr.conc.tested" = con_tested)
as.character(matchToIDTable(ids = sensitivity.info[, "cellid"], tbl = curationCell, column = "GRAY.cellid", returnColumn = "unique.cellid"))

sensitivity.info[, "cellid"] <- as.character(matchToIDTable(ids = sensitivity.info[, "cellid"], tbl = curationCell, column = "GRAY.cellid", returnColumn = "unique.cellid"))
sensitivity.info[, "drugid"] <- as.character(matchToIDTable(ids = sensitivity.info[, "drugid"], tbl = curationDrug, column = "GRAY.drugid", returnColumn = "unique.drugid"))

raw.sensitivity <- raw.sensitivity[, -c(1, 2)]
raw.sensitivity <- array(c(as.matrix(as.numeric(raw.sensitivity[, 1:con_tested])), as.matrix(as.numeric(raw.sensitivity[, (con_tested + 1):(2 * con_tested)]))), c(nrow(raw.sensitivity), con_tested, 2),
  dimnames = list(rownames(raw.sensitivity), colnames(raw.sensitivity[, 1:con_tested]), c("Dose", "Viability"))
)


save(raw.sensitivity, sensitivity.info, tt, con_tested, file = file.path(processed_dir, "drug_norm_post.RData"))


raw.sensitivity.x <- parallel::splitIndices(nrow(raw.sensitivity), floor(nrow(raw.sensitivity) / 1000))

dir.create(file.path(processed_dir, "slices"))

for (i in seq_along(raw.sensitivity.x)) {
  slce <- raw.sensitivity[raw.sensitivity.x[[i]], , ]
  saveRDS(slce, file = paste0(processed_dir, "/slices/gray2013_raw_sens_", i, ".rds"))
}
zip(zipfile=file.path(processed_dir, 'raw_sense_slices.zip'), files=list.files(file.path(processed_dir, 'slices'), full.names = TRUE))
unlink(file.path(processed_dir, 'slices'), recursive=TRUE)