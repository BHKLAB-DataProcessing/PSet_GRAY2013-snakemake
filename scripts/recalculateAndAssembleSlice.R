library(PharmacoGx)

args <- commandArgs(trailingOnly = TRUE)
processed_dir <- paste0(args[1], "processed")

processed_dir <- "/Users/minoru/Code/bhklab/DataProcessing/PSet/getGRAY2013/processed"

unzip(file.path(processed_dir, "raw_sense_slices.zip"), exdir = file.path(processed_dir, "slices"), junkpaths = TRUE)

files <- list.files(file.path(processed_dir, "slices"), full.names = TRUE)
dir.create(file.path(processed_dir, "slices_recomp"))
for (file in files) {
  print(file)
  mybasenm <- basename(file)
  
  slice <- readRDS(file)
  
  res <- PharmacoGx:::.calculateFromRaw(slice)
  
  saveRDS(res, file = file.path(processed_dir, "slices_recomp", gsub(mybasenm, pattern = ".rds", replacement = "_recomp.rds", fixed = TRUE)))
}

recalc_files <- list.files(path = file.path(processed_dir, "slices_recomp"), full.names = TRUE)
recalc_files <- recalc_files[str_detect(recalc_files, "_recomp.rds")]
slices <- list()

for(fn in recalc_files){
  temp <- readRDS(fn)
  parTable <- do.call(rbind,temp[[3]])
  n <- cbind("AAC" = temp[[1]], "IC50" = temp[[2]], parTable) 
  slices[[fn]] <- n
}

res <- do.call(rbind, slices)

save(res, file=file.path(processed_dir, "profiles.RData"))