library(tidyverse)
library(glue)
library(data.table)

library(bigsnpr)
library(BEDMatrix) # read bed from file & avoid loading into RAM

library(devtools)
load_all("~/git/variani/biglmmz") # impute/scale large matrices 

## args
args <- commandArgs(trailingOnly = TRUE)
no_args <- is.na(args[1])

trait <- ifelse(no_args, "bmi", args[1])
file_out <- ifelse(no_args, "tmp.rds", args[2])

basename_out <- gsub("\\.rds$", "", file_out)

## data
file_bed <- switch(trait,
  "bmi" = "dat/bed_top/body_BMIz_orth.top.5000.bed",
  "height" = "dat/bed_top/body_HEIGHTz_orth.top.5000.bed",
  "hip" = "dat/bed_top/body_HIP_orth.top.5000.bed",
  "waist" = "dat/bed_top/body_WAIST_orth.top.5000.bed",
  "weight" = "dat/bed_top/body_WEIGHT_orth.top.5000.bed",
  "whr" = "dat/bed_top/body_WHR_orth.top.5000.bed",
  stop("switch trait"))
file_bim <- gsub(".bed$", ".bim", file_bed)
file_fam <- gsub(".bed$", ".fam", file_bed)

## bed
snps <- fread(file_bim, select = 2) %>% .[[1]] %>% as.character
ids <- fread(file_fam, select = 2) %>% .[[1]] %>% as.character
bed <- BEDMatrix(file_bed, n = length(ids), p = length(snps))

## sync
ids_sync <- read_lines("out/sync.ids")
stopifnot(all(ids_sync %in% ids))

tab_sync <- tibble(id = ids_sync)
tab_bed <- tibble(id = ids, row = seq(length(ids)))
tab_merged <- left_join(tab_sync, tab_bed, by = "id")

rows_fbm <- tab_merged$row

ids_fbm <- ids[rows_fbm]
stopifnot(ids_fbm == ids_sync)

## Z matrix (FBM)
cat(" - convert BED to FBM\n")

file_b <- tempfile()
file_bk <- glue("{file_b}.bk")
file_rds <- snp_readBed2(file_bed, backingfile = file_b, ind.row = rows_fbm)

bigsnp <- snp_attach(file_rds)
snpnames <- bigsnp$map[[2]]
Z0 <- bigsnp$genotypes

stopifnot(nrow(Z0) == length(ids_sync))
stopifnot(all(bigsnp$fam$sample.ID == ids_sync))

## copy to the output file
cat(" - copy to output file\n")
Z <- big_copy(Z0, backingfile = basename_out)
unlink(c(file_bk, file_rds))

cat(" - impute missing with zero\n")
impute_Z(Z, impute = "zero")
Z$save()

cat(" - write meta data\n")
file_ids <- gsub(".rds$", ".ids", file_out)
write_lines(bigsnp$fam$sample.ID, file_ids)
file_variants <- gsub(".rds$", ".variants", file_out)
write_lines(bigsnp$map$marker.ID, file_variants)

cat(" - done!\n")


