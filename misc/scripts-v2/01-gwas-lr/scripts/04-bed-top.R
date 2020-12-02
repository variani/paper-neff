library(tidyverse)
library(glue)
library(data.table)

library(bigsnpr)

library(devtools)
load_all("~/git/variani/matlm") # gwas assoc. model
load_all("~/git/variani/bigcov") # split bed by batches
load_all("~/git/variani/biglmmz") # impute/scale large matrices 

library(BEDMatrix) # read bed from file & avoid loading into RAM

## args
args <- commandArgs(trailingOnly = TRUE)
no_args <- is.na(args[1])

trait <- ifelse(no_args, "height", args[1])
n_top <- ifelse(no_args, 500, as.integer(args[2]))
file_out <- ifelse(no_args, "tmp.bk", args[3])

basename_out <- gsub("\\.bk$", "", file_out)

## testing
testing <- no_args
nrow1 <- 1e3

## data
file_bed <- "output/gen.bed"
file_bim <- gsub(".bed$", ".bim", file_bed)
file_clump <- glue("output/clump-plink/{trait}.tsv.gz")

## read top variants
clump <- read_tsv(file_clump) %>% arrange(pval) 
if(nrow(clump) < n_top) {
  warning(glue("#rows in clump ({nrow(clump)}) < #top ({n_top})"))
}

clump <- head(clump, n_top)
snps_top <- clump$snp

## bed
snps <- read_lines("output/gen.variants")
ids <- read_lines("output/gen.ids")
bed <- BEDMatrix(file_bed, n = length(ids), p = length(snps))
colnames(bed) <- snps

## Z matrix (FBM)
cat(" - convert BED to FBM\n")

file_b <- tempfile()
file_bk <- glue("{file_b}.bk")
cols_top <- which(snps %in% snps_top) %>% sort
if(testing) {
  file_rds <- snp_readBed2(file_bed, backingfile = file_b, ind.row = seq(nrow1), ind.col = cols_top)
} else {
  file_rds <- snp_readBed2(file_bed, backingfile = file_b, ind.col = cols_top)
}

bigsnp <- snp_attach(file_rds)
snpnames_Z <- bigsnp$map[[2]]
Z0 <- bigsnp$genotypes

## copy to the output file
cat(" - copy to output file\n")
Z <- big_copy(Z0, backingfile = basename_out)
unlink(c(file_bk, file_rds))

cat(" - impute missing with zero\n")
impute_Z(Z, impute = "zero")
Z$save()

cat(" - write meta data\n")
file_ids <- gsub(".bk$", ".ids", file_out)
write_lines(bigsnp$fam$sample.ID, file_ids)
file_variants <- gsub(".bk$", ".variants", file_out)
write_lines(bigsnp$map$marker.ID, file_variants)

cat(" - done!\n")

