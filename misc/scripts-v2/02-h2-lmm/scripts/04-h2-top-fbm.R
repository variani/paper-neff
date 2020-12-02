### inc
library(tidyverse)
library(glue)

library(BEDMatrix)
library(bigsnpr)

library(devtools)
load_all("~/git/variani/biglmmz/")

## arguments
args <- commandArgs(trailingOnly = TRUE)
no_args <- is.na(args[1])

trait <- ifelse(no_args, "height", args[1])
n_top <- ifelse(no_args, 1000, as.integer(args[2]))
file_out <- ifelse(no_args, "tmp.tsv.gz", args[3])

## parameters
pat_temp <- "04-h2-top-"

file_bed <- "output/gen.bed"
file_bim <- gsub(".bed$", ".bim", file_bed)
file_phen <- "output/phen.bed.tsv.gz"
file_clump <- glue("output/clump-plink/{trait}.tsv.gz")

## testing
testing <- no_args
nrow1 <- 1e3

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
file_b <- tempfile(paste0(pat_temp, "readBed2-")) %>% basename
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

file_z <- tempfile(paste0(pat_temp, "Z-")) %>% basename
file_z_bk <- glue("{file_z}.bk")
Z <- big_copy(Z0, type = "double", backingfile = file_z)

# clean bigsnp
unlink(c(file_bk, file_rds))

## scale Z
scale_Z(Z, impute = TRUE, M = ncol(Z))

## read phen.
phen <- read_tsv(file_phen, col_types = c("ccddddddddd"))
stopifnot(trait %in% colnames(phen))
y <- phen[[trait]] 
if(testing) {
  y <- head(y, nrow1)
}

## fit LMM
mod <- biglmmz(y, Z = Z, 
  scale = FALSE, impute = FALSE, 
  copy_Z = FALSE,
  REML = TRUE, verbose = 2)

## output table
tab <- mod$ess %>%
  mutate(trait = trait)

## clean file_z_bk
rm(mod)
unlink(file_z_bk)

## save
write_tsv(tab, file_out)
  
