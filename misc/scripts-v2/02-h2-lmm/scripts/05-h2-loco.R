### inc
library(tidyverse)
library(glue)
library(data.table)

library(BEDMatrix)
library(bigstatsr)

library(devtools)
load_all("~/git/variani/biglmmz/")

## arguments
args <- commandArgs(trailingOnly = TRUE)
no_args <- is.na(args[1])

trait <- ifelse(no_args, "height", args[1])
n_top <- ifelse(no_args, 500, as.integer(args[2]))
chr <- ifelse(no_args, 1, as.integer(args[3]))
file_out <- ifelse(no_args, "tmp.tsv.gz", args[4])

## parameters
pat_temp <- "05-h2-loco-"

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

## top variants for LOCO
bim <- fread(file_bim, select = 1:2) %>% as_tibble # columns V1 (chr) and V2 (variant ID)
stopifnot(all(snps_top %in% bim$V2))
bim_loco <- filter(bim, V1 != chr & V2 %in% snps_top)
snps_loco <- bim_loco[[2]]

## pre-computed Z matrix (FBM)
file_rds <- glue("output/top/{trait}.{n_top}.rds")
file_cols <- glue("output/top/{trait}.{n_top}.variants")

Z <- big_attach(file_rds)
cols_z <- read_lines(file_cols)
stopifnot(length(cols_z) == ncol(Z))

cols_loco <- which(cols_z %in% snps_loco) %>% sort

file_z <- tempfile(pat_temp) %>% basename
file_z_bk <- glue("{file_z}.bk")

## read phen.
phen <- read_tsv(file_phen, col_types = c("ccddddddddd"))
stopifnot(trait %in% colnames(phen))
y <- phen[[trait]] 
y_sc <- scale(y)

## load K & convert to K LOCO
K <- readRDS(glue("output/k/{trait}.{n_top}.k.rds"))
stopifnot(all(cols_z == colnames(K)))

Mloco <- length(cols_loco)
M <- ncol(Z)
K <- (M/Mloco) * K[cols_loco, cols_loco]

## fit LMM
mod <- biglmmz(y_sc, 
  Z = Z, cols = cols_loco, # LOCO: only SNPs outside chr 
  scale = TRUE, impute = FALSE, 
  backingfile = file_z,
  K = K,
  REML = TRUE, verbose = 2)

## output table
tab <- mod$ess %>% mutate(trait = trait)
print(tab)

## clean 
rm(mod)
unlink(file_z_bk)

## save
write_tsv(tab, file_out)
  
