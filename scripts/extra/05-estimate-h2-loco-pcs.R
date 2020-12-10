### inc
library(tidyverse)
library(data.table)
library(glue)

library(BEDMatrix)
library(bigsnpr)

library(devtools)
load_all("~/git/variani/biglmmz/")

## arguments
args <- commandArgs(trailingOnly = TRUE)
no_args <- is.na(args[1])

trait <- ifelse(no_args, "height", args[1])
ntop <- ifelse(no_args, 500, as.integer(args[2]))
chr <- ifelse(no_args, 20, as.integer(args[3]))
file_out <- ifelse(no_args, "tmp.tsv.gz", args[4])

## parameters
file_phen <- "out/phen.sync.tsv.gz"
file_pcs <- "out/pcs.sync.tsv.gz"
file_clump <- glue("out/clump/{trait}.top.5000.clump.gz")

## read top variants
cat(" - read clump\n")
clump <- read_tsv(file_clump) %>% arrange(P)
if(nrow(clump) < ntop) {
  stop(glue("#rows in clump ({nrow(clump)}) < #top ({ntop})"))
}

clump <- head(clump, ntop)
snps_top <- clump$SNP

## load Z
cat(" - load Z\n")
Z <- glue("out/fbm/{trait}.rds") %>% big_attach
snps_z <- glue("out/fbm/{trait}.variants") %>% read_lines

## load ZtZ
cat(" - load ZtZ\n")
K <- glue("out/ztz/{trait}.rds") %>% readRDS
snps_k <- colnames(K)

## check snps
stopifnot(all(snps_k %in% snps_z))
stopifnot(all(snps_top %in% snps_z))

## LOCO
snps_chr <- with(clump, SNP[CHR == chr])
snps_loco <- snps_top[!(snps_top %in% snps_chr)]

cols_loco <- which(snps_z %in% snps_loco) 
snps_loco <- snps_z[cols_loco]

## subset K
M_loco <- length(snps_loco)
K_loco <- K[snps_loco, snps_loco] / M_loco

## read phen.
phen <- read_tsv(file_phen)
stopifnot(trait %in% colnames(phen))
y <- phen[[trait]] 
y_sc <- scale(y)

pcs <- file_pcs %>% fread %>% as_tibble
stopifnot(all(pcs$id == phen$id))
Xpcs <- as.matrix(pcs[, -1]) # remove column 1 (ids)
Xpcs <- Xpcs[, 1:20] # subset to only 20 PCs (in total, 40 PCs)

## fit LMM
Xcov <- cbind(matrix(1, nrow = nrow(Z), ncol = 1), Xpcs)
mod <- biglmmg(y_sc, Xcov, G = Z, cols = cols_loco, K = K_loco, verbose = 2)

## output table
tab <- tibble(trait = trait, N = length(y), M = M_loco, 
  gamma = mod$gamma, s2 = mod$s2, 
  trace_factor = mod$trace_factor, mult = mod$mult) 

print(tab)

## save
write_tsv(tab, file_out)
  


