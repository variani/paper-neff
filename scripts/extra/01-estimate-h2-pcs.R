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
ntop <- ifelse(no_args, 500, as.integer(args[2]))
file_out <- ifelse(no_args, "tmp.tsv.gz", args[3])

## parameters
file_phen <- "out/phen.sync.tsv.gz"
file_pcs <- "out/pcs.sync.tsv.gz"
file_clump <- glue("out/clump/{trait}.top.5000.clump.gz")

## testing
testing <- no_args
nrow1 <- 1e3

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

cols_subset <- which(snps_z %in% snps_top) 
snps_subset <- snps_z[cols_subset]

## subset K
M_subset <- length(snps_subset)
K <- K[snps_subset, snps_subset] / M_subset

## read phen.
phen <- read_tsv(file_phen)
stopifnot(trait %in% colnames(phen))
y <- phen[[trait]] 
y_sc <- scale(y)

pcs <- read_tsv(file_pcs)
stopifnot(all(pcs$id == phen$id))
Xpcs <- as.matrix(pcs[, -1]) # remove column 1 (ids)
Xpcs <- Xpcs[, 1:20] # subset to only 20 PCs (in total, 40 PCs)

## fit LMM
X <- cbind(matrix(1, nrow = nrow(Z), ncol = 1), Xpcs)
mod <- biglmmg(y_sc, X = X, G = Z, cols = cols_subset, K = K, verbose = 2)

## output table
tab <- tibble(trait = trait, N = length(y), M = M_subset, 
  gamma = mod$gamma, s2 = mod$s2, 
  trace_factor = mod$trace_factor, mult = mod$mult) 

print(tab)

## save
write_tsv(tab, file_out)
  

