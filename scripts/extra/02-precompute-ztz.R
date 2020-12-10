library(tidyverse)
library(glue)
library(data.table)

library(bigsnpr)

## args
args <- commandArgs(trailingOnly = TRUE)
no_args <- is.na(args[1])

trait <- ifelse(no_args, "height", args[1])
file_out <- ifelse(no_args, "tmp.rds", args[2])

## parameters
ntop <- 5000

## load data
cat(" - load Z\n")
file_rds <- glue("out/fbm/{trait}.rds")
Z <- big_attach(file_rds)

snps_z <- glue("out/fbm/{trait}.variants") %>% read_lines
stopifnot(length(snps_z) == ncol(Z))

## clump
cat(" - read clump\n")
file_clump <- glue("out/clump/{trait}.top.5000.clump.gz")
clump <- read_tsv(file_clump) %>% arrange(P)
if(nrow(clump) < ntop) {
  stop(glue("#rows in clump ({nrow(clump)}) < #top ({ntop})"))
}

clump <- head(clump, ntop)
snps_top <- clump$SNP

## check snps
stopifnot(all(snps_top %in% snps_z))

cols_subset <- which(snps_z %in% snps_top) 
snps_subset <- snps_z[cols_subset]

## compute K = Z'Z, where Z is a matrix of scaled genotypes
cat(" - compute K\n")
K <- big_crossprodSelf(Z, fun.scaling = big_scale(), ind.col = cols_subset)[]
colnames(K) <- snps_subset
rownames(K) <- snps_subset

## save
saveRDS(K, file_out)
cat(" - ZtZ done!\n")


