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
ntop <- ifelse(no_args, 500, as.integer(args[2]))
file_out <- ifelse(no_args, "tmp.rds", args[3])

## top genotypes 
file_rds <- glue("output/top/{trait}.{ntop}.rds")
Z <- big_attach(file_rds)

snps_z <- glue("output/top/{trait}.{ntop}.variants") %>% read_lines
stopifnot(length(snps_z) == ncol(Z))

cat(" - compute K\n")
M <- ncol(Z)
K <- big_crossprodSelf(Z, fun.scaling = big_scale2(M = M))[]
colnames(K) <- snps_z
rownames(K) <- snps_z

## save
saveRDS(K, file_out)
cat(" - done!\n")

