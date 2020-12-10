### inc
library(tidyverse)
library(glue)
library(data.table)

library(BEDMatrix)
library(bigsnpr)

library(devtools)
load_all("~/git/variani/biglmmz/")

## arguments
trait <- "hip"
ntop <- 500

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

## read phen.
phen <- file_phen %>% fread %>% as_tibble
stopifnot(trait %in% colnames(phen))
y <- phen[[trait]] 

pcs <- file_pcs %>% fread %>% as_tibble

## compute corr.
G <- Z[, seq(ntop)]
cor(pcs$PC1, G) %>% as.numeric %>% sort %>% head %>% round(2)
cor(pcs$PC2, G) %>% as.numeric %>% sort %>% head %>% round(2)
