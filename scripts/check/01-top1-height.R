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

cols_subset <- which(snps_z %in% snps_top) 
snps_subset <- snps_z[cols_subset]

## subset K
M_subset <- length(snps_subset)
K_subset <- K[snps_subset, snps_subset] 

## read phen.
phen <- read_tsv(file_phen)
stopifnot(trait %in% colnames(phen))
y <- phen[[trait]] 
y_sc <- scale(y)

## null model
K <- K_subset / M_subset
m0 <- biglmmg(y_sc, 
  G = Z, cols = cols_subset, K = K, 
  compute_mult = FALSE, verbose = 2)

out <- list(h2_m0 = m0$gamma)

## snp top #1
j <- 1
col_j <- cols_subset[j]
snp_j <- snps_subset[j]
 
g <- Z[, col_j]
g_sc <- scale(g)
X_sc <- matrix(g_sc, ncol = 1)

K <- K_subset[-j, -j] / (M_subset - 1) 
cols_j <- cols_subset[cols_subset != col_j]
m1 <- biglmmg(y_sc, X = X_sc, 
  G = Z, cols = cols_j, K = K, 
  compute_mult = FALSE, verbose = 2)

X1_sc <- cbind(matrix(1, nrow = nrow(Z), ncol = 1), X_sc)
m2 <- biglmmg(y_sc, X = X1_sc, 
  G = Z, cols = cols_j, K = K, 
  compute_mult = FALSE, verbose = 2)

## load Z into RAM
Ze <- Z[, cols_j]

X0 <- matrix(1, nrow = nrow(Z), ncol = 1)
X <- cbind(X0, Z[, col_j])

m3 <- biglmmz(y_sc, X, Ze, scale = TRUE, verbose = 2)

# > cor(g, Ze) %>% as.numeric %>% sort %>% head %>% round(2)
# [1] -0.33 -0.14 -0.08 -0.08 -0.07  0.00

## LOCO
snps_top <- clump$SNP
snps_chr <- with(clump, SNP[CHR == 20])

snps_loco <- snps_top[!(snps_top %in% snps_chr)]

cols_loco <- which(snps_z %in% snps_loco) 
snps_loco <- snps_z[cols_loco]

## subset K
K <- glue("out/ztz/{trait}.rds") %>% readRDS
M_loco <- length(snps_loco)
K_loco <- K[snps_loco, snps_loco] / M_loco

m4 <- biglmmg(y_sc, X = X_sc, 
  G = Z, cols = cols_loco, K = K_loco, 
  compute_mult = FALSE, verbose = 2)

# > cor(g, Z[, cols_loco]) %>% as.numeric %>% sort %>% head %>% round(3)
# [1] -0.005 -0.004 -0.004 -0.004 -0.004 -0.004
