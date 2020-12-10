### inc
library(tidyverse)
library(glue)
library(data.table)

library(BEDMatrix)
library(bigsnpr)

library(devtools)
load_all("~/git/variani/biglmmz/")

## parallel
library(parallel)
cores <- detectCores() - 1
# cores <- min(cores, 10)

## arguments
args <- commandArgs(trailingOnly = TRUE)
no_args <- is.na(args[1])

trait <- ifelse(no_args, "height", args[1])
ntop <- ifelse(no_args, 1000, as.integer(args[2]))
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
snps_test <- snps_top[snps_top %in% snps_chr]

# re-order snps
cols_loco <- which(snps_z %in% snps_loco) 
snps_loco <- snps_z[cols_loco]

cols_test <- which(snps_z %in% snps_test) 
snps_test <- snps_z[cols_test]

## subset K
M_loco <- length(snps_loco)
K_loco <- K[snps_loco, snps_loco] / M_loco

## read gamma, s2 from h2 model (LOCO)
tab <- glue("out/h2-loco_pcs/{ntop}/{trait}.{chr}.tsv.gz") %>% read_tsv
# tab <- tibble(gamma = 0.5, s2 = 1)
gamma <- tab$gamma
s2 <- tab$s2

## read phen.
phen <- file_phen %>% fread %>% as_tibble
stopifnot(trait %in% colnames(phen))
y <- phen[[trait]] 
y_sc <- scale(y)

pcs <- file_pcs %>% fread %>% as_tibble
stopifnot(all(pcs$id == phen$id))
Xpcs <- as.matrix(pcs[, -1]) # remove column 1 (ids)
Xpcs <- Xpcs[, 1:20] # subset to only 20 PCs (in total, 40 PCs)

## assoc. by LMM
Xcov <- cbind(matrix(1, nrow = nrow(Z), ncol = 1), Xpcs)

M_test <- length(cols_test)
# for(i in 1:1) {
# out <- mclapply(1:2, function(j) {
out <- mclapply(seq_along(cols_test), function(j) {
  cat(" - snp", j, "/", M_test, "\n")
  col_j <- cols_test[j]
  snp_j <- snps_test[j]
   
  # matrix X (covariates)
  g <- Z[, col_j]
  g_sc <- scale(g)
  X_sc <- cbind(Xcov, matrix(g_sc, ncol = 1))

  # LMM: snp j is excluded from GRM
  est <- biglr_fixef_grm(gamma, y_sc, X_sc,
    Z, cols_loco,
    s2 = s2, K = K_loco)
  coef <- data.frame(beta = est$b, se = sqrt(diag(est$bcov)))
  coef <- tail(coef, 1) # the last covariat is g

  tab <- tibble(snp = snp_j, beta = coef$beta, se = coef$se) %>%
    mutate(zscore = beta / se,
      pval = pchisq(zscore*zscore, df = 1, lower = FALSE),
      gamma = gamma, s2 = s2, 
      N = length(y_sc), M = M_loco)
  print(tab)

  # clean and return
  cat(" - ", j, "/", M_test, "finished\n")
  rm(g, g_sc, X_sc); gc()

  tab
}, mc.cores = cores) 
# }
assoc <- bind_rows(out)

assoc$chr <- chr

assoc %>% arrange(pval) %>% print

## save
write_tsv(assoc, file_out)
  

