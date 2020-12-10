### inc
library(tidyverse)
library(glue)

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
K_subset <- K[snps_subset, snps_subset] / (M_subset - 1)

## read phen.
phen <- read_tsv(file_phen)
stopifnot(trait %in% colnames(phen))
y <- phen[[trait]] 
y_sc <- scale(y)

## assoc. by LMM
# out <- mclapply(1:2, function(j) {
out <- mclapply(seq_along(cols_subset), function(j) {
  cat(" - snp", j, "/", M_subset, "\n")
  col_j <- cols_subset[j]
  snp_j <- snps_subset[j]
   
  # matrix X (covariates)
  g <- Z[, col_j]
  g_sc <- scale(g)
  X_sc <- matrix(g_sc, ncol = 1)

  # Kj
  Kj <- K_subset[-j, -j] # normalized by (M_subset - 1) earlier 

  cols_j <- cols_subset[cols_subset != col_j]

  # LMM: snp j is excluded from GRM
  mod <- biglmmg(y_sc, X = X_sc, 
    G = Z, cols = cols_j, K = Kj, 
    compute_mult = FALSE, verbose = 2)

  ## extract assoc
  coef <- tail(mod$coef, 1) # the last raw for SNP
  tab <- tibble(snp = snp_j, beta = coef$beta, se = coef$se) %>%
    mutate(zscore = beta / se,
      pval = pchisq(zscore*zscore, df = 1, lower = FALSE),
      gamma = mod$gamma, s2 = mod$s2, 
      N = length(y_sc), M = M_subset - 1) 
  print(tab)

  # clean and return
  cat(" - ", j, "/", M_subset, "finished\n")
  rm(g, g_sc, X_sc, Kj, mod); gc()

  tab
}, mc.cores = cores) 
assoc <- bind_rows(out)

assoc %>% arrange(pval) %>% print

## save
write_tsv(assoc, file_out)
  

