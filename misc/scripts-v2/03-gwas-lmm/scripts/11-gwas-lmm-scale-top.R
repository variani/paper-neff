library(tidyverse)
library(glue)
library(data.table)

library(bigsnpr)

library(devtools)
load_all("~/git/variani/matlm") # gwas assoc. model
load_all("~/git/variani/bigcov") # split bed by batches
load_all("~/git/variani/biglmmz") # impute/scale large matrices 

library(BEDMatrix) # read bed from file & avoid loading into RAM

## parallel
library(parallel)
cores <- max(1, detectCores() - 1)

## args
args <- commandArgs(trailingOnly = TRUE)
no_args <- is.na(args[1])

trait <- ifelse(no_args, "height", args[1])
ntop <- ifelse(no_args, 500, as.integer(args[2]))
file_out <- ifelse(no_args, "tmp.tsv.gz", args[3])

## top genotypes 
file_rds <- glue("output/top/{trait}.{ntop}.rds")
Z <- big_attach(file_rds)

snps_z <- glue("output/top/{trait}.{ntop}.variants") %>% read_lines
stopifnot(length(snps_z) == ncol(Z))

## Kj matrix
Kj <- glue("output/kj/{trait}.{ntop}.kj.rds") %>% readRDS
stopifnot(all(colnames(Kj) == snps_z))

# extract estimates of model parameters (variance components)
tab <- glue("output/h2/{trait}.{ntop}.tsv.gz") %>% read_tsv
gamma <- tab$h2
s2 <- tab$s2
comp <- s2 * c(gamma, 1 - gamma)

## phen: load phen sync. with bed ids
file_phen <- "output/phen.bed.tsv.gz"
phen <- read_tsv(file_phen, col_types = c("ccddddddddd"))
stopifnot(trait %in% colnames(phen))
y <- phen[[trait]] 

y_sc <- scale(y)

## assoc. by LMM
snpnames <- snps_z
M <- ncol(Z)
N <- nrow(Z)

out <- mclapply(seq_along(snpnames), function(j) {
  cat(" - snp", j, "/", M, "\n")
  snp <- snpnames[j]
   
  # matrix X (covariates)
  g <- Z[, j]
  g[is.na(g)] <- 0
  g_sc <- scale(g)

  # LMM: snp j is excluded from GRM
  cols <- seq(ncol(Z))
  cols <- cols[cols != j]

  X_sc <- matrix(g_sc, ncol = 1)
  XVt <- biglr_cprodMatInv2(comp, Z, cols = cols, 
    X = X_sc,
    K = Kj[-j, -j], transpose = TRUE) # Vi' X
  XVX <- colSums(XVt * X_sc) 
  b <- as.numeric(crossprod(XVt, y_sc)) / XVX
  se <- 1 / sqrt(XVX)
  assoc <- tibble(predictor = snp, beta = b, se = se) %>%
    mutate(zscore = beta / se, 
      pval = pchisq(zscore*zscore, df = 1, lower = FALSE))

  cat(" - ", j, "/", M, "finished\n")
  rm(g, g_sc, X_sc, XVt, XVX); gc()

  assoc
}, mc.cores = cores) 
assoc <- bind_rows(out)

## save
write_tsv(assoc, file_out)
