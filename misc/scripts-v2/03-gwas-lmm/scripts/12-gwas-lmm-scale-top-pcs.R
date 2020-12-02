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

## PCs
cat(" - read PCs\n")
file_pcs <- "output/pcs.bed.tsv.gz"
pcs <- read_tsv(file_pcs)
C <- select(pcs, starts_with("PC")) %>% as.matrix
stopifnot(ncol(C) == 10)

svd <- svd(C, nv = 0)
eigval_sc <- svd$d/(sqrt(nrow(C)) + sqrt(ncol(C)) - 1)
svd_ncomp <- sum(eigval_sc > 1e-4)
U <- svd$u[, seq(svd_ncomp), drop = FALSE]

## top genotypes 
file_rds <- glue("output/top/{trait}.{ntop}.rds")
Z <- big_attach(file_rds)

snps_z <- glue("output/top/{trait}.{ntop}.variants") %>% read_lines
stopifnot(length(snps_z) == ncol(Z))

## K matrix
cat(" - load K\n")
K <- readRDS(glue("output/k/{trait}.{ntop}.k.rds"))
stopifnot(all(snps_z == colnames(K)))

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

stopifnot(nrow(C) == length(y))
y_orth <- as.numeric(y - U %*% crossprod(U, y))

y_sc <- scale(y_orth)

## assoc. by LMM
snpnames <- snps_z
M <- ncol(Z)
N <- nrow(Z)

out <- mclapply(seq_along(snpnames), function(j) {
  cat(" - snp", j, "/", M, "\n")
  snp <- snpnames[j]
   
  # matrix X (covariates)
  X <- Z[, j, drop = FALSE]
  # impute
  X[is.na(X)] <- 0
  # project PCs out
  X <- X - U %*% crossprod(U, X)
  # scale
  X_sc <- scale(X)

  # LMM: snp j is excluded from GRM
  cols <- seq(ncol(Z))
  cols <- cols[cols != j]

  Kj <- (M/(M-1)) * K[cols, cols]

  XVt <- biglr_cprodMatInv2(comp, Z, cols = cols, X = X_sc, K = Kj, transpose = TRUE) # Vi' X
  XVX <- colSums(XVt * X_sc) 
  b <- as.numeric(crossprod(XVt, y_sc)) / XVX
  se <- 1 / sqrt(XVX)
  assoc <- tibble(predictor = snp, beta = b, se = se) %>%
    mutate(zscore = beta / se, 
      pval = pchisq(zscore*zscore, df = 1, lower = FALSE))

  cat(" - ", j, "/", M, "finished\n")
  rm(Kj, X, X_sc, XVt, XVX); gc()

  assoc
}, mc.cores = cores) 
assoc <- bind_rows(out)

## save
write_tsv(assoc, file_out)
