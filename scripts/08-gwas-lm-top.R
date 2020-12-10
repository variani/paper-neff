### inc
library(tidyverse)
library(glue)
library(data.table)
library(broom)

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
file_out <- ifelse(no_args, "tmp.tsv.gz", args[2])

## parameters
file_phen <- "out/phen.sync.tsv.gz"

## load Z
cat(" - load Z\n")
Z <- glue("out/fbm/{trait}.rds") %>% big_attach
snps_z <- glue("out/fbm/{trait}.variants") %>% read_lines

snps_test <- snps_z
cols_test <- seq_along(snps_test)

## read phen.
phen <- file_phen %>% fread %>% as_tibble
stopifnot(trait %in% colnames(phen))
y <- phen[[trait]] 
y_sc <- scale(y)

## assoc. by LMM
X0 <- matrix(1, nrow = nrow(Z), ncol = 1)

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
  X_sc <- cbind(X0, matrix(g_sc, ncol = 1))

  # LMM: snp j is excluded from GRM
  mod <- lm(y_sc ~ X_sc) 
  s2 <- sigma(mod)^2
  coef <- broom::tidy(mod) %>% tail(1) # the last covariat is g

  tab <- tibble(snp = snp_j, beta = coef$estimate, se = coef$std.error) %>%
    mutate(zscore = beta / se,
      pval = pchisq(zscore*zscore, df = 1, lower = FALSE),
      s2 = s2, N = length(y_sc))
  print(tab)

  # clean and return
  cat(" - ", j, "/", M_test, "finished\n")
  rm(g, g_sc, X_sc); gc()

  tab
}, mc.cores = cores) 
# }
assoc <- bind_rows(out)

assoc %>% arrange(pval) %>% print

## save
write_tsv(assoc, file_out)
  

