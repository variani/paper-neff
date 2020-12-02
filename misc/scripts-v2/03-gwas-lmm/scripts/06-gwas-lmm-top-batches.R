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
cores <- detectCores() - 1
# cores <- min(cores, 20)

## par
pat_temp <- "file-06-top-"

## args
args <- commandArgs(trailingOnly = TRUE)
no_args <- is.na(args[1])

trait <- ifelse(no_args, "height", args[1])
n_top <- ifelse(no_args, 500, as.integer(args[2]))
file_out <- ifelse(no_args, "tmp.tsv.gz", args[3])

## testing
testing <- no_args
nrow1 <- 1e3

## data
file_bed <- "output/gen.bed"
file_bim <- gsub(".bed$", ".bim", file_bed)
file_phen <- "output/phen.bed.tsv.gz"
file_clump <- glue("output/clump/{trait}.tsv.gz")

## read top variants
clump <- read_tsv(file_clump) %>% arrange(pval) 
if(nrow(clump) < n_top) {
  warning(glue("#rows in clump ({nrow(clump)}) < #top ({n_top})"))
}

clump <- head(clump, n_top)
snps_top <- clump$snp

## bed
snps <- read_lines("output/gen.variants")
ids <- read_lines("output/gen.ids")
bed <- BEDMatrix(file_bed, n = length(ids), p = length(snps))
colnames(bed) <- snps

## Z matrix (FBM)
file_b <- tempfile(paste0(pat_temp, "readBed2-")) %>% basename
file_bk <- glue("{file_b}.bk")
cols_top <- which(snps %in% snps_top) %>% sort
if(testing) {
  file_rds <- snp_readBed2(file_bed, backingfile = file_b, ind.row = seq(nrow1), ind.col = cols_top)
} else {
  file_rds <- snp_readBed2(file_bed, backingfile = file_b, ind.col = cols_top)
}

bigsnp <- snp_attach(file_rds)
snpnames_Z <- bigsnp$map[[2]]
Z0 <- bigsnp$genotypes

file_z <- tempfile(paste0(pat_temp, "Z-")) %>% basename
file_z_bk <- glue("{file_z}.bk")
Z <- big_copy(Z0, type = "double", backingfile = file_z)
# clean bigsnp
unlink(c(file_bk, file_rds))

scale_Z(Z, impute = TRUE, M = ncol(Z) - 1)

## phen: load phen sync. with bed ids
phen <- read_tsv(file_phen, col_types = c("ccdddddd"))
stopifnot(trait %in% colnames(phen))
y <- phen[[trait]] 
if(testing) {
  y <- head(y, nrow1)
}

# pre-compute K = Z'Z
K <- crossprod(Z)

## assoc. by LMM
snpnames <- snpnames_Z
M <- ncol(Z)
N <- nrow(Z)

if(testing) {
  snpnames <- head(snpnames, 5)
}

dir_batch <- file.path(dirname(file_out), "batches-06/")
dir.create(dir_batch, show = FALSE)
fun_file_batch <- function(trait, snp) glue("{dir_batch}/{trait}-{snp}.tsv.gz")

out <- mclapply(seq_along(snpnames), function(j) {
  cat(" - snp", j, "/", length(snpnames), "\n")
  snp <- snpnames[j]
  file_batch <- fun_file_batch(trait, snp)

  # early stop
  if(file.exists(file_batch)) {
    cat("  -- batch pre-computed for  snp", j, "/", length(snpnames), "\n")
    return(tibble())
  }
   
  # matrix X (covariates)
  g <- Z[, j]
  g <- scale(g, center = TRUE, scale = TRUE)

  X0 <- matrix(1, nrow = length(g), ncol = 1)
  X <- cbind(X0, g)

  # LMM: snp j is excluded from GRM
  cols <- seq(ncol(Z))
  cols <- cols[cols != j]

  file_j <- tempfile(paste0(pat_temp, "snp-j-")) %>% basename
  file_j_bk <- glue("{file_j}.bk")
  mod <- biglmmz(y, X = X, Z = Z, cols = cols, K = K[-j, -j],
    backingfile = file_j, verbose = 2)
  coef <- tail(mod$coef, 1) # the last raw for SNP
  file.remove(file_j_bk)

  # est <- biglmmz:::biglr_fixef(gamma = gamma, s2 = s2,
  #   y = y, Xmat = X, Z = Zgrm, K = K, REML = TRUE)
  # coef <- data.frame(beta = est$b, se = sqrt(diag(est$bcov)))
  # coef <- tail(coef, 1) # SNP goes last

  tab <- tibble(snp = snp, beta = coef$beta, se = coef$se) %>%
    mutate(zscore = beta / se,
      pval = pchisq(zscore*zscore, df = 1, lower = FALSE),
      h2 = mod$gamma, s2 = mod$s2, N = N, M = M) 
  print(tab)
  write_tsv(tab, file_batch)

  invisible()
}, mc.cores = cores) 

## merge results
files_batch <- fun_file_batch(trait, snpnames)
stopifnot(all(file.exists(files_batch)))

assoc <- mclapply(files_batch, read_tsv, mc.cores = cores) %>%
  bind_rows

## clean
file.remove(file_z_bk)

## save
write_tsv(assoc, file_out)
