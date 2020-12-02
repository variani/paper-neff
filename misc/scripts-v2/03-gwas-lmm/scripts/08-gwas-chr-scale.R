library(tidyverse)
library(glue)
library(data.table)

library(devtools)
load_all("~/git/variani/bigcov") # split bed by batches
load_all("~/git/variani/biglmmz") # impute/scale large matrices 

library(BEDMatrix) # read bed from file & avoid loading into RAM

## parallel
library(parallel)
cores <- detectCores() - 1

## args
args <- commandArgs(trailingOnly = TRUE)
no_args <- is.na(args[1])

trait <- ifelse(no_args, "height", args[1])
ntop <- ifelse(no_args, 500, as.integer(args[2]))
chr <- ifelse(no_args, 22, as.integer(args[3]))
file_out <- ifelse(no_args, "tmp.tsv.gz", args[4])

## testing
testing <- no_args

## parallel
library(parallel)
cores <- max(1, detectCores() - 1)

## par
nsnps_batch <- 50

## data
file_bed <- "output/gen.bed"
file_bim <- gsub(".bed$", ".bim", file_bed)
file_phen <- "output/phen.bed.tsv.gz"

## phen: load phen sync. with bed ids
phen <- read_tsv(file_phen, col_types = c("ccdddddd"))
stopifnot(trait %in% colnames(phen))
y <- phen[[trait]] 

y_sc <- scale(y)

## read Z
file_rds <- glue("output/top/{trait}.{ntop}.rds")
Z <- big_attach(file_rds)

snps_z <- glue("output/top/{trait}.{ntop}.variants") %>% read_lines
stopifnot(length(snps_z) == ncol(Z))

ids_z <- glue("output/top/{trait}.{ntop}.ids") %>% read_lines
stopifnot(length(ids_z) == nrow(Z))
stopifnot(all(ids_z == phen$IID))

# filter top variants by chr
bim <- fread(file_bim, select = 1:2) %>% as_tibble # columns V1 (chr) and V2 (variant ID)
stopifnot(all(snps_z %in% bim$V2))
bim_loco <- filter(bim, V1 != chr & V2 %in% snps_z)
snps_loco <- bim_loco[[2]]

snps_test <- filter(bim, V1 == chr)[[2]]
if(testing) {
  snps_test <- head(snps_test, 100)
}

## bed
snps <- read_lines("output/gen.variants")
bed <- BEDMatrix(file_bed, p = length(snps))
colnames(bed) <- snps

# extract estimates of model parameters (variance components)
# tab <- glue("output/h2/{trait}.{ntop}.tsv.gz") %>% read_tsv
tab <- glue("output/h2_loco/{trait}.{ntop}.{chr}.tsv.gz") %>% read_tsv
gamma <- tab$h2
s2 <- tab$s2
comp <- s2 * c(gamma, 1 - gamma)

## organize tested snps (chr) in batches (batches_test)
bdat <- bigcov::bigdat(bed, batch_size = nsnps_batch)
annot <- tibble(V2 = colnames(bed)) %>% 
  mutate(index = seq(n())) %>% left_join(bim, by = "V2") 
annot_test <- filter(annot, V1 == chr)
min_index <- min(annot_test$index)
max_index <- max(annot_test$index)

batches_test <- filter(bdat$batches, 
  (beg >= min_index & beg <= max_index) |
  (end >= min_index & end <= max_index))
n_batches_test <- batches_test %>% mutate(size = end - beg + 1) %>% .[["size"]] %>% sum
n_chr <- bim %>% filter(V1 == chr) %>% nrow
stopifnot(n_batches_test >= n_chr)

B <- batches_test$batch

# pre-compute
cat(" - pre-compute K\n")
cols <- which(snps_z %in% snps_loco) %>% sort
M <- length(cols)
stopifnot(M > 0)
K <- big_crossprodSelf(Z, fun.scaling = big_scale2(M = M), ind.col = cols)[]

out <- mclapply(seq_along(B), function(bi) {
  cat(" -", bi, "/", length(B), "\n")
  b <- B[bi]
  X <- bigdat_batch(bdat, b)

  ind_test <- (colnames(X) %in% snps_test)
  if(!any(ind_test)) {
    assoc <- tibble()
    return(assoc)
  }

  X <- X[, ind_test, drop = FALSE]
  X_names <- colnames(X)
  X[is.na(X)] <- 0
  X_sc <- scale(X)

  XVt <- biglr_cprodMatInv2(comp, Z, cols = cols, X = X_sc, K = K, transpose = TRUE) # Vi' X
  XVX <- colSums(XVt * X_sc) 
  b <- as.numeric(crossprod(XVt, y_sc)) / XVX
  se <- 1 / sqrt(XVX)
  assoc <- tibble(predictor = X_names, beta = b, se = se) %>%
    mutate(zscore = beta / se, 
      pval = pchisq(zscore*zscore, df = 1, lower = FALSE),
      gamma = gamma, s2 = s2)

  n_hits <- filter(assoc, pval < 5e-8) %>% nrow
  if(n_hits > 0) {
    cat(" -- #hits:", n_hits, "/", nrow(assoc), "\n")
  }

  cat(" - ", bi, "finished\n")
  rm(X, X_sc, XVt, XVX); gc()

  assoc
}, mc.cores = cores)

assoc <- bind_rows(out)
stopifnot(nrow(assoc) == length(snps_test))

## save
write_tsv(assoc, file_out)
