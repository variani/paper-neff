library(tidyverse)
library(glue)
library(data.table)

library(devtools)
load_all("~/git/variani/matlm") # gwas assoc. model
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
n_top <- ifelse(no_args, 500, as.integer(args[2]))
chr <- ifelse(no_args, 22, as.integer(args[3]))
file_out <- ifelse(no_args, "tmp.tsv.gz", args[4])

## testing
testing <- no_args

## parallel
cores <- parallel::detectCores()
cores <- ifelse(cores > 1, cores - 1, cores)
parallel <- (cores > 1)
if(parallel) { doParallel::registerDoParallel(cores = cores) }

## par
nsnps_batch <- 20

## data
file_bed <- "output/gen.bed"
file_bim <- gsub(".bed$", ".bim", file_bed)
file_phen <- "output/phen.bed.tsv.gz"
file_clump <- glue("output/clump/{trait}.tsv.gz")

## read top variants
clump <- read_tsv(file_clump) %>% arrange(pval) 
if(nrow(clump) < n_top) {
  stop(glue("#rows in clump ({nrow(clump)}) < #top ({n_top})"))
}

clump <- head(clump, n_top)
snps_top <- clump$snp

# filter top variants by chr
bim <- fread(file_bim, select = 1:2) %>% as_tibble # columns V1 (chr) and V2 (variant ID)
stopifnot(all(snps_top %in% bim$V2))
bim_loco <- filter(bim, V1 != chr & V2 %in% snps_top)
snps_loco <- bim_loco[[2]]

snps_test <- filter(bim, V1 == chr)[[2]]
if(testing) {
  snps_test <- head(snps_test, 100)
}

## bed
snps <- read_lines("output/gen.variants")
bed <- BEDMatrix(file_bed, p = length(snps))
colnames(bed) <- snps

## phen: load phen sync. with bed ids
phen <- read_tsv(file_phen, col_types = c("ccdddddd"))
stopifnot(trait %in% colnames(phen))
y <- phen[[trait]] 

y_sc <- scale(y)

## fit the null model (no snpts) & get estimate of h2
cols <- which(snps %in% snps_loco)
if(testing) {
  nrow1 <- 1e3
  y <- head(y, nrow1)
  Zg <- bed[seq(nrow1), cols]
} else {
  Zg <- bed[, cols]
}

mod <- biglmmz(y, Z = Zg, impute = TRUE, scale = TRUE, 
  store_mat = TRUE, verbose = 2)

# extract estimates of model parameters (variance components)
gamma <- mod$gamma
s2 <- mod$s2

file_grm <- tempfile()
Zgrm <- big_copy(mod$Z, backingfile = file_grm, is_read_only = TRUE)

# remove model 
rm(mod); rm(Zg)

## organize snps
bdat <- bigdat(bed, batch_size = nsnps_batch)
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
K <- crossprod(Zgrm)

out <- mclapply(seq_along(B), function(bi) {
  cat(" -", bi, "/", length(B), "\n")
  b <- B[bi]
  X <- bigdat_batch(bdat, b)
  if(testing) {
    X <- X[seq(nrow1), ]
  } 

  ind_test <- (colnames(X) %in% snps_test)
  if(!any(ind_test)) {
    assoc <- tibble()
    return(assoc)
  }

  X <- X[, ind_test, drop = FALSE]

  X_names <- colnames(X)
  assoc <- lapply(seq_along(X_names), function(i) {
    # cat("  --", i, "/", ncocl(Xc), "\n")
    
    # scale genotypes: need std. beta
    g <- scale(X[, i], center = TRUE, scale = TRUE) 
    # impute missing
    g[is.na(g)] <- 0

    X0 <- matrix(1, nrow = nrow(X), ncol = 1)
    X <- cbind(X0, g)

    est <- biglmmz:::biglr_fixef(gamma = gamma, s2 = s2,
      y = y, Xmat = X, Z = Zgrm, K = K, REML = TRUE)
    coef <- data.frame(beta = est$b, se = sqrt(diag(est$bcov)))
    coef <- tail(coef, 1) # SNP goes last
    
    tibble(predictor = X_names[i], chr = chr, 
        beta = coef$beta, se = coef$se) %>%
      mutate(zscore = beta / se)
  }) %>% bind_rows
  
  assoc <- mutate(assoc, 
    pval = pchisq(zscore*zscore, df = 1, lower = FALSE),
    gamma = gamma, s2 = s2)

  n_hits <- filter(assoc, pval < 5e-8) %>% nrow
  if(n_hits > 0) {
    cat(" -- #hits:", n_hits, "/", nrow(assoc), "\n")
  }

  cat(" - ", bi, "finished\n")
  rm(X); gc()

  assoc
}, mc.cores = cores) 

assoc <- bind_rows(out)
stopifnot(nrow(assoc) == length(snps_test))

## clean
unlink(file_grm)

## save
write_tsv(assoc, file_out)
