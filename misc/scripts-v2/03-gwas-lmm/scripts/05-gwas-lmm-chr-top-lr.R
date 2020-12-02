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
cores <- min(cores, 10)

## parameters
thr_lr <- 5e-8
pat_temp <- "file-05-top-lr"

## args
args <- commandArgs(trailingOnly = TRUE)
no_args <- is.na(args[1])

trait <- ifelse(no_args, "height", args[1])
n_top <- ifelse(no_args, 500, as.integer(args[2]))
chr <- ifelse(no_args, 22, as.integer(args[3]))
file_out <- ifelse(no_args, "tmp.tsv.gz", args[4])

## testing
testing <- no_args

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

## filter top variants by chr & by LR GWAS
bim <- fread(file_bim, select = 1:2) %>% as_tibble # columns V1 (chr) and V2 (variant ID)
stopifnot(all(snps_top %in% bim$V2))
bim_loco <- filter(bim, V1 != chr & V2 %in% snps_top)
snps_loco <- bim_loco[[2]]

snps_test_chr <- filter(bim, V1 == chr)[[2]]

assoc_lr <- glue("output/gwas-lr/{trait}.tsv.gz") %>% read_tsv %>% filter(pval < thr_lr)
snps_lr <- assoc_lr$predictor
snps_test_lr <- snps_test_chr[snps_test_chr %in% snps_lr]
if(length(snps_test_lr) == 0) {
  snps_test_lr <- arrange(assoc_lr, pval) %>% .[[1]]
}

snps_test <- snps_test_lr

if(testing) {
  snps_test <- head(snps_test, 10)
}

## bed
snps <- read_lines("output/gen.variants")
bed <- BEDMatrix(file_bed, p = length(snps))
colnames(bed) <- snps

## bigSNP
file_b <- tempfile(glue("{pat_temp}-readBed2-")) %>% basename
file_bk <- glue("{file_b}.bk")
stopifnot(!file.exists(file_bk))

cols_test <- which(snps %in% snps_test) %>% sort
file_rds <- snp_readBed2(file_bed, backingfile = file_b, ind.col = cols_test)
bigsnp <- snp_attach(file_rds)
gen <- bigsnp$genotypes

## phen: load phen sync. with bed ids
phen <- read_tsv(file_phen, col_types = c("ccdddddd"))
stopifnot(trait %in% colnames(phen))
y <- phen[[trait]] 

## fit the null model  & get estimate of h2
cols <- which(snps %in% snps_loco)
if(testing) {
  nrow1 <- 1e3
  y <- head(y, nrow1)
  Zg <- bed[seq(nrow1), cols]
} else {
  Zg <- bed[, cols]
}

file_Zmod <- tempfile(glue("{pat_temp}-zmod-")) %>% basename()
mod <- biglmmz(y, Z = Zg, impute = TRUE, scale = TRUE, 
  backingfile = file_Zmod,
  store_mat = TRUE, verbose = 2)

# extract estimates of model parameters (variance components)
gamma <- mod$gamma
s2 <- mod$s2

file_grm <- tempfile(glue("{pat_temp}-grm-")) %>% basename
Zgrm <- big_copy(mod$Z, backingfile = file_grm, is_read_only = TRUE)
# mod2 <- biglmmz(yc, Z = Zgrm, verbose = 2)

# pre-compute
K <- mod$K

# remove model 
rm(mod); rm(Zg); file.remove(paste0(file_Zmod, ".bk"))

gnames <- bigsnp$map[[2]]
n <- ncol(gen)
out <- mclapply(seq(n), function(i) {
  cat(" - ", i, "/", n,  "started\n")
  
  g <- gen[, i]
  if(testing) {
    g <- g[seq(nrow1)]
  } 

  g <- scale(g, center = TRUE, scale = TRUE) 
  # impute missing
  g[is.na(g)] <- 0

  X0 <- matrix(1, nrow = length(g), ncol = 1)
  X <- cbind(X0, g)

  est <- biglmmz:::biglr_fixef(gamma = gamma, s2 = s2,
    y = y, Xmat = X, Z = Zgrm, K = K, REML = TRUE)
  coef <- data.frame(beta = est$b, se = sqrt(diag(est$bcov)))
  coef <- tail(coef, 1) # SNP goes last
  
  assoc <- tibble(predictor = gnames[i], chr = chr, 
      beta = coef$beta, se = coef$se) %>%
    mutate(zscore = beta / se) %>%
    mutate(
      pval = pchisq(zscore*zscore, df = 1, lower = FALSE),
      gamma = gamma, s2 = s2)

  cat(" - ", i, "/", n,  "finished\n")
  rm(X); gc()

  assoc
}, mc.cores = cores) 

assoc <- bind_rows(out)
stopifnot(nrow(assoc) == length(snps_test))

## clean
file.remove(c(paste0(file_grm, ".bk"), file_rds, file_bk))

## print info
assoc %>% arrange(pval) %>% print

## save
write_tsv(assoc, file_out)
