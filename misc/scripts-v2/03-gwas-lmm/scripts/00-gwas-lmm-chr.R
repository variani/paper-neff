### inc
library(magrittr)
library(dplyr)
library(readr)

library(devtools)
load_all("~/git/hemostat/ukbjass/")
load_all("~/git/hemostat/ukbpower/biglmm/")

library(gaston)
options(gaston.auto.set.stats = FALSE)

### par from Snakemake
if(exists("snakemake")) {
  test <- snakemake@params[["test"]] %>% as.logical
  dir_bed <- snakemake@params[["dir_bed"]]
  dir_top <- snakemake@params[["dir_top"]]
  
  trait <- snakemake@params[["trait"]]
  chr <- snakemake@params[["chr"]]
  num_sel <- snakemake@params[["num_sel"]] %>% as.integer
  
  output <- snakemake@output[["file"]]
} else {
  dir_bed <- "out/bed_top/"
  dir_top <- "out/top/"
  
  test <- TRUE
  trait <- "body_BMIz_orth"
  chr <- 22
  num_sel <- 500
  
  output <- "tmp.tsv"
}

cat(" - test: ", test, "\n")
cat(" - dir_bed: ", dir_bed, "\n")
cat(" - dir_top: ", dir_top, "\n")
cat(" - num_sel: ", num_sel, "\n")
cat(" - trait: ", trait, "\n")
cat(" - chr: ", chr, "\n")
cat(" - output: ", output, "\n")

### top
stopifnot(file.exists(dir_top))
file_top <- list.files(dir_top, full = TRUE) %>% grep(trait, ., value = TRUE)
print(file_top)
stopifnot(length(file_top) == 1)

top <- read_tsv(file_top) %>% 
  head(num_sel) # select top snps to be selected, e.g. 500 out of 5,000

snps_all <- top$SNP

snps_test <- filter(top, CHR == chr) %$% SNP

# early stop
if(length(snps_test) == 0) {
  write_tsv(tibble(), output)
}

snps_grm <- filter(top, CHR != chr) %$% SNP
stopifnot(length(snps_grm) > 0)

### bed files
stopifnot(file.exists(dir_bed))
file_bed <- list.files(dir_bed, "\\.bed$", full = TRUE) %>% grep(trait, ., value = TRUE)
print(file_bed)
stopifnot(length(file_bed) == 1)

y <- ukbjass_val_trait(trait, impute_mean = TRUE)
  
bed <- read.bed.matrix(file_bed)
ids_bed <- paste(bed@ped$famid, bed@ped$id, sep = "_")
stopifnot(all(ids_bed == names(y)))

if(test) {
  nrow1 <- 1e3
  y <- head(y, nrow1)
  Zg <- as.matrix(bed[seq(nrow1), seq(num_sel)])
} else {
  Zg <- as.matrix(bed[, seq(num_sel)])
}

stopifnot(all(snps_all %in% colnames(Zg)))

# impute missing entries
Zg <- impute_mean(Zg)

Z <- scale_z(Zg) # need for assoc

# (2) asoc by lmm
N <- nrow(Z)
M <- ncol(Z)
S <- length(snps_grm)

Zgrm <- Z[, snps_grm] / sqrt(S)

if(test) {
  snps_test <- head(snps_test, 5)
}


assoc_lmm <- lapply(seq_along(snps_test), function(j) {
  cat(" - snp", j, "/", length(snps_test), "\n")
  snp <- snps_test[j]
 
  # snp j is excluded from GR
  X0 <- matrix(1, nrow = nrow(Z), ncol = 1)
  X <- cbind(X0, Z[, snp])

  mod <- biglmm(y, X, Zgrm, scale = FALSE, verbose = 2)
  
  coef <- tail(mod$coef, 1) # the last raw for SNP
  
  tibble(trait = trait, 
      snp = snp, beta = coef$beta, se = coef$se) %>%
    mutate(
      zscore = beta / se,
      pval = pchisq(zscore*zscore, df = 1, lower = FALSE),
      h2 = mod$gamma, N = N, M = M, S = S) 
}) %>% bind_rows

### write

print(assoc_lmm)

write_tsv(assoc_lmm, output)


  
