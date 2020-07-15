### inc
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

  trait <- snakemake@params[["trait"]]
  num_sel <- snakemake@params[["num_sel"]] %>% as.integer
  
  output <- snakemake@output[["file"]]
} else {
  dir_bed <- "out/bed_top/"
  test <- TRUE
  trait <- "body_BMIz_orth"
  num_sel <- 500
  
  output <- "tmp.tsv"
}

cat(" - test: ", test, "\n")
cat(" - dir_bed: ", dir_bed, "\n")
cat(" - num_sel: ", num_sel, "\n")
cat(" - trait: ", trait, "\n")
cat(" - output: ", output, "\n")

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

# impute missing entries
Zg <- impute_mean(Zg)

Zgrm <- scale_zg(Zg) # need for EVD-based calculation of the trace factor

### lmm
#mod <- biglmm(y, Z = Zg, scale = TRUE, verbose = 2)
mod <- biglmm(y, Z = Zgrm, scale = FALSE, verbose = 2)

### evd
K <- crossprod(Zgrm)
lamdas <- eigen(K)$values

N <- length(y)
M <- ncol(Zgrm)
h2 <- mod$gamma

trace_factor <- (sum(1/(h2*lamdas + (1-h2))) + (N-M)/(1-h2)) / N

### write
tab <- tibble(trait = trait, N = N, M = M, h2 = h2, trace_factor = trace_factor)

write_tsv(tab, output)

  
