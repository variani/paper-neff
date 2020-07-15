library(tidyverse)

library(devtools)
load_all("~/git/variani/matlm") # gwas assoc. model
load_all("~/git/variani/bigcov") # split bed by batches
load_all("~/git/variani/biglmmz") # impute/scal large matrices 

library(BEDMatrix) # read bed from file & avoid loading into RAM

library(doParallel)
library(parallel)
library(plyr)

## args
args <- commandArgs(trailingOnly = TRUE)
trait <- ifelse(is.na(args[1]), "height", args[1])
file_out <- ifelse(is.na(args[2]), "tmp.tsv.gz", args[2])

## parallel
cores <- parallel::detectCores()
cores <- ifelse(cores > 1, cores - 1, cores)
parallel <- (cores > 1)
if(parallel) { doParallel::registerDoParallel(cores = cores) }

## par
nsnps_batch <- 150

## data
file_bed <- "output/gen.bed"
file_phen <- "output/phen.bed.tsv.gz"

## bed
bed <- BEDMatrix(file_bed)

## phen: load phen sync. with bed ids
phen <- read_tsv(file_phen, col_types = c("ccdddddd"))

stopifnot(trait %in% colnames(phen))
y <- phen[[trait]] 

## assoc
bdat <- bigdat(bed, batch_size = nsnps_batch)
B <- bigdat_nbatch(bdat)

assoc <- llply(seq(B), function(b) {
  cat(" -", b, "/", B, "\n")
  X <- bigdat_batch(bdat, b)

  Z <- impute_mean(X)
  Z <- scale_z(Z)

  assoc <- matlm::matreg0(y, Z, verbose = 2)

  n_hits <- filter(assoc, pval < 5e-8) %>% nrow
  if(n_hits > 0) {
    cat(" -- #hits:", n_hits, "/", nrow(assoc), "\n")
  }

  assoc
}, .parallel = parallel) %>% bind_rows

## save
write_tsv(assoc, file_out)
