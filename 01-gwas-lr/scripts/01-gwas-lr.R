library(tidyverse)

library(devtools)
load_all("~/git/variani/matlm") # gwas assoc. model
load_all("~/git/variani/bigcov") # split bed by batches

library(BEDMatrix) # read bed from file & avoid loading into RAM

library(doParallel)
library(parallel)
library(plyr)

## args
args <- commandArgs(trailingOnly = TRUE)
trait <- ifelse(is.na(args[1]), "bmi", args[1])
file_out <- ifelse(is.na(args[2]), "tmp.tsv.gz", args[2])

## parallel
cores <- parallel::detectCores()
cores <- ifelse(cores > 1, cores - 1, cores)
parallel <- (cores > 1)
if(parallel) { doParallel::registerDoParallel(cores = cores) }

## par
nsnps_batch <- 250

## data
file_bed <- "output/ukb.bed"
file_phen <- "output/phen.tsv.gz"

## bed
bed <- BEDMatrix(file_bed)
ids_bed <- rownames(bed)

# fix id names such "0_ID_ID"
ids_bed <- ids_bed %>% strsplit("_") %>% sapply(function(x) tail(x, 1))

## phen: load phen, align ids with bed, extract y
phen <- read_tsv(file_phen, col_types = c("ccdddddd"))
ids_phen <- phen[["IID"]]

ind_phen <- match(ids_phen, ids_bed) 
stopifnot(any(!is.na(ind_phen)))

stopifnot(trait %in% colnames(phen))
y <- phen[[trait]] %>% .[ind_phen] # now y is aligned with bed by ids

## assoc
bdat <- bigdat(bed, batch_size = nsnps_batch)
B <- bigdat_nbatch(bdat)

assoc <- llply(seq(B), function(b) {
  cat(" -", b, "/", B, "\n")
  X <- bigdat_batch(bdat, b)
  matlm::matreg0(y, X, verbose = 2)
}, .parallel = parallel) %>% bind_rows

## save
write_tsv(assoc, file_out)
