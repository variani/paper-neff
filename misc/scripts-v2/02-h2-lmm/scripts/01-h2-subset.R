### inc
library(tidyverse)

library(BEDMatrix)

library(devtools)
load_all("~/git/variani/biglmmz/")

## parameters
trait <- "height"
file_bed <- "output/gen.bed"
file_phen <- "output/phen.bed.tsv.gz"

## data
bed <- BEDMatrix(file_bed)
phen <- read_tsv(file_phen, col_types = c("ccdddddd"))

stopifnot(trait %in% colnames(phen))
y <- phen[[trait]] 

## subset 
n <- 500
p <- 100

y <- y[1:n]
Zg <- bed[1:n, 1:p]

## impute missing entries
Zg <- impute_mean(Zg)

cols_sd <- which(apply(Zg, 2, sd) > 0)
Zg <- Zg[, cols_sd]

## Zgrm -> scaled Z / sqrt(M)
Zgrm <- scale_zg(Zg) 

## lmm
mod <- biglmmz(y, Z = Zgrm, scale = FALSE, verbose = 2)

## evd
K <- crossprod(Zgrm)
lamdas <- eigen(K)$values

N <- length(y)
M <- ncol(Zgrm)
h2 <- mod$gamma

trace_factor <- (sum(1/(h2*lamdas + (1-h2))) + (N-M)/(1-h2)) / N

## output table
tab <- tibble(trait = trait, N = N, M = M, h2 = h2, trace_factor = trace_factor)


  
