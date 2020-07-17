### inc
library(tidyverse)
library(glue)

library(BEDMatrix)

library(devtools)
load_all("~/git/variani/biglmmz/")

## arguments
args <- commandArgs(trailingOnly = TRUE)
trait <- ifelse(is.na(args[1]), "whr", args[1])
n_top <- ifelse(is.na(args[2]), 500, as.integer(args[2]))
file_out <- ifelse(is.na(args[3]), "tmp.tsv.gz", args[3])

## parameters
file_bed <- "output/gen.bed"
file_phen <- "output/phen.bed.tsv.gz"
file_clump <- glue("output/clump/{trait}.tsv.gz")

## data
# read genotypes from bed
snps <- read_lines("output/gen.variants")
bed <- BEDMatrix(file_bed, p = length(snps))
colnames(bed) <- snps

# read phenotypes
phen <- read_tsv(file_phen, col_types = c("ccdddddd"))
stopifnot(trait %in% colnames(phen))
y <- phen[[trait]] 

# read top variants
clump <- read_tsv(file_clump) %>% arrange(pval) 
if(nrow(clump) < n_top) {
  stop(glue("#rows in clump ({nrow(clump)}) < #top ({n_top})"))
}

clump <- head(clump, n_top)
snps_top <- clump$snp

## impute missing entries
Zg <- bed[, snps_top]

# Zg <- impute_mean(Zg)
# stopifnot(!any(apply(Zg, 2, sd) == 0))
# Zgrm <- scale_zg(Zg) # Zgrm -> scaled Z / sqrt(M)

## lmm
mod <- biglmmz(y, Z = Zg, impute = TRUE, scale = TRUE, verbose = 2)
# - impute: impute missing values in Zg by means
# - scale: scale genotypes in Zg to have Zgrm = scaled(Z) / sqrt(M)

# ## evd
# K <- crossprod(Zgrm)
# lamdas <- eigen(K)$values

# N <- length(y)
# M <- ncol(Zgrm)
# h2 <- mod$gamma

# trace_factor <- (sum(1/(h2*lamdas + (1-h2))) + (N-M)/(1-h2)) / N

# ## output table
# tab <- tibble(trait = trait, N = N, M = M, s2 = mod$s2, h2 = h2,
#   trace_factor = trace_factor)

tab <- mod$ess

## save
write_tsv(tab, file_out)
  
