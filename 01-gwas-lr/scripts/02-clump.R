library(tidyverse)
library(glue)

library(bigsnpr)

## par
r2 <- 0.01 # plink default: 0.01
size <- 250 # plink default: 250 (distance in kb)

file_bed <- "output/ukb.bed"

## args
args <- commandArgs(trailingOnly = TRUE)
trait <- ifelse(is.na(args[1]), "weight", args[1])
basename_out <- ifelse(is.na(args[2]), "tmp", args[2])

## parallel
cores <- parallel::detectCores()
cores <- ifelse(cores > 1, cores - 1, cores)

## read GWAS LR summary statistics
assoc <- glue("output/gwas-lr/{trait}.tsv.gz") %>% read_tsv
#  predictor           beta      se zscore  pval
#  <chr>              <dbl>   <dbl>  <dbl> <dbl>
#1 1:818802:A:G_A -0.00472  0.00365 -1.29  0.197
#2 1:833068:G:A_A -0.000885 0.00397 -0.223 0.824

# ceheck sync. of snps names

snps_bed <- read_lines("output/gen.variants")
#stopifnot(all(assoc[[1]] == snps_bed))

scores <- tibble(snp = snps_bed) %>%
  left_join(tibble(snp = assoc[[1]], pval = assoc[["pval"]],
    nlogp = -log10(assoc[["pval"]])))

# fill missing scores for snps not present in `assoc`
scores <- mutate(scores, nlogp = ifelse(is.na(nlogp), 0, nlogp))

## clump
bed <- bed(file_bed)
ind_keep <- bed_clumping(bed, thr.r2 = r2, size = size,
  S = scores$nlogp,
  ncores = cores)

snps_keep <- snps_bed[ind_keep]

## filter snps with p-value > 0.1
tab <- filter(scores, snp %in% snps_keep) %>% select(snp, pval)
stopifnot(nrow(tab) == length(snps_keep))

tab <- filter(tab, pval < 0.1)
snps_keep <- tab$snp

## save
write_lines(snps_keep, glue("{basename_out}.txt.gz"))

write_tsv(tab, glue("{basename_out}.tsv.gz"))
