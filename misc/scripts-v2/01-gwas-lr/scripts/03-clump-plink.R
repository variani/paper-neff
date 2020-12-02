library(tidyverse)
library(glue)

## par
# https://www.cog-genomics.org/plink/1.9/postproc#clump
clump_p1 <- 0.001 # plink default: 0.0001
clump_p2 <- 0.1 # plink default: 0.01
clump_r2 <- 0.1 # plink default: 0.5; Young 2018: 0.1 for unrelated UKB; 0.01 + MAF>0.4 for related UKB
clump_kb <- 250 # plink default: 250

file_bed <- "output/gen.bed"
bfile_bed <- "output/gen"

## args
args <- commandArgs(trailingOnly = TRUE)
trait <- ifelse(is.na(args[1]), "whr", args[1])
basename_out <- ifelse(is.na(args[2]), "tmp", args[2])

## parallel
cores <- parallel::detectCores()
cores <- ifelse(cores > 1, cores - 1, cores)

## read GWAS LR summary statistics
snps <- read_lines("output/gen.variants")
assoc <- glue("output/gwas-lr/{trait}.tsv.gz") %>% read_tsv
#  predictor           beta      se zscore  pval
#  <chr>              <dbl>   <dbl>  <dbl> <dbl>
#1 1:818802:A:G_A -0.00472  0.00365 -1.29  0.197
#2 1:833068:G:A_A -0.000885 0.00397 -0.223 0.824

assoc_clump <- mutate(assoc,
  predictor = str_remove(predictor, "_.$")) %>%
  select(predictor, pval) %>% rename(SNP = predictor, P = pval)
stopifnot(all(assoc_clump$SNP %in% snps))

## clump
tmp <- tempfile()
# input file for clumping
write.table(assoc_clump, file = paste0(tmp, ".assoc"), row.names = FALSE, quote = FALSE)     

# run plink for clumping
cmd <- glue("plink --bfile {bfile_bed} --out {tmp}",
  " --threads {cores} ",
  " --clump {tmp}.assoc",
  " --clump-p1 {clump_p1} --clump-p2 {clump_p2} --clump-r2 {clump_r2}")
ret <- system(cmd)

# read output
clumped <- read.table(glue("{tmp}.clumped"), header = TRUE) %>% 
  as_tibble %>% select(CHR, SNP, BP, P) %>%
  select(SNP, P) %>% rename(snp = SNP, pval = P)
stopifnot(!duplicated(clumped$snp))

# sort clumped snps by p-value
clumped <- arrange(clumped, pval)

## filter snps with p-value > 0.1
# tab <- filter(tab, pval < 0.1)

## save
write_lines(clumped$snp, glue("{basename_out}.txt.gz"))
write_tsv(clumped, glue("{basename_out}.tsv.gz"))
