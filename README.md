 paper-neff

Scripts for the manuscript, 
Ziyatdinov, et al. "Estimating the effective sample size in association studies of quantitative traits." bioRxiv (2019),
https://www.biorxiv.org/content/10.1101/2019.12.15.877217v2.full.

## Analysis steps

### Step 0. UKB data preparation

- 336K unrelated White British individuals
- &gt;600K genotypes SNPs, autosomes, MAF >0.1%, QC
- 6 traits: bmi, weight, waist, hip, height, whr
  - mean-impute missing trait values (<1% of missingness)
  - project out covariates: age/sex + PC1-20
  - apply rank-based inverse normal transformation

Scripts to process raw UKB data are not shared.

Output:

```
output/
├── gen.bed
├── gen.bim
├── gen.fam
├── gen.ids
├── gen.log
├── gen.variants
├── phen.bed.tsv.gz
└── phen.tsv.gz
```

### Step 1. GWAS by Linear Regression (GWAS-LR)

GWAS by LR

- script [01-gwas-lr/scripts/01-gwas-lr.R](01-gwas-lr/scripts/01-gwas-lr.R)

```
snakemake -s snakemake.py gwas_lr 
```

Output:

```
output/gwas-lr/
├── bmi.tsv.gz
├── height.tsv.gz
├── hip.tsv.gz
├── waist.tsv.gz
├── weight.tsv.gz
└── whr.tsv.gz
```

### Step 3: Clumping with p-values from GWAS-LR

Clumping by 
[bigsnpr](https://privefl.github.io/bigsnpr/reference/snp_clumping.html)
using the default plink settings (r2 0.01, window 250kb)

- script: [01-gwas-lr/scripts/02-clump.R](01-gwas-lr/scripts/02-clump.R)

Output:

```
output/clump/
├── height.tsv.gz
└── height.txt.gz
```

### Step 4: Heritability estimation by low-rank LMM

- script: [02-h2-lmm/scripts/01-h2.R](02-h2-lmm/scripts/01-h2.R)

### Step 5: GWAS by low-rank LMM (GWAS-LMM)

- script: [03-gwas-lmm/scripts/01-gwas-lmm-chr.R](03-gwas-lmm/scripts/01-gwas-lmm-chr.R)

## Common R code

```r
library(tidyverse)
library(BEDMatrix)

file_bed <- "output/gen.bed"
file_phen <- "output/phen.bed.tsv.gz"

# read genotypes from bed
snps <- read_lines("output/gen.variants")
bed <- BEDMatrix(file_bed, p = length(snps))
colnames(bed) <- snps

# read phenotypes
phen <- read_tsv(file_phen, col_types = c("ccdddddd"))
stopifnot(trait %in% colnames(phen))
y <- phen[[trait]] 
```
