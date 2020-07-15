# paper-neff

Scripts for the manuscript, 
Ziyatdinov, et al. "Estimating the effective sample size in association studies of quantitative traits." bioRxiv (2019),
https://www.biorxiv.org/content/10.1101/2019.12.15.877217v2.full.

## Analysis steps

0. UKB data preparation

- 336K unrelated White British individuals
- >600K genotypes SNPs, autosomes, MAF >0.1%, QC
- 6 traits: bmi, weight, waist, hip, height, whr
  - mean-impute missing trait values (<1% of missingness)
  - project out covariates: age/sex + PC1-20
  - apply rank-based inverse normal transformation

Scripts to process raw UKB data are not shared.

Output:

```
output/
├── ids2.txt
├── ids.txt
├── phen.tsv.gz
├── phen.sync.tsv.gz
├── ukb.bed
├── ukb.bim
├── ukb.fam
└── ukb.log
```

1. GWAS by Linear Regression (LR)

GWAS by LR:

- script [01-gwas-lr/scripts/01-gwas-lr.R](01-gwas-lr/scripts/01-gwas-lr.R)

```
snakemake -s snakemake.py gwas_lr 
```

Clumping by 
[bigsnpr](https://privefl.github.io/bigsnpr/reference/snp_clumping.html):

- script: [01-gwas-lr/scripts/02-clump.R](01-gwas-lr/scripts/02-clump.R)


