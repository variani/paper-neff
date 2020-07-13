# paper-neff
Scripts for the Neff paper 

## Analysis steps

0. UKB data preparation

- 336K unrelated White British individuals
- HapMap3 imputed SNPs, autosomes
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


