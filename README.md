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

### Step 1. GWAS by Linear Regression (GWAS-LR)

### Step 3: Clumping with p-values from GWAS-LR

### Step 4: Heritability estimation by low-rank LMM

- script: [scripts/04-estimate-h2.R](scripts/04-estimate-h2.R)

![](figures/h2.png)


### Step 5: GWAS by low-rank LMM (GWAS-LMM)

- script: [scripts/05-gwas-lmm-top.R](scripts/05-gwas-lmm-top.R)

