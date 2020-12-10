# Analysis steps
# - GWAS by Linear Regression (LR)
# - Clumping variants using p-values from GWAS-LR

# TRAITS = ["bmi", "weight", "waist", "hip", "height", "whr", "pdw", "rbc", "rbcdw"]
# TRAITS = ["pdw", "rbc", "rbcdw"]
TRAITS = "pdw"
# NTOP_H2 = [100, 500, 1000, 2000, 3000, 4000]
# NTOP_H2 = [500, 2000]
NTOP_H2 = [2000]
# CHR = list(range(1, 23))
# CHR = list(range(1, 4))
CHR = 1

rule output:
    output: "output/"
    shell: "mkdir -p {output}"

#-----------------------
# FBM top (Top LR SNPs)
#-----------------------

rule k:
    input: expand("output/k/{trait}.{ntop}.k.rds", trait = TRAITS, ntop = NTOP_H2)

rule k_trait:
    input: "output/top/{trait}.{ntop}.bk"
    output: "output/k/{trait}.{ntop}.k.rds"
    shell: "Rscript 01-gwas-lr/scripts/06-K-top.R {wildcards.trait} {wildcards.ntop} {output}"

rule kj:
    input: expand("output/kj/{trait}.{ntop}.kj.rds", trait = TRAITS, ntop = NTOP_H2)

rule kj_trait:
    input: "output/top/{trait}.{ntop}.bk"
    output: "output/kj/{trait}.{ntop}.kj.rds"
    shell: "Rscript 01-gwas-lr/scripts/05-Kj-top.R {wildcards.trait} {wildcards.ntop} {output}"

rule fbm:
    input: expand("output/top/{trait}.{ntop}.{ext}", trait = TRAITS, ntop = NTOP_H2, ext = ["bk", "rds"])

rule fbm_trait:
    input: "output/h2/{trait}.{ntop}.tsv.gz"
    output: bk = "output/top/{trait}.{ntop}.bk", rds = "output/top/{trait}.{ntop}.rds"
    shell: "Rscript 01-gwas-lr/scripts/04-bed-top.R {wildcards.trait} {wildcards.ntop} {output.bk}"

#-----------------------
# GWAS LMM (Top LR SNPs)
#-----------------------

rule gwas_lmm_top_lr_pcs:
    input: expand("output/gwas-lmm-top-lr-pcs/{ntop}/chr/{trait}.{ntop}.{chr}.tsv.gz", trait = TRAITS, ntop = NTOP_H2, chr = CHR)

rule gwas_lmm_top_lr_pcs_trait:
    threads: 1000
    input: "output/h2_loco/{trait}.{ntop}.{chr}.tsv.gz"
    output: "output/gwas-lmm-top-lr-pcs/{ntop}/chr/{trait}.{ntop}.{chr}.tsv.gz"
    shell: "Rscript 03-gwas-lmm/scripts/10-gwas-lmm-scale-chr-top-lr-pcs.R {wildcards.trait} {wildcards.ntop} {wildcards.chr} {output}"

rule gwas_lmm_top_lr:
    input: expand("output/gwas-lmm-top-lr/{ntop}/chr/{trait}.{ntop}.{chr}.tsv.gz", trait = TRAITS, ntop = NTOP_H2, chr = CHR)

rule gwas_lmm_top_lr_trait:
    # input: "output/h2/{trait}.{ntop}.tsv.gz", "output/gwas-lr/{trait}.tsv.gz"
    output: "output/gwas-lmm-top-lr/{ntop}/chr/{trait}.{ntop}.{chr}.tsv.gz"
    shell: "Rscript 03-gwas-lmm/scripts/09-gwas-lmm-scale-chr-top-lr.R {wildcards.trait} {wildcards.ntop} {wildcards.chr} {output}"

#-----------------------
# GWAS LMM (Top SNPs)
#-----------------------

rule gwas_lmm_top_pcs:
    input: expand("output/gwas-lmm-top-pcs/{ntop}/{trait}.{ntop}.tsv.gz", trait = TRAITS, ntop = NTOP_H2)

rule gwas_lmm_top_pcs_trait:
    input: "output/k/{trait}.{ntop}.k.rds"
    output: "output/gwas-lmm-top-pcs/{ntop}/{trait}.{ntop}.tsv.gz"
    shell: "Rscript 03-gwas-lmm/scripts/12-gwas-lmm-scale-top-pcs.R {wildcards.trait} {wildcards.ntop} {output}"

rule gwas_lmm_top:
    input: expand("output/gwas-lmm-top/{ntop}/{trait}.{ntop}.tsv.gz", trait = TRAITS, ntop = NTOP_H2)

rule gwas_lmm_top_trait:
    input: "output/kj/{trait}.{ntop}.kj.rds"
    output: "output/gwas-lmm-top/{ntop}/{trait}.{ntop}.tsv.gz"
    # shell: "Rscript 03-gwas-lmm/scripts/04-gwas-lmm-top.R {wildcards.trait} {wildcards.ntop} {output}"
    # shell: "Rscript 03-gwas-lmm/scripts/06-gwas-lmm-top-batches.R {wildcards.trait} {wildcards.ntop} {output}"
    shell: "Rscript 03-gwas-lmm/scripts/11-gwas-lmm-scale-top.R {wildcards.trait} {wildcards.ntop} {output}"

#-----------------------
# GWAS LMM (CHR)
#-----------------------

rule gwas_lmm:
    input: expand("output/gwas-lmm/{ntop}/chr/{trait}.{ntop}.{chr}.tsv.gz", trait = TRAITS, ntop = NTOP_H2, chr = CHR)

rule gwas_lmm_trait_chr:
    input: "output/h2_loco/{trait}.{ntop}.{chr}.tsv.gz"
    output: "output/gwas-lmm/{ntop}/chr/{trait}.{ntop}.{chr}.tsv.gz"
    shell: "Rscript 03-gwas-lmm/scripts/08-gwas-chr-scale.R {wildcards.trait} {wildcards.ntop} {wildcards.chr} {output}"

#-----------------------
# Estimating h2
#-----------------------

rule h2_loco:
    input: expand("output/h2_loco/{trait}.{ntop}.{chr}.tsv.gz", trait = TRAITS, ntop = NTOP_H2, chr = CHR)

rule h2_loco_trait:
    input: "output/clump-plink/{trait}.tsv.gz", "output/top/{trait}.{ntop}.bk",
        "output/k/{trait}.{ntop}.k.rds"
    output: "output/h2_loco/{trait}.{ntop}.{chr}.tsv.gz"
    shell: "Rscript 02-h2-lmm/scripts/05-h2-loco.R {wildcards.trait} {wildcards.ntop} {wildcards.chr} {output}"

rule h2:
    input: expand("output/h2/{trait}.{ntop}.tsv.gz", trait = TRAITS, ntop = NTOP_H2)

rule h2_trait:
    input: "output/clump-plink/{trait}.tsv.gz"
    output: "output/h2/{trait}.{ntop}.tsv.gz"
    params: trait = "{trait}", ntop = "{ntop}"
    shell: "Rscript 02-h2-lmm/scripts/04-h2-top-fbm.R {params.trait} {params.ntop} {output}"


#-----------------------
# Clumping (plink)
#-----------------------

rule clump:
    input: expand("output/clump-plink/{trait}.{ext}", ext = ["tsv.gz", "txt.gz"], trait = TRAITS)

rule clump_trait:
    input: "output/gwas-lr/{trait}.tsv.gz"
    output: "output/clump-plink/{trait}.tsv.gz", "output/clump-plink/{trait}.txt.gz"
    params: trait = "{trait}"
    shell: "Rscript 01-gwas-lr/scripts/03-clump-plink.R {params.trait} output/clump-plink/{params.trait}"

#-----------------------
# Clumping (bigsnpr)
#-----------------------

rule clump_bigsnpr:
    input: expand("output/clump/{trait}.{ext}", ext = ["tsv.gz", "txt.gz"], trait = TRAITS)

rule clump_trait_bigsnpr:
    input: "output/gwas-lr/{trait}.tsv.gz"
    output: "output/clump/{trait}.tsv.gz", "output/clump/{trait}.txt.gz"
    params: trait = "{trait}"
    shell: "Rscript 01-gwas-lr/scripts/02-clump.R {params.trait} output/clump/{params.trait}"

#-----------------------
# GWAS LR
#-----------------------

rule gwas_lm:
    input: expand("output/gwas-lm/{trait}.tsv.gz", trait = TRAITS)

rule gwas_lm_trait:
    output: "output/gwas-lm/{trait}.tsv.gz"
    params: trait = "{trait}"
    shell: "Rscript 01-gwas-lr/scripts/01c-gwas-lm.R {params.trait} {output}"

rule gwas_lr2:
    input: expand("output/gwas-lr2/{trait}.tsv.gz", trait = TRAITS)

rule gwas_lr2_trait:
    output: "output/gwas-lr2/{trait}.tsv.gz"
    params: trait = "{trait}"
    shell: "Rscript 01-gwas-lr/scripts/01b-gwas-lr-matlm.R {params.trait} {output}"

rule gwas_lr_trait:
    output: "output/gwas-lr/{trait}.tsv.gz"
    params: trait = "{trait}"
    shell: "Rscript 01-gwas-lr/scripts/01-gwas-lr.R {params.trait} {output}"

