# Analysis steps
# - GWAS by Linear Regression (LR)
# - Clumping variants using p-values from GWAS-LR

TRAITS = ["bmi", "weight", "waist", "hip", "height", "whr"]
#TRAITS = ["height"]
NTOP_H2 = [100, 500, 1000, 2000]
# NTOP_H2 = [100, 500]

rule output:
    output: "output/"
    shell: "mkdir -p {output}"

#-----------------------
# Estimating h2
#-----------------------

rule h2:
    input: expand("output/h2/{trait}.{ntop}.tsv.gz", trait = TRAITS, ntop = NTOP_H2)

rule h2_trait:
    input: "output/clump/{trait}.tsv.gz"
    output: "output/h2/{trait}.{ntop}.tsv.gz"
    params: trait = "{trait}", ntop = "{ntop}"
    shell: "Rscript 02-h2-lmm/scripts/02-h2-top.R {params.trait} {params.ntop} {output}"

#-----------------------
# Clumping
#-----------------------

rule clump:
    input: expand("output/clump/{trait}.{ext}", ext = ["tsv.gz", "txt.gz"], trait = TRAITS)

rule clump_trait:
    input: "output/gwas-lr/{trait}.tsv.gz"
    output: "output/clump/{trait}.tsv.gz", "output/clump/{trait}.txt.gz"
    params: trait = "{trait}"
    shell: "Rscript 01-gwas-lr/scripts/02-clump.R {params.trait} output/clump/{params.trait}"

#-----------------------
# GWAS LR
#-----------------------

rule gwas_lr_trait:
    output: "output/gwas-lr/{trait}.tsv.gz"
    params: trait = "{trait}"
    shell: "Rscript 01-gwas-lr/scripts/01-gwas-lr.R {params.trait} {output}"

