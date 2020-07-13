# Analysis steps
# - GWAS by Linear Regression (LR)
# - Clumping variants using p-values from GWAS-LR

TRAITS = ["bmi", "weight", "waist", "hip", "height", "whr"]
#TRAITS = ["height", "bmi"]

rule output:
    output: "output/"
    shell: "mkdir -p {output}"

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

