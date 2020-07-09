# Analysis steps
# - GWAS by Linear Regression (LR)
# - Clumping variants using p-values from GWAS-LR

TRAITS = ["bmi", "weight", "waist", "hip", "height", "whr"]

rule output:
    output: "output/"
    shell: "mkdir -p {output}"

#-----------------------
# GWAS LR
#-----------------------

rule gwas_lr:
    input: expand("output/gwas-lr/{trait}.tsv.gz", trait = TRAITS)

rule gwas_lr_trait:
    output: "output/gwas-lr/{trait}.tsv.gz"
    params: trait = "{trait}"
    shell: "Rscript 01-gwas-lr/scripts/01-gwas-lr.R {params.trait} {output}"

