TRAITS = ["bmi", "weight", "waist", "hip", "height", "whr"]
# TRAITS = ["bmi", "height", "hip"]
NTOP_H2 = [500, 1000, 2000]

#-----------------------
# figures
#-----------------------

rule fig:
    input: "figures/h2.png"

rule fig_h2:
    input: "scripts/fig/01-fig-h2.R"
    output: "figures/h2.png"
    shell: "Rscript {input} && cp tmp.png {output}"

#-----------------------
# h2 top 
#-----------------------
rule h2:
    input: expand("out/h2/{ntop}/{trait}.tsv.gz", trait = TRAITS, ntop = NTOP_H2)

rule h2_trait:
    input: "out/ztz/{trait}.rds"
    output: "out/h2/{ntop}/{trait}.tsv.gz"
    shell: "Rscript scripts/04-estimate-h2.R {wildcards.trait} {wildcards.ntop} {output}"

#-----------------------
# ZtZ top 
#-----------------------
rule ztz:
    input: expand("out/ztz/{trait}.rds", trait = TRAITS)

rule ztz_trait:
    input: "out/fbm/{trait}.rds"
    output: "out/ztz/{trait}.rds"
    shell: "Rscript scripts/03-precompute-ztz.R {wildcards.trait} {output}"

#-----------------------
# FBM top 
#-----------------------
rule fbm:
    input: expand("out/fbm/{trait}.rds", trait = TRAITS)

rule fbm_trait:
    output: "out/fbm/{trait}.rds"
    shell: "Rscript scripts/02-convert-bed2fbm.R {wildcards.trait} {output}"

