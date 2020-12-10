TRAITS = ["bmi", "weight", "waist", "hip", "height", "whr"]
# TRAITS = "height"
NTOP = [500, 1000]#, 2000]
# NTOP = 1000
CHR = list(range(1, 23))
# CHR = [1, 2, 3, 20, 21, 22]

#-----------------------
# figures
#-----------------------

rule fig:
    input: ["figures/h2.png", "figures/clump.png"]

rule fig_h2:
    input: "scripts/fig/01-fig-h2.R"
    output: "figures/h2.png"
    shell: "Rscript {input} && cp tmp.png {output}"

rule fig_clump:
    input: "scripts/fig/03-clump-pvals.R"
    output: "figures/clump.png"
    shell: "Rscript {input} && cp tmp.png {output}"

#-----------------------
# GWAS LMM (LOCO)
#-----------------------
rule loco:
    input: expand("out/lmm_loco_top/{ntop}/{trait}.{chr}.tsv.gz", trait = TRAITS, ntop = NTOP, chr = CHR)

rule loco_trait:
    input: "out/ztz/{trait}.rds", "out/h2-loco/{ntop}/{trait}.{chr}.tsv.gz"
    output: "out/lmm_loco_top/{ntop}/{trait}.{chr}.tsv.gz"
    threads: 1000
    shell: "Rscript scripts/07-gwas-lmm-loco-top.R {wildcards.trait} {wildcards.ntop} {wildcards.chr} {output}"

#-----------------------
# GWAS LMM
#-----------------------
rule lmm:
    input: expand("out/lmm_top/{ntop}/{trait}.tsv.gz", trait = TRAITS, ntop = NTOP)

rule lmm_trait:
    input: "out/ztz/{trait}.rds"
    output: "out/lmm_top/{ntop}/{trait}.tsv.gz"
    shell: "Rscript scripts/05-gwas-lmm-top.R {wildcards.trait} {wildcards.ntop} {output}"

#-----------------------
# h2 top (LOCO)
#-----------------------
rule h2loco:
    input: expand("out/h2-loco/{ntop}/{trait}.{chr}.tsv.gz", trait = TRAITS, ntop = NTOP, chr = CHR)

rule h2loco_trait:
    input: "out/ztz/{trait}.rds"
    output: "out/h2-loco/{ntop}/{trait}.{chr}.tsv.gz"
    shell: "Rscript scripts/06-estimate-h2-loco.R {wildcards.trait} {wildcards.ntop} {wildcards.chr} {output}"

#-----------------------
# h2 top 
#-----------------------
rule h2:
    input: expand("out/h2/{ntop}/{trait}.tsv.gz", trait = TRAITS, ntop = NTOP)

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

#-----------------------
# GWAS LR
#-----------------------
rule lm:
    input: expand("out/lm_top/{trait}.tsv.gz", trait = TRAITS)

rule lm_trait:
    output: "out/lm_top/{trait}.tsv.gz"
    threads: 1000
    shell: "Rscript scripts/08-gwas-lm-top.R {wildcards.trait} {output}"

#-----------------------
# Extra analyses
#-----------------------
rule loco_pcs:
    input: expand("out/lmm_loco_pcs_top/{ntop}/{trait}.{chr}.tsv.gz", trait = TRAITS, ntop = NTOP, chr = CHR)

rule loco_pcs_trait:
    input: "out/ztz/{trait}.rds", "out/h2-loco_pcs/{ntop}/{trait}.{chr}.tsv.gz"
    output: "out/lmm_loco_pcs_top/{ntop}/{trait}.{chr}.tsv.gz"
    threads: 1000
    shell: "Rscript scripts/extra/03-gwas-lmm-loco-pcs-top.R {wildcards.trait} {wildcards.ntop} {wildcards.chr} {output}"

rule h2loco_pcs:
    input: expand("out/h2-loco_pcs/{ntop}/{trait}.{chr}.tsv.gz", trait = TRAITS, ntop = NTOP, chr = CHR)

rule h2loco_pcs_trait:
    input: "out/ztz/{trait}.rds"
    output: "out/h2-loco_pcs/{ntop}/{trait}.{chr}.tsv.gz"
    shell: "Rscript scripts/extra/05-estimate-h2-loco-pcs.R {wildcards.trait} {wildcards.ntop} {wildcards.chr} {output}"

rule h2pcs:
    input: expand("out/h2-pcs/{ntop}/{trait}.tsv.gz", trait = TRAITS, ntop = NTOP)

rule h2pcs_trait:
    input: "out/ztz/{trait}.rds"
    output: "out/h2-pcs/{ntop}/{trait}.tsv.gz"
    shell: "Rscript scripts/extra/01-estimate-h2-pcs.R {wildcards.trait} {wildcards.ntop} {output}"

rule ztz5k:
    input: expand("out/ztz5k/{trait}.rds", trait = TRAITS)

rule ztz5k_trait:
    input: "out/fbm/{trait}.rds"
    output: "out/ztz5k/{trait}.rds"
    shell: "Rscript scripts/extra/02-precompute-ztz.R {wildcards.trait} {output}"

rule lm_pcs:
    input: expand("out/lm_pcs_top/{trait}.tsv.gz", trait = TRAITS)

rule lm_pcs_trait:
    output: "out/lm_pcs_top/{trait}.tsv.gz"
    threads: 1000
    shell: "Rscript scripts/extra/04-gwas-lm-pcs-top.R {wildcards.trait} {output}"
