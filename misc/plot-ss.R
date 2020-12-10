library(tidyverse)
library(glue)

library(cowplot)
theme_set(theme_cowplot(14))

vals_traits <- c("body_HEIGHTz_orth", "body_BMIz_orth",
  "body_HIP_orth", "body_WAIST_orth", "body_WEIGHT_orth", "body_WHR_orth")  
trait <- vals_traits[1]
vals_chr <- 1:22
ntop <- 1000

size <- 2
alpha <- 0.75
thr <- 5e-8

gamma <- glue("misc/out_o2/h2/{trait}.{ntop}.tsv.gz") %>% read_tsv %>% .[["trace_factor"]]

t1 <- glue("misc/out_o2/assoc_lm/{trait}.{ntop}.tsv.gz") %>% read_tsv

t2 <- lapply(vals_chr, function(chr) {
  f <- glue("misc/out_o2/assoc_lmm_chr/{chr}.{trait}.{ntop}.tsv.gz") 
  `if`(file.exists(f), read_tsv(f), tibble())
}) %>% bind_rows

t <- left_join(
  dplyr::rename(t2, b_lmm = beta, se_lmm = se, p_lmm = pval) %>% 
    select(snp, b_lmm, se_lmm, p_lmm),
  dplyr::rename(t1, b_lr = beta, se_lr = se, p_lr = pval) %>% 
    select(snp, b_lr, p_lr, se_lr)) 

p1 <- ggplot(t, aes(abs(b_lr), abs(b_lmm))) + 
  geom_point(size = size, alpha = alpha) +
  geom_abline(linetype = 3)
  # facet_wrap(~ chr, scales = "free")

p2 <- ggplot(filter(t, p_lr < thr), aes((1/se_lr^2), (1/se_lmm^2))) + 
  geom_point(size = size, alpha = alpha) +
  geom_abline(linetype = 3)
  # facet_wrap(~ chr, scales = "free")

p3 <- ggplot(t, aes(-log10(p_lr), -log10(p_lmm))) +
  geom_point(size = size, alpha = alpha) +
  geom_abline(linetype = 3) +
  geom_hline(yintercept = -log10(thr), linetype = 3) +
  geom_vline(xintercept = -log10(thr), linetype = 3) 
  # facet_wrap(~ chr, scales = "free")

g <- plot_grid(p1, p2, p3, labels = "auto")
ggsave("tmp.png", plot = g, dpi = 100)
