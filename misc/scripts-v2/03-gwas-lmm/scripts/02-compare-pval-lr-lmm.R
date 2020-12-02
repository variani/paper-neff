library(tidyverse)
library(glue)

library(cowplot)
theme_set(theme_cowplot(14))

# trait <- "bmi" # "height"
trait <- "height"
# vals_chr <- 21:22 # c(4, 7, 9:15, 18, 20, 21)
vals_chr <- c(4, 7, 9:15, 18, 20, 21)
ntop <- 1000

thr <- 5e-8
size <- 2
alpha <- 0.75

top <- glue("output/clump/{trait}.txt.gz") %>% 
  read_lines(n_max = n_top)

t1 <- glue("output/gwas-lr/{trait}.tsv.gz") %>% read_tsv
dups <- filter(t1, duplicated(predictor)) %>% .[["predictor"]]
stopifnot(length(dups) == 0)

f2 <- glue("output/gwas-lmm-top-lr-1000/chr/{trait}.{ntop}.{vals_chr}.tsv.gz")
t2 <- lapply(f2, read_tsv) %>% bind_rows
dups <- filter(t2, duplicated(predictor)) %>% .[["predictor"]]
stopifnot(length(dups) == 0)

t <- left_join(
  dplyr::rename(t2, se_lmm = se, p_lmm = pval) %>% select(predictor, se_lmm, p_lmm, chr),
  dplyr::rename(t1, se_lr = se, p_lr = pval) %>% select(predictor, p_lr, se_lr))
t <- mutate(t, top = predictor %in% top)

p <- ggplot(t, aes(-log10(p_lr), -log10(p_lmm))) +
  geom_point(size = size, alpha = alpha) +
  geom_abline(linetype = 3) +
  geom_hline(yintercept = -log10(thr), linetype = 3) +
  geom_vline(xintercept = -log10(thr), linetype = 3) 
  # facet_wrap(~ chr, scales = "free")

ggsave("tmp.png", plot = p, dpi = 100)

