library(tidyverse)
library(glue)

library(cowplot)
theme_set(theme_cowplot(14))

trait <- "pdw"
ntop <- 2000

size <- 2
alpha <- 0.75
thr <- 5e-8

# t0 <- glue("output/gwas-lr/{trait}.tsv.gz") %>% read_tsv
# t0 <- read_tsv(glue("output/gwas-lm/{trait}.tsv.gz")) %>% dplyr::rename(beta = b)
t1 <- glue("output/gwas-lr/{trait}.tsv.gz") %>% read_tsv
t0 <- t1
# t1 <- read_tsv("output/gwas-lr2/height.tsv.gz") %>% dplyr::rename(beta = b)

# t2 <- glue("output/gwas-lmm-top/{ntop}/{trait}.{ntop}.tsv.gz") %>% read_tsv
# t2 <- read_tsv("output/gwas-lmm/2000/chr/height.2000.3.tsv.gz")
# t2 <- glue("output/gwas-lmm-4000/chr/height.4000.{1:3}.tsv.gz") %>% lapply(read_tsv) %>% bind_rows
# t2 <- glue("bck-output/gwas-lmm-top-lr/{ntop}/chr/height.2000.{1:22}.tsv.gz") %>% 
#   lapply(function(f) `if`(file.exists(f), read_tsv(f), tibble())) %>% bind_rows
# t2 <- glue("output/gwas-lmm-top-lr-pcs/{ntop}/chr/{trait}.{ntop}.{1:22}.tsv.gz") %>% 
#   lapply(function(f) `if`(file.exists(f), read_tsv(f), tibble())) %>% bind_rows
# t2 <- glue("output/gwas-lmm-top/{ntop}/{trait}.{ntop}.tsv.gz") %>% read_tsv
t2 <- glue("output/gwas-lmm-top-pcs/{ntop}/{trait}.{ntop}.tsv.gz") %>% read_tsv

t <- left_join(
  # dplyr::rename(t2, se_lmm = se, p_lmm = pval, predictor = snp) %>% select(predictor, se_lmm, p_lmm),
  dplyr::rename(t2, b_lmm = beta, se_lmm = se, p_lmm = pval) %>% 
    select(predictor, b_lmm, se_lmm, p_lmm),
  dplyr::rename(t1, b_lr = beta, se_lr = se, p_lr = pval) %>% 
    select(predictor, b_lr, p_lr, se_lr)) %>%
    left_join(dplyr::rename(t0, b_lr0 = beta, se_lr0 = se, p_lr0 = pval) %>%
      select(predictor, b_lr0, p_lr0, se_lr0))
# t <- mutate(t, top = predictor %in% top)

p0 <- ggplot(t, aes(abs(b_lr0), abs(b_lr))) + 
  geom_point(size = size, alpha = alpha) +
  geom_abline(linetype = 3)

p1 <- ggplot(t, aes(abs(b_lr), abs(b_lmm))) + 
  geom_point(size = size, alpha = alpha) +
  geom_abline(linetype = 3) 
  # facet_wrap(~ chr, scales = "free")

b_max <- with(t, max(abs(c(b_lmm, b_lr))))
lims <- c(0.02, b_max)
p1b <- p1 + xlim(lims) + ylim(lims) + geom_smooth(method = "lm", color = "red") 
p1 <- p1 + geom_smooth(method = "lm")

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

nlogp_max <- with(t, max(-log10(c(p_lmm, p_lr))))
lims <- c(7, 50)
p3b <- p3 + xlim(lims) + ylim(lims) + geom_smooth(method = "lm", color = "red")
p3 <- p3 + geom_smooth(method = "lm")
  
# g <- plot_grid(p0, p1, p2, p3, labels = "auto")
g <- plot_grid(p1, p1b, p3, p2, labels = "auto")
ggsave("tmp.png", plot = g, dpi = 100)

