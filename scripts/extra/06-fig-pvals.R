library(tidyverse)
library(glue)
library(data.table)
library(RColorBrewer) 

library(cowplot)
library(latex2exp)

### settings
theme_set(theme_cowplot(8))

# plotting par.
# trait <- c("bmi", "weight", "waist", "hip", "height", "whr")[2]
trait <- "height"
ntop <- 1000
vals_chr <- 1:22

cat("trait", trait, "\n")

f1 <- glue("out/lm_pcs_top/{trait}.tsv.gz")
t1 <- fread(f1) %>% as_tibble

t1 <- select(t1, snp, beta, se, zscore, pval) %>% 
  dplyr::rename(se_lr = se, b_lr = beta, z_lr = zscore, p_lr = pval)

t2 <- glue("out/lmm_loco_pcs_top/{ntop}/{trait}.{vals_chr}.tsv.gz") %>%
  lapply(fread) %>% bind_rows %>% as_tibble
t2 <- select(t2, snp, beta, se, zscore, pval) %>% 
  dplyr::rename(se_lmm = se, b_lmm = beta, z_lmm = zscore, p_lmm = pval)

t <- left_join(t2, t1, by = "snp")
thr <- 5e-8
ts <- filter(t, p_lr < thr & p_lmm < thr)

## plot
up_p <- function(p) p + geom_point() + geom_abline()

p1 <- ggplot(t, aes(abs(b_lr), abs(b_lmm))) %>% up_p
p1 <- p1 + geom_smooth()

tb <- filter(t, abs(b_lmm) > 0.01 & abs(b_lmm) < 0.05) 
p11 <- ggplot(tb, aes(abs(b_lr), abs(b_lmm))) %>% up_p
p11 <- p11 + geom_smooth(color = "red")

p2 <- ggplot(ts, aes(1/se_lr^2, 1/se_lmm^2)) %>% up_p

p3 <- ggplot(t, aes(-log10(p_lr), -log10(p_lmm))) %>% up_p
p3 <- p3 + geom_smooth()

p33 <- ggplot(ts, aes(-log10(p_lr), -log10(p_lmm))) %>% up_p
p33 <- p33 + geom_smooth(color = "red")

g <- plot_grid(p1, p11, p2, p3)
ggsave("tmp.png", g)
