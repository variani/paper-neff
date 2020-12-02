library(tidyverse)
library(glue)

library(cowplot)
theme_set(theme_cowplot(14))

# vals_traits <- c("bmi", "height")
# vals_chr <- 21:22 
# vals_traits <- c("height", "pdw")
vals_traits <- c("bmi", "weight", "waist", "hip", "height", "whr", "pdw", "rbc", "rbcdw")
# vals_chr <- 1:3
ntop <- 4000

chr <- 1

vals_est <- c("mean", "median")
vals_filt <- paste0("f", 0:3)

thr2 <- c(1e-3, 1e-5, 5e-8)[3]
thr1_lmm <- thr2
thr1_lr <- 0.05
thr3 <- thr2

tab <- lapply(vals_traits, function(trait) {
  cat("trait", trait, "\n")

  # h2 <- glue("output/h2/{trait}.{ntop}.tsv.gz") %>% read_tsv 
  h2 <- glue("output/h2_loco/{trait}.{ntop}.{chr}.tsv.gz") %>% read_tsv 
  gamma0 <- h2$trace_factor
  gamma <- with(h2, trace_factor / s2)
  cat(" - gamma", gamma, "; gamma (not scaled by s2) ", gamma0, "\n")

  t1 <- glue("output/gwas-lr/{trait}.tsv.gz") %>% read_tsv
  t1 <- select(t1, predictor, se, zscore, pval) %>% 
    dplyr::rename(se_lr = se, z_lr = zscore, p_lr = pval)

  # f2 <- glue("output/gwas-lmm-1000/chr/{trait}.{ntop}.{vals_chr}.tsv.gz")
  # f2 <- glue("output/gwas-lmm-top-lr-1000/chr/{trait}.{ntop}.{vals_chr}.tsv.gz")
  # t2 <- lapply(f2, read_tsv) %>% bind_rows
  t2 <- glue("output/gwas-lmm-top-lr-pcs/{ntop}/chr/{trait}.{ntop}.{1:22}.tsv.gz") %>% 
    lapply(function(f) `if`(file.exists(f), read_tsv(f), tibble())) %>% bind_rows
  # t2 <- glue("output/gwas-lmm/{ntop}/chr/height.{ntop}.3.tsv.gz") %>% read_tsv
  t2 <- select(t2, predictor, se, zscore, pval) %>% 
    dplyr::rename(se_lmm = se, z_lmm = zscore, p_lmm = pval)

  t <- left_join(t2, t1)

  lapply(vals_filt, function(filt) {
    cat("filter", filt, "\n")

    t <- switch(filt, 
      "f0" = t,
      "f1" = filter(t, p_lr < thr1_lr & p_lmm < thr1_lmm),
      "f2" = filter(t, p_lr < thr2 & p_lmm < thr2),
      "f3" = filter(t, p_lr < thr3 & p_lmm < thr3),
      stop("filt"))

    print(t)

    vals_se2 <- with(t, (se_lr / se_lmm)^2)
    vals_z2 <- with(t, (z_lmm / z_lr)^2)
    # vals_z2 <- with(t, (z_lr / z_lmm)^2)

    lapply(vals_est, function(est) {
      se2 <- switch(est, "median" = median(vals_se2), "mean" = mean(vals_se2), stop("est"))
      z2 <- switch(est, "median" = median(vals_z2), "mean" = mean(vals_z2), stop("est"))

      tibble(
          gamma = c(gamma, se2, z2), 
          estimator = c("trace", "se2", "z2"),
          q25 = c(NA, quantile(vals_se2, 0.25), quantile(vals_z2, 0.25)),
          q75 = c(NA, quantile(vals_se2, 0.75), quantile(vals_z2, 0.75))) %>%
        mutate(trait = trait, filter = filt, m = nrow(t), est = est)
      }) %>% bind_rows
  }) %>% bind_rows
}) %>% bind_rows

## plot
# NB: filter se2 out
# ptab <- filter(tab, estimator != "se2")
ptab <- tab
ptab <- filter(ptab, !(estimator == "z2" & filter == "f0"))

ptab <- filter(ptab, filter != "f0" & est == "median")

offset <- 0.9
# ylims <- c(0, 0.5)
p <- ggplot(ptab, aes(trait, gamma - offset, fill = estimator, group = estimator)) + 
  geom_bar(stat = "identity", position = position_dodge(0.9)) + 
  geom_errorbar(aes(ymin = q25 - offset, ymax = q75 - offset, group = estimator), 
    width = 0.3, position = position_dodge(0.9)) +
  geom_hline(yintercept = 1 - offset, linetype = 3, color = "grey")

p <- p + facet_grid(est ~ filter)

p <- p + 
  scale_y_continuous(labels = function(x) x + offset) +
  theme(legend.position = "top") +
  labs(x = NULL, y = NULL)

ggsave("tmp.png", plot = p, dpi = 100, width = 12, height = 6)

