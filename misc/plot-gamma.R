library(tidyverse)
library(glue)

library(cowplot)
theme_set(theme_cowplot(14))

vals_traits <- c("body_HEIGHTz_orth", "body_BMIz_orth",
  "body_HIP_orth", "body_WAIST_orth", "body_WEIGHT_orth", "body_WHR_orth")  
vals_chr <- 1:22
ntop <- 1000

vals_est <- c("mean", "median")
vals_filt <- paste0("f", 0:3)

thr2 <- c(1e-3, 1e-5, 5e-8)[3]
thr1_lmm <- thr2
thr1_lr <- 0.05
thr3 <- thr2

tab <- lapply(vals_traits, function(trait) {
  cat("trait", trait, "\n")

  # gamma <- glue("output/h2/{trait}.{ntop}.tsv.gz") %>% read_tsv %>% .[["trace_factor"]]
  gamma <- glue("misc/out_o2/h2/{trait}.{ntop}.tsv.gz") %>% read_tsv %>% .[["trace_factor"]]
  cat(" - gamma", gamma, "\n")

  t1 <- glue("misc/out_o2/assoc_lm/{trait}.{ntop}.tsv.gz") %>% read_tsv
  t1 <- select(t1, snp, se, zscore, pval) %>% 
    dplyr::rename(se_lr = se, z_lr = zscore, p_lr = pval)

  t2 <- lapply(vals_chr, function(chr) {
    f <- glue("misc/out_o2/assoc_lmm_chr/{chr}.{trait}.{ntop}.tsv.gz") 
    `if`(file.exists(f), read_tsv(f), tibble())
  }) %>% bind_rows
  # stopifnot(nrow(t2) == ntop)

  t2 <- select(t2, snp, se, zscore, pval) %>% 
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
ptab <- tab %>% filter(est == "median")

offset <- 0.9
# ylims <- c(0, 0.5)
p <- ggplot(ptab, aes(filter, gamma - offset, fill = estimator)) + 
  geom_bar(stat = "identity", position = position_dodge()) + 
  geom_errorbar(aes(ymin = q25 - offset, ymax = q75 - offset), 
    width = 0.3, position = position_dodge(0.9)) +
  geom_text(aes(x = filter, y = 0, label = m), vjust = -1.5, size = 3.5) +
  geom_hline(yintercept = 1 - offset, linetype = 3, color = "grey20")

p <- p + facet_wrap(est ~ trait, scales = "free", ncol = 3)

p <- p + 
  scale_y_continuous(labels = function(x) x + offset) +
  theme(legend.position = "top") +
    # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(x = NULL, y = NULL)

ggsave("tmp.png", plot = p, dpi = 100) #, width = 10, height = 6)


