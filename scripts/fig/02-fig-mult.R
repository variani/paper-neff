library(tidyverse)
library(glue)

library(cowplot)
theme_set(theme_cowplot(14))

vals_traits <- c("bmi", "weight", "waist", "hip", "height", "whr")
# vals_traits <- "height"
ntop <- 1000
vals_chr <- 1:22

vals_est <- c("mean", "median")
vals_filt <- paste0("f", 0:3)

thr2 <- c(1e-3, 1e-5, 5e-8)[3]
thr1_lmm <- thr2
thr1_lr <- 0.05
thr3 <- thr2

tab <- lapply(vals_traits, function(trait) {
  cat("trait", trait, "\n")

  h2 <- glue("out/h2/{ntop}/{trait}.tsv.gz") %>% read_tsv 
  # h2 <- glue("out/h2-loco/{ntop}/{trait}.{vals_chr}.tsv.gz") %>% read_tsv
  gamma <- h2$mult
  cat(" - gamma", gamma, "\n")

  # f1 <- list.files("misc/out_o2/assoc_lm/", trait, ignore.case = TRUE, full = TRUE)
  f1 <- glue("out/lm_pcs_top/{trait}.tsv.gz")
  t1 <- read_tsv(f1)
  t1 <- mutate(t1, se = 1/sqrt(N))
  t1 <- select(t1, snp, beta, se, zscore, pval) %>% 
    dplyr::rename(b_lr = beta, se_lr = se, z_lr = zscore, p_lr = pval)

  # t2 <- glue("out/lmm_loco_top/{ntop}/{trait}.{vals_chr}.tsv.gz") %>%
  t2 <- glue("out/lmm_loco_pcs_top/{ntop}/{trait}.{vals_chr}.tsv.gz") %>%
    lapply(read_tsv) %>% bind_rows
  t2 <- select(t2, snp, beta, se, zscore, pval) %>% 
    dplyr::rename(b_lmm = beta, se_lmm = se, z_lmm = zscore, p_lmm = pval)

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
          q25 = c(NA, quantile(vals_se2, 0.25, na.rm = TRUE), quantile(vals_z2, 0.25, na.rm = TRUE)),
          q75 = c(NA, quantile(vals_se2, 0.75, na.rm = TRUE), quantile(vals_z2, 0.75, na.rm = TRUE))) %>%
        mutate(trait = trait, filter = filt, m = nrow(t), est = est)
      }) %>% bind_rows
  }) %>% bind_rows
}) %>% bind_rows

## plot
# NB: filter se2 out
# ptab <- filter(tab, estimator != "se2")
ptab <- tab
# ptab <- filter(ptab, !(estimator == "z2" & filter == "f0"))
# ptab <- filter(ptab, filter != "f0" & est == "median")
ptab <- filter(ptab, est == "median")

offset <- 0.9
# ylims <- c(0, 0.5)
p <- ggplot(ptab, aes(filter, gamma - offset, fill = estimator)) + 
  geom_bar(stat = "identity", position = position_dodge()) + 
  geom_errorbar(aes(ymin = q25 - offset, ymax = q75 - offset), 
    width = 0.3, position = position_dodge(0.9)) +
  geom_text(aes(x = filter, y = 0, label = m), vjust = -1.5, size = 3.5) +
  geom_hline(yintercept = 1 - offset, linetype = 3, color = "grey20")

p <- p + facet_wrap(est ~ trait, scales = "free", ncol = 3)

# p <- ggplot(ptab, aes(trait, gamma - offset, fill = estimator, group = estimator)) + 
#   geom_bar(stat = "identity", position = position_dodge(0.9)) + 
#   geom_errorbar(aes(ymin = q25 - offset, ymax = q75 - offset, group = estimator), 
#     width = 0.3, position = position_dodge(0.9)) +
#   geom_hline(yintercept = 1 - offset, linetype = 3, color = "grey")

# p <- p + facet_grid(est ~ filter)

p <- p + 
  scale_y_continuous(labels = function(x) x + offset) +
  theme(legend.position = "top") +
  labs(x = NULL, y = NULL)

ggsave("tmp.png", plot = p, dpi = 100, width = 12, height = 6)

