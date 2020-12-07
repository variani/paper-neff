library(tidyverse)
library(glue)
library(data.table)
library(RColorBrewer) 

library(cowplot)
library(latex2exp)

### settings
theme_set(theme_minimal())

# plotting par.
offset <- 0.9
ylims <- c(0, 0.35) 
cols_accent8 <- brewer.pal(8, "Accent")

filters_plot <- c("f0", "f2")
labs_filters <- c("Empirical (all selected variants)",
  "Empirical (top selected variants)", "Analytical")
vals_traits <- c("bmi", "height", "hip", "waist", "weight", "whr")
labs_traits <- c("BMI", "Height", "HIP", "Waist", "Weight", "WHR")
# vals_traits <- c("height", "bmi")
ntop <- 1000
vals_chr <- 1:22

vals_est <- c("mean", "median")
vals_filt <- paste0("f", 0:3)

thr2 <- 1e-5 
# thr2 <- 5e-8
thr1_lmm <- thr2
thr1_lr <- 0.05
thr3 <- thr2

tab <- lapply(vals_traits, function(trait) {
  cat("trait", trait, "\n")

  h2 <- glue("out/h2/{ntop}/{trait}.tsv.gz") %>% read_tsv 
  if(length(vals_chr) == 1) {
    h2 <- glue("out/h2-loco/{ntop}/{trait}.{vals_chr}.tsv.gz") %>% read_tsv
  }
  gamma <- h2$mult
  cat(" - gamma", gamma, "\n")

  # f1 <- list.files("misc/out_o2/assoc_lm/", trait, ignore.case = TRUE, full = TRUE)
  # t1 <- fread(f1) %>% as_tibble
  # t1 <- mutate(t1, se = 1/sqrt(N))

  f1 <- glue("out/lm_pcs_top/{trait}.tsv.gz")
  t1 <- fread(f1) %>% as_tibble

  t1 <- select(t1, snp, se, zscore, pval) %>% 
    dplyr::rename(se_lr = se, z_lr = zscore, p_lr = pval)
  print(t1)

  t2 <- glue("out/lmm_loco_top/{ntop}/{trait}.{vals_chr}.tsv.gz") %>%
    lapply(fread) %>% bind_rows %>% as_tibble
  t2 <- select(t2, snp, se, zscore, pval) %>% 
    dplyr::rename(se_lmm = se, z_lmm = zscore, p_lmm = pval)

  stopifnot(all(t2$snp %in% t1$snp))
  t <- left_join(t2, t1, by = "snp")

  lapply(vals_filt, function(filt) {
    cat("filter", filt, "\n")

    t <- switch(filt, 
      "f0" = t,
      "f1" = filter(t, p_lr < thr1_lr & p_lmm < thr1_lmm),
      "f2" = filter(t, p_lr < thr2 & p_lmm < thr2),
      "f3" = filter(t, p_lr < thr3 & p_lmm < thr3),
      stop("filt"))

    vals_se2 <- with(t, (se_lr / se_lmm)^2)
    vals_z2 <- with(t, (z_lmm / z_lr)^2)

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

# plot panel a: gamma_se
ptab <- filter(tab, filter %in% filters_plot & est == "median") %>%
  mutate(filter = ifelse(estimator == "trace", "f", filter)) %>%
  filter(estimator %in% c("trace", "se2")) %>%
  group_split(trait) %>% lapply(function(t) filter(t, !duplicated(filter))) %>% bind_rows %>%
  mutate(
    trait = factor(trait, vals_traits, labs_traits),
    filter = factor(filter, c(filters_plot, "f"), labs_filters))

pa <- ggplot(ptab, aes(trait, gamma - offset, fill = filter, group = filter)) + 
  geom_bar(stat = "identity", position = position_dodge(0.9)) + 
  geom_errorbar(aes(ymin = q25 - offset, ymax = q75 - offset, group = filter), 
    width = 0.3, position = position_dodge(0.9))

# format y labels as x1.1
pa <- pa + scale_y_continuous(limits = ylims, labels = function(x) 
  paste0(formatC(x + offset, format = "f", digits = 1), "x"))
    
# configure colors, labels, theme
cols <- cols_accent8[c(8, 3, 6)]
labx <- NULL
laby <-  "ESS Multiplier"
title <- TeX("$\\gamma_{\\beta}^{se}$ $\\mathit{vs.}$ $\\gamma_{\\beta}$")

pa <- pa + scale_fill_manual(values = cols) + 
  theme(legend.position = "none", strip.background = element_blank(), 
    strip.text.x = element_blank(), panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), 
    axis.title.x = element_blank()) + 
  labs(x = labx, y = laby, subtitle = title)

## plot panel b: gamma_zw
ptab <- filter(tab, filter %in% filters_plot & est == "median") %>%
  mutate(filter = ifelse(estimator == "trace", "f", filter)) %>%
  filter(estimator %in% c("trace", "z2")) %>%
  group_split(trait) %>% lapply(function(t) filter(t, !duplicated(filter))) %>% bind_rows %>%
  mutate(
    trait = factor(trait, vals_traits, labs_traits),
    filter = factor(filter, c(filters_plot, "f"), labs_filters))

pb <- ggplot(ptab, aes(trait, gamma - offset, fill = filter, group = filter)) + 
  geom_bar(stat = "identity", position = position_dodge(0.9)) + 
  geom_errorbar(aes(ymin = q25 - offset, ymax = q75 - offset, group = filter), 
    width = 0.3, position = position_dodge(0.9))

pb <- pb + scale_y_continuous(limits = ylims, labels = function(x) 
  paste0(formatC(x + offset, format = "f", digits = 1), "x"))

title <- TeX("$\\gamma_{\\beta}^{s}$ $\\mathit{vs.}$ $\\gamma_{\\beta}$")
pb <- pb + theme(legend.position = "none", strip.background = element_blank(), strip.text.x = element_blank(), panel.grid.minor = element_blank(), panel.grid.major.x = element_blank()) + 
  labs(x = NULL, y = "", fill = "Multiplier:  ", subtitle = title)
  
pb <- pb + scale_fill_manual(values = cols) + 
  theme(legend.position = "bottom", legend.justification = "top", axis.title.x = element_blank(), 
    axis.text.y = element_blank()) + 
  guides(fill = guide_legend(nrow = 1)) 

## combine panels a + b
pab <- plot_grid(pa, pb + theme(legend.position = "none"), 
  labels = "auto", rel_widths = c(1.1, 1))

## legend
pl <- pb
pl <- pl + theme(legend.text = element_text(margin = margin(l = -1, r = 15)),
  legend.key.height = unit(0.5, "cm"), legend.key.width = unit(0.3, "cm"))
legend <- get_legend(pl + theme(legend.position = "bottom"))

## combine all plots & save
g <- plot_grid(pab, legend, ncol = 1, rel_widths = c(1, 1), rel_heights = c(1, 0.05))

width <- 8
height <- (7/16) * width
ggplot2::ggsave("tmp.png", g, width = width, height = height)
