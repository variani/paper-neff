library(tidyverse)
library(glue)

library(cowplot)
theme_set(theme_cowplot())

library(matlm)

library(devtools)
load_all("~/git/variani/biglmmz")

N <- 5000; M <- 2000; Mc <- 500
Mnull <- 100
Mgrm <- Mc + Mnull
h2 <- 0.8

vals_est <- c("median", "mean")
vals_filt <- paste0("f", 0:3)
thr2 <- c(1e-3, 1e-5, 5e-8)[3]
thr1_lmm <- thr2
thr1_lr <- 0.05
thr3 <- thr2

# simulate genotypes
# set.seed(33)
freqs <- rep(0.5, M) # allele freq. = 0.5
Z <- sapply(freqs, function(f) rbinom(N, 2, f)) 

Z_names <- paste0("z", seq(M))
colnames(Z) <- Z_names 

z_c <- seq(Mc)
z_grm <- seq(Mgrm)
z_test <- seq(M)

# simulate data
Z_sc <- scale(Z)

b <- rnorm(Mc, 0, sqrt(h2/Mc))
y <- Z_sc[, z_c] %*% b + rnorm(N, 0, sqrt(1 - h2))

y_sc <- scale(y)

## assoc by LR
assoc1 <- matlm(y ~ 1, data.frame(y = y_sc), pred = Z_sc, stats_full = TRUE)$tab %>% arrange(pval)

## assoc by LMM
# step 1: null model
mod <- biglmmz(y_sc, Z = Z[, z_grm], impute = FALSE, scale = TRUE, verbose = 2)
 
# extract estimates of model parameters (variance components)
gamma <- mod$gamma
s2 <- mod$s2

## assoc by LMM
# step 2: assoc
Zgrm <- Z_sc[, z_grm] / sqrt(length(z_grm) - 1)
K <- crossprod(Zgrm)

assoc2 <- lapply(z_grm, function(j) {
  cat(" -", j, "/", length(z_grm), "\n")
  # matrix X (covariates)
  g <- Z_sc[, j]
  g <- scale(g, center = TRUE, scale = TRUE)
  
  X0 <- matrix(1, nrow = length(g), ncol = 1)
  X <- cbind(X0, g)
  
  # LMM: snp j is excluded from GRM
  est <- biglmmz:::biglr_fixef(gamma = gamma, s2 = s2,
    y = y_sc, Xmat = X, Z = Zgrm[, -j], K = K[-j, -j], REML = TRUE)
  coef <- data.frame(beta = est$b, se = sqrt(diag(est$bcov)))
  coef <- tail(coef, 1) # SNP goes last

  # mod <- biglmmz(y_sc, X = X, Z = Zgrm[, -j], K = K[-j, -j])
  # coef <- tail(mod$coef, 1) # the last raw for SNP
  
  tibble(predictor = Z_names[j], 
    b = coef$beta, se = coef$se) %>%
    mutate(zscore = b / se)
}) %>% bind_rows %>%
  mutate(pval = pchisq(zscore*zscore, df = 1, lower = FALSE)) %>% 
  arrange(pval)

## scatter plot of p-values: assoc1 (LR) vs assoc2 (LMM)
tab <- assoc1 %>% select(predictor, se, zscore, pval) %>% 
    rename(se_lr = se, z_lr = zscore, pval_lr = pval) %>%
  left_join(assoc2 %>% select(predictor, se, zscore, pval) %>% 
    rename(se_lmm = se, z_lmm = zscore, pval_lmm = pval)) 

tab <- filter(tab, !is.na(pval_lmm))

## ESS
gamma <- mod$ess$trace_factor
h2_hat <- mod$gamma

ess <- lapply(vals_filt, function(filt) {
  tab <- switch(filt,
    "f0" = tab,
    "f1" = filter(tab, pval_lr < thr1_lr & pval_lmm < thr1_lmm),
    "f2" = filter(tab, pval_lr < thr2 & pval_lmm < thr2),
    "f3" = filter(tab, pval_lr < thr3 & pval_lmm < thr3),
    stop("filt"))
  mc <- sum(as.integer(gsub("^z", "", tab$predictor)) %in% z_c)

  lapply(vals_est, function(est) {
    vals_se2 <- with(tab, (se_lr / se_lmm)^2)
    vals_z2 <- with(tab, (z_lmm / z_lr)^2)

    se2 <- switch(est, "median" = median(vals_se2), "mean" = mean(vals_se2), stop("est"))
    z2 <- switch(est, "median" = median(vals_z2), "mean" = mean(vals_z2), stop("est"))
    
    tibble(gamma = c(gamma, se2, z2),
        estimator = c("trace", "se2", "z2"),
        q25 = c(NA, quantile(vals_se2, 0.25), quantile(vals_z2, 0.25)),
        q75 = c(NA, quantile(vals_se2, 0.75), quantile(vals_z2, 0.75))) %>%
      mutate(filt = filt, h2_hat = h2_hat, m = nrow(tab), mc = mc, 
        est = est, trait = "simulated")
  }) %>% bind_rows
}) %>% bind_rows

## plot 1: p-values
p1 <- ggplot(tab, aes(-log10(pval_lr), -log10(pval_lmm))) + 
  geom_point() + geom_abline(linetype = 3) +
  geom_vline(xintercept = -log10(thr2), linetype = 3) +
  geom_hline(yintercept = -log10(thr2), linetype = 3) +
  coord_equal()

p1 <- p1 + labs(subtitle = glue("h2 = {h2}, N = {N}, M = {M}, Mc = {Mc}", "\n",
    "Mnull in GRM = {Mnull}, Mgrm = {Mgrm}, Mtest = Mgrm"))

## plot 2: ESS
# ptab <- filter(ess, est == "median")
ptab <- ess

ptab <- mutate(ptab, gamma = ifelse(estimator == "z2" & est == "mean" & filt == "f0",
    NA, gamma))

ptab <- mutate(ptab, 
  h2_hat = as.factor(paste0("h2_hat = ", round(h2_hat, 2))),
  filt = as.factor(paste0(filt, " [", mc, "/", m, "]")))

p2 <- ggplot(ptab, aes(h2_hat, gamma, fill = estimator, group = estimator)) +
  geom_bar(stat = "identity", position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = q25, ymax = q75, group = estimator),
    width = 0.3, position = position_dodge(0.9)) +
  geom_hline(yintercept = 1, linetype = 2, color = "grey40") +
  # facet_wrap(est ~ filt, scale = "free_y", nrow = 2)
  facet_wrap(est ~ filt, nrow = 2) +
  labs(x = NULL) +
  theme(legend.position = "top")
 
## combind plots
g <- plot_grid(p1, p2, labels = "auto")

ggsave("tmp.png", plot = g, dpi = 100)
