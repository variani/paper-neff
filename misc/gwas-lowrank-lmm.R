library(biglmmz)
library(matlm)
library(tidyverse)
library(cowplot)
theme_set(theme_minimal())

N <- 1500; M <- 200; h2 <- 0.8

# simulate genotypes
set.seed(33)
freqs <- rep(0.5, M) # allele freq. = 0.5
Z <- sapply(freqs, function(f) rbinom(N, 2, f)) 

Z_names <- paste0("z", seq(M))
colnames(Z) <- Z_names 

# select 2/3 for GRM, 1/3 to testing
z_grm <- seq(floor(0.75*M))
z_test <- seq(max(z_grm), M)

# simulate data
Z_means <- colMeans(Z, na.rm = TRUE)
Z_freq <- Z_means / 2  # Z_means = 2 * Z_freq
Z_sd <- sqrt(2 * Z_freq * (1 - Z_freq))

Z_sc <- sweep(Z, 2, Z_means, "-")
Z_sc <- sweep(Z_sc, 2, Z_sd , "/")

b <- rnorm(M, 0, sqrt(h2/M))
y <- Z_sc %*% b + rnorm(N, 0, sqrt(1 - h2))

# fit model on scaled genotypes and normalized by sqrt(M)
m1 <- biglmmz(y, Z = Z, scale = TRUE)

Zgrm <- Z[, z_grm]
Zgrm <- sweep(Zgrm, 2, Z_means[z_grm], "-")
Zgrm <- sweep(Zgrm, 2, Z_sd[z_grm], "/")
Zgrm <- Zgrm / sqrt(M)

m2 <- biglmmz(y, Z = Zgrm, scale = FALSE)

h2 <- data.frame(model = c("all casual", "75% causal"),
  h2 = c(m1$gamma, m2$gamma))

## assoc by LR
Z_test <- Z[, z_test]
assoc1 <- matlm(y ~ 1, data.frame(y = y), pred = Z_test)$tab %>% arrange(pval)

## assoc by LMM (EVD) 
# - prohibitive at large N, but here N = 1500, so we can stor NxN GRM
GRM <- tcrossprod(Zgrm)
h2 <- m2$gamma
V <- h2*GRM + (1-h2)*diag(N)

assoc2 <- matlm(y ~ 1, data.frame(y = y), varcov = V, pred = Z_test)$tab %>% arrange(pval)

## assoc by low-rank LMM (no EVD)
yc <- y - mean(y)
K <- crossprod(Zgrm)
assoc3 <- lapply(z_test, function(i) {
  Xc <- matrix(Z[, i] - Z_means[i], ncol = 1)
  est <- biglmmz:::biglr_fixef(gamma = m2$gamma, s2 = m2$s2,
    y = y, Xmat = Xc, Z = Zgrm, K = K, REML = TRUE)
  tibble(predictor = Z_names[i], zscore = est$b / sqrt(diag(est$bcov)))
}) %>% bind_rows %>%
  mutate(pval = pchisq(zscore*zscore, df = 1, lower = FALSE)) %>% 
  arrange(pval)

## scatter plot of p-values: assoc1 (LR) vs assoc2 (LMM)
ptab <- assoc1 %>% select(predictor, pval) %>% rename(pval_LR = pval) %>%
  left_join(assoc2 %>% select(predictor, pval) %>% rename(pval_LMM = pval)) %>%
  left_join(assoc3 %>% select(predictor, pval) %>% rename(pval_LMM_LowRank = pval))

p1 <- ggplot(ptab, aes(-log10(pval_LR), -log10(pval_LMM))) + geom_point() + geom_abline(linetype = 3)
p2 <- ggplot(ptab, aes(-log10(pval_LR), -log10(pval_LMM_LowRank))) + geom_point() + geom_abline(linetype = 3)
p3 <- ggplot(ptab, aes(-log10(pval_LMM), -log10(pval_LMM_LowRank))) + geom_point() + geom_abline(linetype = 3)

g <- plot_grid(p1, p2, p3)

