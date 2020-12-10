library(tidyverse)
library(cowplot)
theme_set(theme_minimal(8))

library(devtools)
load_all("~/git/variani/biglmmz/")
load_all("~/git/variani/matlm/")

N <- 2000; M <- 100; h2 <- 0.8

# simulate genotypes
set.seed(33)
freqs <- rep(0.5, M) # allele freq. = 0.5
Z <- sapply(freqs, function(f) rbinom(N, 2, f)) 

Z_names <- paste0("z", seq(M))
colnames(Z) <- Z_names 

# select 2/3 for GRM, 1/3 to testing
z_grm <- seq(floor(0.75*M))
z_test <- seq(max(z_grm), M)
# z_test <- z_grm

# simulate data
Z_sc <- scale(Z)

b_true <- rnorm(M, 0, sqrt(h2/M))
y <- Z_sc %*% b_true + rnorm(N, 0, sqrt(1 - h2))

y_sc <- scale(y)

# fit model on scaled genotypes and normalized by sqrt(M)
m1 <- biglmmz(y, Z = Z, scale = TRUE)

Zgrm <- Z_sc[, z_grm]
Zgrm <- Zgrm / sqrt(ncol(Zgrm))

m2 <- biglmmz(y, Z = Zgrm, scale = FALSE)

h2 <- data.frame(model = c("all casual", "75% causal"),
  h2 = c(m1$gamma, m2$gamma))

s2 <- m2$s2
gamma <- m2$gamma

## assoc by LR
Z_test <- Z_sc[, z_test]
assoc1 <- matlm(y ~ 1, data.frame(y = y_sc), pred = Z_test, stats_full = TRUE)$tab %>% arrange(pval)

# assoc10 <- lapply(colnames(Z_test), function(z) {
#   x_sc <- scale(Z_test[, z])
#   mod <- lm(y_sc ~ x_sc)
#   coef <- broom::tidy(mod) %>% tail(1) 
#   with(coef, tibble(predictor = z, b = estimate, se = std.error, pval = p.value))
# }) %>% bind_rows
# assoc1 <- assoc10

## assoc by low-rank LMM (no EVD)
yc <- y - mean(y)
K <- crossprod(Zgrm)

assoc3 <- lapply(z_test, function(i) {
  X_sc <- scale(Z[, i, drop = FALSE])
  est <- biglmmz:::biglr_fixef(gamma = gamma, s2 = s2,
    y = y_sc, Xmat = X_sc, Z = Zgrm, K = K) # , REML = FALSE)
  tibble(predictor = Z_names[i], 
    b = est$b, se = sqrt(diag(est$bcov)))
}) %>% bind_rows %>%
  mutate(zscore = b / se, pval = pchisq(zscore*zscore, df = 1, lower = FALSE)) %>% 
  arrange(pval)

## assoc by low-rank LMM by matrix computation: explicit scaling 
X <- Z_test
X_sc <- scale(X)
comp <- s2 * c(gamma, 1 - gamma)

XVt <- biglr_cprodMatInv(comp, Zgrm, X_sc, K, transpose = TRUE) # Vi' X
XVX <- colSums(XVt * X_sc) # https://stackoverflow.com/a/21708690
b <- as.numeric(crossprod(XVt, y_sc)) / XVX
se <- 1 / sqrt(XVX)
assoc4 <- tibble(predictor = colnames(Z_test), b = b, se = se) %>%
  mutate(zscore = b / se, pval = pchisq(zscore*zscore, df = 1, lower = FALSE)) %>% 
  arrange(pval)

## assoc by low-rank LMM by matrix computation: auto scaling 
Zgrm_fbm <- as_FBM(Z[, z_grm])
XVt <- biglr_cprodMatInv2(comp, Zgrm_fbm, X_sc, transpose = TRUE) # Vi' X
XVX <- colSums(XVt * X_sc) 
b <- as.numeric(crossprod(XVt, y_sc)) / XVX
se <- 1 / sqrt(XVX)
assoc5 <- tibble(predictor = colnames(Z_test), b = b, se = se) %>%
  mutate(zscore = b / se, pval = pchisq(zscore*zscore, df = 1, lower = FALSE)) %>% 
  arrange(pval)

## assoc by low-rank LMM by matrix computation: auto scaling + pre-computed K
K <- crossprod(Zgrm)
XVt <- biglr_cprodMatInv2(comp, Zgrm_fbm, X_sc, K, transpose = TRUE) # Vi' X
XVX <- colSums(XVt * X_sc) 
b <- as.numeric(crossprod(XVt, y_sc)) / XVX
se <- 1 / sqrt(XVX)
assoc6 <- tibble(predictor = colnames(Z_test), b = b, se = se) %>%
  mutate(zscore = b / se, pval = pchisq(zscore*zscore, df = 1, lower = FALSE)) %>% 
  arrange(pval)
  
## scatter plot of p-values: assoc1 (LR) vs assoc2 (LMM)
ptab <- assoc1 %>% select(predictor, b, pval) %>% rename(pval_LR = pval, b_LR = b) %>%
  left_join(assoc3 %>% select(predictor, pval) %>% rename(pval_LMM_LowRank = pval)) %>%
  left_join(assoc4 %>% select(predictor, pval) %>% rename(pval_LMM_Matrix = pval)) %>% 
  left_join(assoc5 %>% select(predictor, pval) %>% rename(pval_LMM_Matrix_Sc = pval)) %>%
  left_join(assoc6 %>% select(predictor, b, pval) %>% rename(pval_LMM_Matrix_Sc_K = pval, b_LMM = b))

ptab <- left_join(ptab,
  tibble(predictor = colnames(Z), b_true = b_true))

pb1 <- ggplot(ptab, aes(abs(b_true), abs(b_LR))) + geom_point() + geom_abline(linetype = 3)
pb2 <- ggplot(ptab, aes(abs(b_true), abs(b_LMM))) + geom_point() + geom_abline(linetype = 3)
pb3 <- ggplot(ptab, aes(abs(b_LR), abs(b_LMM))) + geom_point() + geom_abline(linetype = 3)

p1 <- ggplot(ptab, aes(-log10(pval_LR), -log10(pval_LMM_LowRank))) + geom_point() + geom_abline(linetype = 3)
p2 <- ggplot(ptab, aes(-log10(pval_LR), -log10(pval_LMM_Matrix))) + geom_point() + geom_abline(linetype = 3)
p3 <- ggplot(ptab, aes(-log10(pval_LMM_LowRank), -log10(pval_LMM_Matrix))) + geom_point() + geom_abline(linetype = 3)
p4 <- ggplot(ptab, aes(-log10(pval_LMM_Matrix), -log10(pval_LMM_Matrix_Sc))) + geom_point() + geom_abline(linetype = 3)
p5 <- ggplot(ptab, aes(-log10(pval_LMM_Matrix), -log10(pval_LMM_Matrix_Sc_K))) + geom_point() + geom_abline(linetype = 3)

g <- plot_grid(
  pb1, pb2, pb3,
  p1, p2, NULL,
  p3, p4, p5, labels = "auto")
ggsave("tmp.png", plot = g, dpi = 100)
