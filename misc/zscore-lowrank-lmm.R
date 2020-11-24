library(matlm)
library(tidyverse)
library(cowplot)
theme_set(theme_minimal())

library(devtools)
load_all("~/git/variani/biglmmz/")

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
# z_test <- z_grm

# simulate data
Z_sc <- scale(Z)

b <- rnorm(M, 0, sqrt(h2/M))
y <- Z_sc %*% b + rnorm(N, 0, sqrt(1 - h2))

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
Z_test <- Z[, z_test]
assoc1 <- matlm(y ~ 1, data.frame(y = y), pred = Z_test)$tab %>% arrange(pval)

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
ptab <- assoc1 %>% select(predictor, pval) %>% rename(pval_LR = pval) %>%
  left_join(assoc3 %>% select(predictor, pval) %>% rename(pval_LMM_LowRank = pval)) %>%
  left_join(assoc4 %>% select(predictor, pval) %>% rename(pval_LMM_Matrix = pval)) %>% 
  left_join(assoc5 %>% select(predictor, pval) %>% rename(pval_LMM_Matrix_Sc = pval)) %>%
  left_join(assoc6 %>% select(predictor, pval) %>% rename(pval_LMM_Matrix_Sc_K = pval))

p1 <- ggplot(ptab, aes(-log10(pval_LR), -log10(pval_LMM_LowRank))) + geom_point() + geom_abline(linetype = 3)
p2 <- ggplot(ptab, aes(-log10(pval_LR), -log10(pval_LMM_Matrix))) + geom_point() + geom_abline(linetype = 3)
p3 <- ggplot(ptab, aes(-log10(pval_LMM_LowRank), -log10(pval_LMM_Matrix))) + geom_point() + geom_abline(linetype = 3)
p4 <- ggplot(ptab, aes(-log10(pval_LMM_Matrix), -log10(pval_LMM_Matrix_Sc))) + geom_point() + geom_abline(linetype = 3)
p5 <- ggplot(ptab, aes(-log10(pval_LMM_Matrix), -log10(pval_LMM_Matrix_Sc_K))) + geom_point() + geom_abline(linetype = 3)

g <- plot_grid(p1, p2, p3, p4, p5, labels = "auto")
