library(tidyverse)
library(glue)
library(cowplot)
theme_set(theme_minimal())

library(devtools)
load_all("~/git/variani/bigcov") # split bed by batches
load_all("~/git/variani/matlm") # split bed by batches
load_all("~/git/variani/biglmmz") # split bed by batches

N <- 10e3; h2 <- 0.5
Mc <- 500; M <- 1000
Mchr <- M / 2

# simulate genotypes
# set.seed(1)
freqs <- rep(0.5, M) # allele freq. = 0.5
Z <- sapply(freqs, function(f) rbinom(N, 2, f)) 

Z_names <- paste0("z", seq(M))
colnames(Z) <- Z_names 

Z_sc <- scale(Z)

# simulate data
z_c <- sample(seq(M), Mc)

b <- rnorm(Mc, 0, sqrt(h2/Mc))
y <- Z_sc[, z_c] %*% b + rnorm(N, 0, sqrt(1 - h2))

y_sc <- scale(y)

# fit model on scaled genotypes and normalized by sqrt(M)
z_grm <- seq(Mchr)
z_grm <- z_grm[z_grm %in% z_c]

Mgrm <- length(z_grm)
Zgrm <- Z[, z_grm] 

mod <- biglmmz(y, Z = Zgrm, scale = TRUE)

## assoc by LR
z_test <- seq(M)
z_test <- z_test[!(z_test %in% z_grm)]

Mtest <- length(z_test)
Z_test_sc <- scale(Z[, z_test])

assoc1 <- matlm(y ~ 1, data.frame(y = y_sc), pred = Z_test_sc, stats_full = TRUE) %>%
  .[["tab"]] %>% arrange(pval) %>% mutate(causal = (predictor %in% Z_names[z_c]))

## Z as FBM
file_z <- tempfile("misc-") %>% basename
file_z_bk <- glue("{file_z}.bk")
Zgrm_fbm <- as_FBM(Zgrm, type = "integer", backingfile = file_z) 

# pre-compute K
cat(" - pre-compute K\n")
K <- big_crossprodSelf(Zgrm_fbm, fun.scaling = big_scale2(M = Mgrm))[]

gamma <- mod$gamma
s2 <- mod$s2
comp <- s2 * c(gamma, 1 - gamma)
# comp <- c(gamma, 1 - gamma)

## assoc by LMM
beg <- seq(1, Mtest, 100)
end <- c(beg[-1] - 1, Mtest)
nb <- length(beg)

assoc2 <- lapply(seq(nb), function(b) {
  cat(" -", b, "/", nb, "\n")

  cols <- seq(beg[b], end[b])
  Xmat <- Z_test_sc[, cols, drop = FALSE]

  XVt <- biglr_cprodMatInv2(comp, Zgrm_fbm, X = Xmat, K = K, transpose = TRUE) # Vi' X
  XVX <- colSums(XVt * Xmat) 
  b <- as.numeric(crossprod(XVt, y_sc)) / XVX
  se <- 1 / sqrt(XVX)
  assoc <- tibble(predictor = colnames(Xmat), beta = b, se = se) %>%
    mutate(zscore = beta / se, 
      pval = pchisq(zscore*zscore, df = 1, lower = FALSE))
  assoc
}) %>% bind_rows
unlink(file_z_bk)

assoc2 <- arrange(assoc2, pval) %>% mutate(causal = (predictor %in% Z_names[z_c]))

## scatter plot of p-values: assoc1 (LR) vs assoc2 (LMM)
ptab1 <- rename(assoc1, b_LR = b, se_LR = se, p_LR = pval)
ptab2 <- rename(assoc2, b_LMM = beta, se_LMM = se, p_LMM = pval)
ptabc <- tibble(predictor = Z_names[z_c], b_causal = b)

ptab <- left_join(ptab1, ptab2, by = "predictor") %>% left_join(ptabc, by = "predictor")
ptab_signif <- filter(ptab, p_LR < 5e-8 & p_LMM < 5e-8)

p01 <- ggplot(ptab, aes(abs(b_causal), abs(b_LR))) + geom_point() + geom_abline(linetype = 3) + geom_smooth(method = "lm")
p02 <- ggplot(ptab, aes(abs(b_causal), abs(b_LMM))) + geom_point() + geom_abline(linetype = 3) + geom_smooth(method = "lm")

lims <- c(0.05, 0.1)
up_p <- function(p) p + xlim(lims) + ylim(lims) + geom_smooth(method = "lm", color = "red") 

p1 <- ggplot(ptab, aes(abs(b_LR), abs(b_LMM))) + geom_point() + geom_abline(linetype = 3) + geom_smooth(method = "lm")
p2 <- ggplot(ptab, aes(1/se_LR^2, 1/se_LMM^2)) + geom_point() + geom_abline(linetype = 3)
p3 <- ggplot(ptab, aes(-log10(p_LR), -log10(p_LMM))) + geom_point() + geom_abline(linetype = 3)

g <- plot_grid(p01, p02, p1, 
  up_p(p01), up_p(p02), up_p(p1),
  p2, p3, labels = "auto")
ggsave("tmp.png", plot = g, dpi = 100) 

