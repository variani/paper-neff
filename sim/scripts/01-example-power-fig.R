## inc
library(devtools)
load_all("~/git/hemostat/gxefam")
load_all("~/git/variani/matlm")
load_all("~/git/variani/qq")
library(dplyr)
library(magrittr)
library(tidyverse)
library(ggplot2)
library(cowplot)
theme_set(theme_minimal())

compute_power <- function(ncp, alpha = 0.05, df = 1, lower = FALSE)
{
  q_alpha <- qchisq(alpha, df = df, lower = lower)
  power <- pchisq(q_alpha, df = df, ncp = ncp, lower = lower)
}

## par
M <- 1e3 # the number of SNPs
N <- as.integer(c(1e2, 5e2, 1e3)) # the number of individuals
p <- 0.3 # the MAF of SNPs
f <- 3 / 5 # the frequency of binary exposure
betaG <- 0.1 # the effect size (simulations on the main genetic effect)
var_explained <- 0.6 # the variance explained by polygenic effect (h2)
ibetaG <- 0.1 # the marginal genetic effect size (GxE simulations)
ibetaE <- 0.1 # the marginal environmental effect size
ibetaGE <- 0.1 # the interaction effect size
S <- 1e3 # number of replications
alpha <- 0.05

## parallel
cores <- parallel::detectCores()
cores <- ifelse(cores > 1, cores - 1, cores)
parallel <- (cores > 1)
if(parallel) { doParallel::registerDoParallel(cores) }

## main function
run_sim <- function(sim_data, est_power)
{
  cat(" - single simulation\n")
  tab <- do.call(sim_data, list())
  
  cat(" - multipler simulations in parallel (", cores, ")\n")
  tab_s <- plyr::llply(seq(S), function(s) {
    cat(" - s:", s, "\n")
    tab <- do.call(sim_data, list()) %>% mutate(s = s)
  }, .parallel = parallel) %>% bind_rows
  tab_power <- group_by(tab_s, N, s) %>% summarize(power = mean(pval < alpha))

  cat(" - analytical multiplier\n")
  tab_expected <- do.call(est_power, list())

  tab <- mutate(tab, N = factor(N))
  tab_power <- mutate(tab_power, N = factor(N))
  tab_expected <- mutate(tab_expected, N = factor(N))
  
  p_beta <- ggplot(tab, aes(N, beta)) + geom_boxplot() + geom_boxplot(data = tab_expected, aes(N, beta), color = "red")
  p_se <- ggplot(tab, aes(N, se)) + geom_boxplot() + geom_boxplot(data = tab_expected, aes(N, se), color = "red") + scale_y_log10() + labs(y = "se (log 10 scale)")
  p_z <- ggplot(tab, aes(N, (beta/se)^2)) + geom_boxplot() + geom_boxplot(data = tab_expected, aes(N, ncp), color = "red")
  p_power <- ggplot(tab_power, aes(N, power)) + geom_boxplot() + geom_boxplot(data = tab_expected, aes(N, power), color = "red") + ylim(c(0, 1))
  p_s2 <- ggplot(tab, aes(N, sigma2)) + geom_boxplot() + geom_boxplot(data = tab_expected, aes(N, sigma2), color = "red")
  p_gamma <- ggplot(tab, aes(N, gamma)) + geom_boxplot() + geom_boxplot(data = tab_expected, aes(N, gamma), color = "red")

  list(tab = tab, tab_power = tab_power, tab_expected = tab_expected,
   p = plot_grid(p_se, p_beta, p_z, p_s2, p_gamma, p_power, labels = "auto"))
}

# simulation 6: unrelated individuals with grouping- marginal genetic effect (betaG)
get_par <- function(n) gxefam_par(num_fam = n / 5, sim_gen = "rbinom", 
  sim_phen = "fam", freq_g = p)

sim_data <- function()
{
  plyr::llply(N, function(n) {
    par <- get_par(n)
    dat <- gxefam_sim_dat(par, n = M, beta = c(betaG, 0, 0))

    tr <- as.matrix(gxefam_transform(par))
    assoc <- with(dat, matreg(Y = Y, X = G, transform = tr))

    mutate(assoc$tab, N = as.integer(n), 
      Neff = sigma2 / (2*p*(1-p)*se^2), gamma = Neff/N, sim = "rel_gen")
  }) %>% bind_rows
}

est_power <- function()
  tibble(N = N, beta = betaG, sigma2 = 1) %>%
    mutate(gamma = sapply(N, function(n) gxefam_trace_factor(get_par(n), mult = "I")),
      se = 1 / sqrt(2 * p * (1 - p) * N * gamma), ncp = (beta / se)^2,
      power = compute_power(ncp))

out <- run_sim(sim_data, est_power)
ggsave("fam-marginal.png", out$p, dpi = 150)

# simulation 5: related in families, siblings exposedi (f = 3/5 = 0.6),
# interaction effect (ibetaGE), 
# phen. variance explained by GxE random effect is 50% of total varince explained
get_par <- function(n) gxefam_par(num_fam = n / 5, sim_gen = "rbinom_fam",
  sim_phen = "gen", sim_env = "exposed_sib", freq_g  = p, freq_e = f,
  var_explained = var_explained, var_geni = 0.5 * var_explained)

sim_data <- function()
{
  plyr::llply(N, function(n) {
    par <- get_par(n)
    dat <- gxefam_sim_dat(par, n = M, center = TRUE, 
      beta = c(ibetaG, ibetaE, ibetaGE))

    tr <- as.matrix(gxefam_transform(par))
    assoc <- with(dat, matreg(Y = Y, X = GE, Xlist = list(G, E), transform = tr))

    mutate(assoc$tab, N = as.integer(n), 
      Neff = sigma2 / (2*p*(1-p)*f*(1-f)*se^2), gamma = Neff/N, sim = "unrel_int")
  }) %>% bind_rows
}

est_power <- function()
  tibble(N = N, beta = ibetaGE, sigma2 = 1) %>%
    mutate(gamma = sapply(N, function(n) gxefam_trace_factor(get_par(n), mult = "Kdc")),
      se = 1 / sqrt(2 * p * (1 - p) * f * (1 - f) * N * gamma), 
      ncp = (beta / se)^2, power = compute_power(ncp))

out <- run_sim(sim_data, est_power)
ggsave("rel-int2.png", out$p, dpi = 150)

# simulation 4: related in families, siblings exposedi (f = 3/5 = 0.6),
# interaction effect (ibetaGE), 
# phen. variance explained by GxE random effect is zero
get_par <- function(n) gxefam_par(num_fam = n / 5, sim_gen = "rbinom_fam",
  sim_phen = "gen", sim_env = "exposed_sib", freq_g  = p, freq_e = f,
  var_explained = var_explained, var_geni = 0)

sim_data <- function()
{
  plyr::llply(N, function(n) {
    par <- get_par(n)
    dat <- gxefam_sim_dat(par, n = M, center = TRUE, 
      beta = c(ibetaG, ibetaE, ibetaGE))

    tr <- as.matrix(gxefam_transform(par))
    assoc <- with(dat, matreg(Y = Y, X = GE, Xlist = list(G, E), transform = tr))

    mutate(assoc$tab, N = as.integer(n), 
      Neff = sigma2 / (2*p*(1-p)*f*(1-f)*se^2), gamma = Neff/N, sim = "unrel_int")
  }) %>% bind_rows
}

est_power <- function()
  tibble(N = N, beta = ibetaGE, sigma2 = 1) %>%
    mutate(gamma = sapply(N, function(n) gxefam_trace_factor(get_par(n), mult = "Kdc")),
      se = 1 / sqrt(2 * p * (1 - p) * f * (1 - f) * N * gamma), 
      ncp = (beta / se)^2, power = compute_power(ncp))

out <- run_sim(sim_data, est_power)
ggsave("rel-int.png", out$p, dpi = 150)

# simulation 3: unrelated, interaction effect (ibetaGE)
get_par <- function(n) gxefam_par(num_fam = n / 5, sim_gen = "rbinom", 
  sim_phen = "no_rand", freq_g  = p, freq_e = f)

sim_data <- function()
{
  plyr::llply(N, function(n) {
    par <- get_par(n)
    dat <- gxefam_sim_dat(par, n = M, center = TRUE, 
      beta = c(ibetaG, ibetaE, ibetaGE))

    assoc <- with(dat, matreg(Y = Y, X = GE, Xlist = list(G, E)))

    mutate(assoc$tab, N = as.integer(n), 
      Neff = sigma2 / (2*p*(1-p)*f*(1-f)*se^2), gamma = Neff/N, sim = "unrel_int")
  }) %>% bind_rows
}

est_power <- function()
  tibble(N = N, beta = ibetaGE, sigma2 = 1) %>%
    mutate(gamma = 1,
      se = 1 / sqrt(2 * p * (1 - p) * f * (1 - f) * N * gamma), 
      ncp = (beta / se)^2, power = compute_power(ncp))

out <- run_sim(sim_data, est_power)
ggsave("unrel-int.png", out$p, dpi = 150)

# simulation 2: related in families - marginal genetic effect (betaG)
get_par <- function(n) gxefam_par(num_fam = n / 5, sim_gen = "rbinom_fam", 
  sim_phen = "gen", var_explained = var_explained, freq_g = p,
  var_explained = var_explained)

sim_data <- function()
{
  plyr::llply(N, function(n) {
    par <- get_par(n)
    dat <- gxefam_sim_dat(par, n = M, beta = c(betaG, 0, 0))

    tr <- as.matrix(gxefam_transform(par))
    assoc <- with(dat, matreg(Y = Y, X = G, transform = tr))

    mutate(assoc$tab, N = as.integer(n), 
      Neff = sigma2 / (2*p*(1-p)*se^2), gamma = Neff/N, sim = "rel_gen")
  }) %>% bind_rows
}

est_power <- function()
  tibble(N = N, beta = betaG, sigma2 = 1) %>%
    mutate(gamma = sapply(N, function(n) gxefam_trace_factor(get_par(n), mult = "K")),
      se = 1 / sqrt(2 * p * (1 - p) * N * gamma), ncp = (beta / se)^2,
      power = compute_power(ncp))

out <- run_sim(sim_data, est_power)
ggsave("rel-marginal.png", out$p, dpi = 150)

# simulation 1: unrel - marginal genetic effect (betaG)
sim_data <- function()
{
  plyr::llply(N, function(n) {
    par <- gxefam_par(num_fam = n / 5, sim_gen = "rbinom", 
      sim_phen = "no_rand", freq_g  = p)
    dat <- gxefam_sim_dat(par, n = M, beta = c(betaG, 0, 0))
    assoc <- with(dat, matreg(Y = Y, X = G)) 
    mutate(assoc$tab, N = as.integer(n), 
      Neff = sigma2 / (2*p*(1-p)*se^2), gamma = Neff/N, sim = "unrel")
  }) %>% bind_rows
}

est_power <- function()
  tibble(N = N, beta = betaG, sigma2 = 1) %>%
    mutate(se = 1 / sqrt(2 * p * (1 - p) * N), ncp = (beta / se)^2,
      power = compute_power(ncp), gamma = 1)

out <- run_sim(sim_data, est_power)
ggsave("unrel-marginal.png", out$p, dpi = 150)
