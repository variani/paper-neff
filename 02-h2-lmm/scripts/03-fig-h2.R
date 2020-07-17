### inc
library(tidyverse)
library(glue)

library(ggplot2)
library(cowplot)
theme_set(theme_minimal())

## arguments
traits <- c("bmi", "weight", "waist", "hip", "height", "whr")
n_tops <- c(100, 500, 1000, 2000)

## read h2 results
tab <- lapply(traits, function(t) {
  lapply(n_tops, function(n) {
    f <- glue("output/h2/{t}.{n}.tsv.gz")
    read_tsv(f)
  }) %>% bind_rows
}) %>% bind_rows

## plot
ptab <- mutate(tab, 
  M = as.factor(M), 
  gamma = trace_factor,
  Neff = N * gamma)

p1 <- ggplot(ptab, aes(M, h2, group = trait, color = trait)) + 
  geom_point() + geom_line()

p2 <- ggplot(ptab, aes(M, gamma, group = trait, color = trait)) + 
  geom_point() + geom_line() 

p3 <- ggplot(ptab, aes(M, Neff, group = trait, color = trait)) + 
  geom_point() + geom_line() + geom_hline(aes(yintercept = N), linetype = 3)

p <- plot_grid(p1, p2, p3, ncol = 1)

ggsave("ukb-h2.png", p, dpi = 100)

