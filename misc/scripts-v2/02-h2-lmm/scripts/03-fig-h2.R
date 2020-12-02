### inc
library(tidyverse)
library(glue)

library(ggplot2)
library(cowplot)
theme_set(theme_minimal())

## arguments
traits <- c("bmi", "weight", "waist", "hip", "height", "whr",
  "pdw", "rbc", "rbcdw")
# n_tops <- c(100, 500, 1000, 2000, 3000, 4000)
n_tops <- c(500, 2000, 4000)

## read h2 results
tab <- lapply(traits, function(t) {
  lapply(n_tops, function(n) {
    f <- glue("output/h2/{t}.{n}.tsv.gz")
    if(file.exists(f)) {
      read_tsv(f) %>% mutate(trait = t)
    } else {
      tibble()
    }
  }) %>% bind_rows
}) %>% bind_rows

## plot
ptab <- mutate(tab, 
  M = as.factor(M), 
  gamma = trace_factor,
  gamma_sc = trace_factor / s2,
  Neff = N * gamma_sc)

p1 <- ggplot(ptab, aes(M, h2, group = trait, color = trait)) + 
  geom_point() + geom_line()

p2 <- ggplot(ptab, aes(M, gamma_sc, group = trait, color = trait)) + 
  geom_point() + geom_line() 

p3 <- ggplot(ptab, aes(M, Neff, group = trait, color = trait)) + 
  geom_point() + geom_line() + geom_hline(aes(yintercept = N), linetype = 3)

p <- plot_grid(p1, p2, p3, ncol = 1)

ggsave("tmp.png", p, dpi = 100, width = 6, height = 8)

