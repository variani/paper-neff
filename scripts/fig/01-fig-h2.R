### inc
library(tidyverse)
library(glue)

library(ggplot2)
library(cowplot)
# theme_set(theme_cowplot())
theme_set(theme_minimal_hgrid())

library(wesanderson)
cols <- wes_palette("Darjeeling1") %>% as.character
cols <- c("grey50", cols)

## arguments
traits <- c("bmi", "weight", "waist", "hip", "height", "whr")
ntop <- c(500, 1000, 2000)

## read h2 results
tab <- lapply(traits, function(t) {
  lapply(ntop, function(n) {
    glue("out/h2/{n}/{t}.tsv.gz") %>% read_tsv
  }) %>% bind_rows
}) %>% bind_rows

## plot
ptab <- mutate(tab, 
  M = as.factor(M), 
  h2_hat = gamma,
  Neff = N * mult)

up_p <- function(p) 
  p + geom_point(size = 3, shape = 1) + geom_line(size = 1) +
    scale_color_manual(values = cols) +
    labs(color = NULL)

p1 <- ggplot(ptab, aes(M, h2_hat, group = trait, color = trait)) %>% up_p

p2 <- ggplot(ptab, aes(M, mult, group = trait, color = trait)) %>% up_p

p3 <- ggplot(ptab, aes(M, Neff, group = trait, color = trait)) %>% up_p
p3 <- p3 + geom_hline(aes(yintercept = N), linetype = 3)

p <- plot_grid(p1, p2, p3, ncol = 1, align = "v")

ggsave("tmp.png", p, dpi = 100)#, width = 6, height = 8)

