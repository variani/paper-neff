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
ntop <- c(500, 1000, 2000)[-1]

## read h2 results
tab <- lapply(traits, function(t) {
  lapply(ntop, function(n) {
    bind_rows(
      glue("out/h2/{n}/{t}.tsv.gz") %>% read_tsv %>% mutate(method = "no_pcs"),
      glue("out/h2-pcs/{n}/{t}.tsv.gz") %>% read_tsv %>% mutate(method = "pcs_as_cov"))
  }) %>% bind_rows
}) %>% bind_rows

## plot
ptab <- mutate(tab, 
  M = as.factor(M), 
  h2_hat = gamma)

up_p <- function(p, i) 
  p + geom_point(size = 3, shape = 1) + 
    # geom_line(size = 1) +
    scale_color_manual(values = cols[i]) +
    labs(color = NULL) +
    theme(legend.position = "top")

plots <- lapply(seq_along(traits), function(i) {
  t <- traits[i]
  ggplot(filter(ptab, trait == t), aes(M, h2_hat, group = method, color = trait)) %>% up_p(i)
})

p <- plot_grid(plotlist = plots, align = "v")

ggsave("tmp.png", p, dpi = 100)#, width = 6, height = 8)

