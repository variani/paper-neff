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

## read h2 results
tab <- lapply(traits, function(t) {
  glue("out/clump/{t}.top.5000.clump.gz") %>% read_tsv %>%
    mutate(trait = t, order = seq(n()))
}) %>% bind_rows

## plot
ptab <- tab

up_p <- function(p) 
  p + 
    # geom_point(size = 3, shape = 1) + 
    geom_line(size = 1) +
    scale_color_manual(values = cols) +
    scale_y_log10() +
    labs(color = NULL, y = "-log10(p) at log scale")

p0 <- ggplot(ptab, aes(order, -log10(P), group = trait, color = trait)) %>% up_p

p <- plot_grid(p0, ncol = 1, align = "v")

ggsave("tmp.png", p, dpi = 100)#, width = 6, height = 8)

