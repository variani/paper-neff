library(tidyverse)

library(BEDMatrix) # read bed from file & avoid loading into RAM

## phen_synca
file_bed <- "dat/bed_top/body_BMIz_orth.top.5000.bed" 
file_phen <- "dat/ukb.jass.tab.gz"

## bed
bed <- BEDMatrix(file_bed)
ids_bed <- rownames(bed)

# fix id names such "ID_ID"
ids_bed <- ids_bed %>% strsplit("_") %>% sapply(function(x) tail(x, 1))

## phen: load phen, align ids with bed, extract y
phen <- read_tsv(file_phen) %>% 
  select(-id2, -famid) %>% mutate(id = as.character(id))
ids_phen <- phen$id

ids_common <- intersect(ids_bed, ids_phen)

phen_sync <- tibble(id = ids_common) %>% left_join(phen) %>%
  select(id, everything())

## process columns of phen_sync
phen_sync <- rename(phen_sync, c(
  height = "body_HEIGHTz_orth",
  hip = "body_HIP_orth",
  bmi = "body_BMIz_orth",
  waist = "body_WAIST_orth",
  weight = "body_WEIGHT_orth",
  whr = "body_WHR_orth"))

# names of traits
traits <- c("height", "hip", "bmi", "waist", "weight", "whr")

# final transformations: impute and scale
impute <- function(x) { x[is.na(x)] <- mean(x, na.rm = TRUE); x }
inormal <- function(x) qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))

phen_sync <- mutate_at(phen_sync, traits, impute) 
phen_sync <- mutate_at(phen_sync, traits, inormal) 

## pcs
pcs <- list.files("dat/", "pcs", full = TRUE) %>%
  lapply(read_tsv) %>% bind_cols
names(pcs)[1] <- "id"
pcs <- select(pcs, id, starts_with("PC")) %>%
  mutate(id = as.character(id))

stopifnot(all(ids_common %in% pcs$id))
pcs_sync <- tibble(id = ids_common) %>% left_join(pcs) %>%
  select(id, everything())
stopifnot(all(!is.na(pcs_sync$PC1)))

## save
write_lines(phen_sync$id, "sync.ids")
write_tsv(phen_sync, "phen.sync.tsv.gz")
write_tsv(pcs_sync, "pcs.sync.tsv.gz")

