#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Compare NDEVs and serum
# Create time-estimate plots, colored by metabolite class
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(tidyverse)
library(VennDetail)

hilic_pos_ndev <- read_tsv("results/NDEV/hilic_pos_serum_features_NDEV.txt")
hilic_pos_serum <- read_tsv("results/NDEV/hilic_pos_serum_features_serum.txt")

ven <- venndetail(list(ndev = hilic_pos_ndev$feature, serum = hilic_pos_serum$feature))
plot(ven)
