#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal library annotations
# 2024, Februrary 20
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(tidyverse)
library(openxlsx)
library(xMSanalyzer)

# C18 neg
# Call in the annotation files
c18_neg_lib <- read_tsv("results/C18neg_annotations_level1-5.txt") %>% 
    dplyr::rename(feature = mz_time)

# Call in the feature list of interest
ndev_C18neg <- read_tsv("results/NDEV/c18_neg_serum_features_NDEV.txt")

c18neg_ann <- left_join(ndev_C18neg, c18_neg_lib, by = "feature")

c18neg_ann %>% write_tsv("results/NDEV/level1_c18_neg_features_NDEV.txt")

# C18 pos 
# Call in the annotation files
c18_pos_lib <- read_tsv("results/C18pos_annotations_level1-5.txt") %>% 
    dplyr::rename(feature = mz_time)

# Call in the feature list of interest
ndev_C18pos <- read_tsv("results/NDEV/c18_pos_serum_features_NDEV.txt")

c18pos_ann <- left_join(ndev_C18pos, c18_pos_lib, by = "feature")

c18pos_ann %>% write_tsv("results/NDEV/level1_c18_pos_features_NDEV.txt")

# HILIC pos

# HILIC neg