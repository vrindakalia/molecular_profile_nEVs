#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Update miRNA names for pathway analysis
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(tidyverse)

# Call in the file with miR information
mir <- read_csv("HTG/results/mirpath/miR_high_corr_nEVs_brain.csv")

mir_low <- read_csv("HTG/results/mirpath/miR_LOW_corr_nEVs_brain.csv")

mir <- mir %>% 
    mutate(use = gsub("mi_r_", "hsa-miR-", .$miR)) %>% 
    mutate(use_2 = gsub("_", "-", .$use)) %>% 
    select(-miR, -use)

mir_low <- mir_low %>% 
    mutate(use = gsub("mi_r_", "hsa-miR-", .$miRNA)) %>% 
    mutate(use_2 = gsub("_", "-", .$use)) %>% 
    select(-miRNA, -use)

mir %>% write_tsv("HTG/results/mirpath/miR_high_corr_nEVs_brain_updated.csv")
mir_low %>% write_tsv("HTG/results/mirpath/miR_LOW_corr_nEVs_brain_updated.csv")
