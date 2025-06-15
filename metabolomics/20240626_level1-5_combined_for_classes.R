#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Combining annotations for main class box plots
# 2024.06.26
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(tidyverse)

ndev_all <- read_tsv("results/NDEV/ndev_hmdbids.txt") %>% 
    unite("mz_time", c(mz, time), sep = "_", remove = FALSE)
ev_all <- read_tsv("results/NDEV/TotalEV_hmdbids.txt") %>% 
    unite("mz_time", c(mz, time), sep = "_", remove = FALSE)
serum_all <- read_tsv("results/NDEV/serum_hmdbids.txt") %>% 
    unite("mz_time", c(mz, time), sep = "_", remove = FALSE)

serum_not_in_ndev <- serum_all %>% 
    filter(!mz_time %in% ndev_all$mz_time)

tev_not_in_ndev <- ev_all %>% 
    filter(!mz_time %in% ndev_all$mz_time)

tev_not_in_serum<- ev_all %>% 
    filter(!mz_time %in% serum_all$mz_time)

all_unique <- rbind(ndev_all, serum_not_in_ndev, tev_not_in_ndev, tev_not_in_serum) %>% 
    distinct(mz_time, .keep_all = TRUE)

all_unique %>% write_tsv("results/NDEV/all_matrix_hmdbids.txt")
