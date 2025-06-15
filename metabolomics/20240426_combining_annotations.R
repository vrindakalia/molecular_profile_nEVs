#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Combine all annotation results into one file
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(tidyverse)

# Call in the files
# NDEV ----
c18neg_ndev <- read_tsv("results/NDEV/level1-5_c18_neg_features_NDEV.txt") %>% 
    mutate(col = "C18pos")
c18pos_ndev <- read_tsv("results/NDEV/level1-5_c18_pos_features_NDEV.txt") %>% 
    mutate(col = "C18neg")
hilpos_ndev <- read_tsv("results/NDEV/level1-5_hilic_pos_features_NDEV.txt") %>% 
    mutate(col = "Hilpos")
hilneg_ndev <- read_tsv("results/NDEV/level1-5_hilic_neg_features_NDEV.txt") %>% 
    mutate(col = "Hilneg")

ndev_all <- rbind(c18neg_ndev, c18pos_ndev, hilpos_ndev, hilneg_ndev) %>% 
    separate(Name, into = c("use_name", "other_1", "other_2"), sep = ";", remove = FALSE) %>% 
    mutate(pubchem_CID = case_when(grepl("Pub", ID) ~ gsub("PubChem_", "", ID),
                                   !grepl("Pub", ID) ~ ""))

# EV ----
c18neg_ev <- read_tsv("results/NDEV/level1-5_c18_neg_features_totalEV.txt") %>% 
    mutate(col = "C18pos")
c18pos_ev <- read_tsv("results/NDEV/level1-5_c18_pos_features_totalEV.txt") %>% 
    mutate(col = "C18neg")
hilpos_ev <- read_tsv("results/NDEV/level1-5_hilic_pos_features_totalEV.txt") %>% 
    mutate(col = "Hilpos")
hilneg_ev <- read_tsv("results/NDEV/level1-5_hilic_neg_features_totalEV.txt") %>% 
    mutate(col = "Hilneg")

ev_all <- rbind(c18neg_ev, c18pos_ev, hilpos_ev, hilneg_ev)%>% 
    separate(Name, into = c("use_name", "other_1", "other_2"), sep = ";", remove = FALSE) %>% 
    mutate(pubchem_CID = case_when(grepl("Pub", ID) ~ gsub("PubChem_", "", ID),
                                   !grepl("Pub", ID) ~ ""))

# Serum ----
c18neg_serum <- read_tsv("results/NDEV/level1-5_c18_neg_features_serum.txt") %>% 
    mutate(col = "C18pos")
c18pos_serum <- read_tsv("results/NDEV/level1-5_c18_pos_features_serum.txt") %>% 
    mutate(col = "C18neg")
hilpos_serum <- read_tsv("results/NDEV/level1-5_hilic_pos_features_serum.txt") %>% 
    mutate(col = "Hilpos")
hilneg_serum <- read_tsv("results/NDEV/level1-5_hilic_neg_features_serum.txt") %>% 
    mutate(col = "Hilneg")

serum_all <- rbind(c18neg_serum, c18pos_serum, hilpos_serum, hilneg_serum)%>% 
    separate(Name, into = c("use_name", "other_1", "other_2"), sep = ";", remove = FALSE) %>% 
    mutate(pubchem_CID = case_when(grepl("Pub", ID) ~ gsub("PubChem_", "", ID),
                                   !grepl("Pub", ID) ~ ""))

# Save the results
ndev_all %>% write_tsv("results/NDEV/ndev_all_col_annotations.txt")
ev_all %>% write_tsv("results/NDEV/TotalEV_all_col_annotations.txt")
serum_all %>% write_tsv("results/NDEV/Serum_all_col_annotations.txt")

# Merge HMDBIDs for level 1 
ndev_name_map <- read_csv("results/NDEV/ndev_level1_name_map.csv") %>% 
    mutate_at("Query", as.character) %>% 
    dplyr::rename(pubchem_CID = Query)
    
ev_name_map <- read_csv("results/NDEV/TotalEV_level1_name_map.csv") %>% 
    mutate_at("Query", as.character) %>% 
    dplyr::rename(pubchem_CID = Query)

serum_name_map <- read_csv("results/NDEV/serum_level1_name_map.csv") %>% 
    mutate_at("Query", as.character) %>% 
    dplyr::rename(pubchem_CID = Query)

ndev_all_map <- left_join(ndev_all, ndev_name_map, by = "pubchem_CID") %>% 
    mutate(hmdb_id_new = case_when(confidence == "Level 1" ~ HMDB,
                                   confidence != "Level 1" ~ ID)) %>% 
    dplyr::select(mz, time, mean.int, Name, use_name, confidence, ID, pubchem_CID, hmdb_id_new)

ev_all_map <- left_join(ev_all, ev_name_map, by = "pubchem_CID") %>% 
    mutate(hmdb_id_new = case_when(confidence == "Level 1" ~ HMDB,
                                   confidence != "Level 1" ~ ID)) %>% 
    dplyr::select(mz, time, mean.int, Name, use_name, confidence, ID, pubchem_CID, hmdb_id_new)

serum_all_map <- left_join(serum_all, serum_name_map, by = "pubchem_CID") %>% 
    mutate(hmdb_id_new = case_when(confidence == "Level 1" ~ HMDB,
                                   confidence != "Level 1" ~ ID)) %>% 
    dplyr::select(mz, time, mean.int, Name, use_name, confidence, ID, pubchem_CID, hmdb_id_new)

hist(log2(ndev_all_map$mean.int))
hist(log2(ev_all_map$mean.int))
hist(log2(serum_all_map$mean.int))

ndev_all_map %>% write_tsv("results/NDEV/ndev_hmdbids.txt")
ev_all_map %>% write_tsv("results/NDEV/TotalEV_hmdbids.txt")
serum_all_map %>% write_tsv("results/NDEV/serum_hmdbids.txt")

ndev_all_map %>% 
    filter(mean.int > median(mean.int)) %>% 
    write_tsv("results/NDEV/gt_median_ndev_hmdbids.txt") 
ev_all_map %>% 
    filter(mean.int > median(mean.int)) %>% 
    write_tsv("results/NDEV/gt_median_TotalEV_hmdbids.txt")
serum_all_map %>% 
    filter(mean.int > median(mean.int)) %>% 
    write_tsv("results/NDEV/gt_median_serum_hmdbids.txt")

ndev_all_map %>% 
    filter(mean.int > quantile(mean.int,0.75)) %>% 
    write_tsv("results/NDEV/gt_q3_ndev_hmdbids.txt") 
ev_all_map %>% 
    filter(mean.int > quantile(mean.int,0.75)) %>% 
    write_tsv("results/NDEV/gt_q3_TotalEV_hmdbids.txt")
serum_all_map %>% 
    filter(mean.int > quantile(mean.int,0.75)) %>% 
    write_tsv("results/NDEV/gt_q3_serum_hmdbids.txt")

ndev_all_map %>% 
    filter(mean.int > quantile(mean.int,0.9)) %>% 
    write_tsv("results/NDEV/gt_90perc_ndev_hmdbids.txt") 
ev_all_map %>% 
    filter(mean.int > quantile(mean.int,0.9)) %>% 
    write_tsv("results/NDEV/gt_90perc_TotalEV_hmdbids.txt")
serum_all_map %>% 
    filter(mean.int > quantile(mean.int,0.9)) %>% 
    write_tsv("results/NDEV/gt_90perc_serum_hmdbids.txt")
