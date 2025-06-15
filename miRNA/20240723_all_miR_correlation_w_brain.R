#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Exploring HTG miRNA data
# 04/27/2023
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Call in data 
library(tidyverse)
library(openxlsx)
key <- read.xlsx("HTG/VLP00932 Sample ID Key.xlsx") %>% 
    rename(sample_id = 4)

key$samples <- paste0("x", key$sample_id)
key$samples <- tolower(key$samples)
key$samples <- gsub("\\-", "\\_", key$samples)

tmm_cpm <- read_csv("HTG/clean_data/log2cpm_TMM_normalized_HW.csv")

clean <- tmm_cpm %>% 
    t() %>% 
    as.data.frame() %>% 
    janitor::row_to_names(1) %>% 
    janitor::clean_names() %>% 
    rownames_to_column(var = "samples") %>% 
    filter(!grepl("mtc", samples))

clean_merge <- merge(key, clean, by = "samples")

# Objectives 
# 1. miRNA measured in brain (eRNA) are positively correlated with NDEV miRNA
# 2. How are these different in their correlation when compared to miRNA in total EVs or serum?

data_long <- clean_merge %>%
    select(ID = Alternate.Customer.ID, type = Sample.Description, everything()) %>% 
    select(-Accession.ID, -HTG.ID, -samples, -sample_id) %>% 
    gather(key = "mirna", value = "l2cpm", -ID, -type) %>% 
    mutate_at("l2cpm", as.character) %>% 
    mutate_at("l2cpm", as.numeric) 

# Most abundant miRNA in the four sample types
data_summary <- clean_merge %>%
    select(ID = Alternate.Customer.ID, type = Sample.Description, everything()) %>% 
    select(-Accession.ID, -HTG.ID, -samples, -sample_id) %>% 
    gather(key = "mirna", value = "l2cpm", -ID, -type) %>% 
    mutate_at("l2cpm", as.character) %>% 
    mutate_at("l2cpm", as.numeric) %>% 
    group_by(mirna, type) %>% 
    summarise(mean.count = mean(l2cpm, na.rm = TRUE)) %>% 
    ungroup()


# Correlation among the miRNA in common
mirna_brain_ev_ndev <- data_long 

mirna_brain_ev_ndev$type <- factor(mirna_brain_ev_ndev$type, levels = c("eRNA", "NDEVs ", "EVs", "Serum"))

# Correlation of the 28 miRNA between brain and other compartments
# Brain and NDEVs
brain_ndev_plot <- data_long %>% 
    #filter(mirna %in% brainrna_top30$mirna) %>% 
    filter(type %in% c("eRNA", "NDEVs ")) %>% 
    spread(key = "type", value = "l2cpm")

brain_ndev_data <- data_long %>% 
    filter(type %in% c("eRNA", "NDEVs ")) %>% 
    spread(key = "type", value = "l2cpm") %>% 
    mutate_at("mirna", as.character)

all_mirna <- brain_ndev_data %>% 
    distinct(mirna)

mirna_name <- as.character(all_mirna$mirna) 

ndev_brain_cor <- list()

for(i in 1:length(mirna_name)){
    brain_ndev_subset <- brain_ndev_data %>% 
        filter(mirna %in% mirna_name[i])
    
    ndev_brain_cor[[i]] <- cor.test(brain_ndev_subset$eRNA, brain_ndev_subset$`NDEVs `, method = "spearman")
}

ndev_brain_cor_results <- ndev_brain_cor %>% 
    map_df(broom::tidy) %>% 
    mutate(mirna = mirna_name) %>% 
    mutate(matrix = "nEVs")

# Brain and EVs
brain_ev_data <- data_long %>% 
    filter(type %in% c("eRNA", "EVs")) %>% 
    spread(key = "type", value = "l2cpm") %>% 
    mutate_at("mirna", as.character)

ev_brain_cor <- list()

for(i in 1:length(mirna_name)){
    brain_ev_subset <- brain_ev_data %>% 
        filter(mirna %in% mirna_name[i])
    
    ev_brain_cor[[i]] <- cor.test(brain_ev_subset$eRNA, brain_ev_subset$EVs, method = "spearman")
}

ev_brain_cor_results <- ev_brain_cor %>% 
    map_df(broom::tidy) %>% 
    mutate(mirna = mirna_name) %>% 
    mutate(matrix = "Total EVs")

# Brain and Serum 
brain_serum_data <- data_long %>% 
    filter(type %in% c("eRNA", "Serum")) %>% 
    spread(key = "type", value = "l2cpm") %>% 
    mutate_at("mirna", as.character)

serum_brain_cor <- list()

for(i in 1:length(mirna_name)){
    brain_serum_subset <- brain_serum_data %>% 
        filter(mirna %in% mirna_name[i])
    
    serum_brain_cor[[i]] <- cor.test(brain_serum_subset$eRNA, brain_serum_subset$Serum, method = "spearman")
}

serum_brain_cor_results <- serum_brain_cor %>% 
    map_df(broom::tidy) %>% 
    mutate(mirna = mirna_name) %>% 
    mutate(matrix = "Serum")

serum_ev_ndev_cor <- rbind(ndev_brain_cor_results, ev_brain_cor_results, serum_brain_cor_results)

serum_ev_ndev_cor$matrix <- factor(serum_ev_ndev_cor$matrix, levels = c("Serum", "Total EVs", "nEVs"),
                                   labels = c("Serum", "TotalEV", "nEV"))

# Save file to disc
serum_ev_ndev_cor %>% write_tsv("HTG/results/20240723_all_miR_correlations_w_brain.txt")

ndev_brain_cor_results_small <- ndev_brain_cor_results %>% 
    select(mirna, estimate, p.value) %>% 
    rename(nEV.estimate = estimate, nEV.p = p.value)
ev_brain_cor_results_small <- ev_brain_cor_results %>% 
    select(mirna, estimate, p.value) %>% 
    rename(TotalEV.estimate = estimate, TotalEV.p = p.value)
serum_brain_cor_results_small <- serum_brain_cor_results %>% 
    select(mirna, estimate, p.value) %>% 
    rename(serum.estimate = estimate, serum.p = p.value)

cor_all_wide <- merge(ndev_brain_cor_results_small, ev_brain_cor_results_small, by = "mirna") %>% 
    merge(., serum_brain_cor_results_small, by = "mirna")
cor_all_wide %>% write_tsv("HTG/results/20240724_all_miR_correlations_w_brain_wide.txt")
