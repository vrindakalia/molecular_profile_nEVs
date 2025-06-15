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

clean_merge %>% write_tsv("HTG/clean_data/miRNA_clean_merged.txt")
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

brainrna_top30 <- data_summary %>% 
    filter(type == "eRNA") %>% 
    arrange(-mean.count) %>% 
    slice(1:30)

brainrna_top100 <- data_summary %>% 
    filter(type == "eRNA") %>% 
    arrange(-mean.count) %>% 
    slice(1:100)

brainrna_top100 %>% write_tsv("HTG/results/top100_brain_mirna_mean.txt")

# Correlation among the miRNA in common
mirna_brain_ev_ndev <- data_long %>% 
    filter(mirna %in% brainrna_top30$mirna)

mirna_brain_ev_ndev$type <- factor(mirna_brain_ev_ndev$type, levels = c("eRNA", "NDEVs ", "EVs", "Serum"))

mirna_brain_ev_ndev %>% 
    ggplot(aes(x = type, y = l2cpm)) +
    geom_boxplot(width = 0.4) +
    theme_bw() +
    facet_wrap(~mirna, scales = "free_y")

# Correlation of the 28 miRNA between brain and other compartments
# Brain and NDEVs
brain_ndev_plot <- data_long %>% 
    filter(mirna %in% brainrna_top30$mirna) %>% 
    filter(type %in% c("eRNA", "NDEVs ")) %>% 
    spread(key = "type", value = "l2cpm")

png("figures/20240125_individual_correlations_top30_brain_NDEV_miRNA.png", res = 300, units = "in", h = 6, w = 8.5)
brain_ndev_plot %>% 
    ggplot(aes(x = eRNA, y = `NDEVs `)) +
    geom_point() +
    geom_smooth(method = "lm") +
    facet_wrap(~mirna, scales = "free") +
    theme_bw() +
    labs(x = "Log2(CPM Tissue)", y = "Log2(CPM NDEV)")
dev.off()

brain_ndev_data <- data_long %>% 
    filter(mirna %in% brainrna_top30$mirna) %>% 
    filter(type %in% c("eRNA", "NDEVs ")) %>% 
    spread(key = "type", value = "l2cpm") %>% 
    mutate_at("mirna", as.character)

mirna_name <- as.character(brainrna_top30$mirna) 

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
brain_ev_plot <- data_long %>% 
    filter(mirna %in% brainrna_top30$mirna) %>% 
    filter(type %in% c("eRNA", "EVs")) %>% 
    spread(key = "type", value = "l2cpm")

png("figures/20240125_individual_correlations_top30_brain_TotalEVs_miRNA.png", res = 300, units = "in", h = 6, w = 8.5)
brain_ev_plot %>% 
    ggplot(aes(x = eRNA, y = EVs)) +
    geom_point() +
    geom_smooth(method = "lm") +
    facet_wrap(~mirna, scales = "free") +
    theme_bw() +
    labs(x = "Log2(CPM Tissue)", y = "Log2(CPM Total EVs)")
dev.off()

brain_ev_data <- data_long %>% 
    filter(mirna %in% brainrna_top30$mirna) %>% 
    filter(type %in% c("eRNA", "EVs")) %>% 
    spread(key = "type", value = "l2cpm") %>% 
    mutate_at("mirna", as.character)

mirna_name <- as.character(brainrna_top30$mirna)

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
brain_serum_plot <- data_long %>% 
    filter(mirna %in% brainrna_top30$mirna) %>% 
    filter(type %in% c("eRNA", "Serum")) %>% 
    spread(key = "type", value = "l2cpm")

png("figures/20240125_individual_correlations_top30_brain_serum_miRNA.png", res = 300, units = "in", h = 6, w = 8.5)
brain_serum_plot %>% 
    ggplot(aes(x = eRNA, y = Serum)) +
    geom_point() +
    geom_smooth(method = "lm") +
    facet_wrap(~mirna, scales = "free") +
    theme_bw() +
    labs(x = "Log2(CPM Tissue)", y = "Log2(CPM Serum)")
dev.off()

brain_serum_plot %>% 
    ggplot(aes(x = eRNA, y = Serum, color = mirna)) +
    geom_point() +
    geom_smooth(method = "lm") +
    theme_bw()

brain_serum_data <- data_long %>% 
    filter(mirna %in% brainrna_top30$mirna) %>% 
    filter(type %in% c("eRNA", "Serum")) %>% 
    spread(key = "type", value = "l2cpm") %>% 
    mutate_at("mirna", as.character)

mirna_name <- as.character(brainrna_top30$mirna)

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

png("figures/20240624_correlations_top30_brain_miRNA.png", res = 300, units = "in", h = 7, w = 4.5)
serum_ev_ndev_cor %>% 
    mutate(mir_name = gsub("mi_r", "miR", mirna)) %>% 
    mutate(mir_name = gsub("_", "-", mir_name)) %>% 
    ggplot(aes(x = matrix, y = reorder(mir_name, estimate), fill = estimate, label = round(estimate, 3))) +
    geom_tile(color = "white") +
    geom_text(size = 3) +
    theme_bw() +
    scale_fill_gradient2(high = "coral", low = "royalblue", mid = "white") +
    labs(x = "", y = "", fill = "coefficient") +
    theme(axis.text.x = element_text(size = 10)) 
dev.off()

mir_cor_coef <- serum_ev_ndev_cor %>% 
    mutate(mir_name = gsub("mi_r", "miR", mirna)) %>% 
    mutate(mir_name = gsub("_", "-", mir_name)) %>% 
    ggplot(aes(x = matrix, y = reorder(mir_name, estimate), fill = estimate, label = round(estimate, 3))) +
    geom_tile(color = "white") +
    geom_text(size = 3) +
    theme_bw() +
    scale_fill_gradient2(high = "coral", low = "royalblue", mid = "white") +
    labs(x = "", y = "", fill = "Coefficient") +
    theme(axis.text.x = element_text(size = 10)) 
