#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# mean levels across compartments
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(tidyverse)

# exposomics
chem <- read_tsv("exposomics/clean_data/in_four_replaced_atl3.txt")
chem_key <- read_tsv("exposomics/clean_data/chemical_key.txt")
chem_ids <- read_tsv("exposomics/clean_data/gc_ids.txt")

chem_plot <- chem %>% 
    merge(chem_ids, by = "ID") %>% 
    separate(ID, into = c("matrix", "sample"), sep = "_") %>% 
    gather(key = "index", value = "concentration", -matrix, -sample, -`Sample ID`) %>% 
    merge(., chem_key, by = "index") %>% 
    filter(matrix != "Depleted") %>% 
    mutate(conc.ppb = ifelse(matrix == "Tissue", concentration*0.001, concentration))

chem_plot$matrix <- factor(chem_plot$matrix, levels = c("Tissue", "NDEV", "TotalEV", "Serum"),
                           labels = c("Brain \ntissue", "nEV", "TotalEV", "Serum"))

png("exposomics/results/20240724_boxplot_matrix_levels.png", res = 300,
    units = "in", h = 5, w = 13)
chem_plot %>% 
    ggplot(aes(x = matrix, y = log2(conc.ppb))) +
    geom_boxplot() +
    theme_bw() +
    facet_wrap(~Chemical, nrow = 3, scales = "free_y") +
    labs(y = "Log2(Concentration [ppb])", x = "")
dev.off()

png("exposomics/results/20240724_lineplot_samples_matrix_levels.png", res = 300,
    units = "in", h = 5, w = 13)
chem_plot %>% 
    ggplot(aes(x = matrix, y = log2(conc.ppb), color = factor(`Sample ID`))) +
    geom_point() +
    geom_line(aes(group = `Sample ID`)) +
    theme_bw() + 
    facet_wrap(~Chemical, scales = "free_y", nrow = 3) +
    theme(legend.position = "bottom") +
    labs(color = "ID", x = "", y = "Log2(Concentration [ppb])")
dev.off()

chem_plot %>% 
    group_by(Chemical, matrix) %>% 
    summarise(mean.conc = mean(concentration)) %>% 
    spread(key = "matrix", value = "mean.conc") %>% 
    write_tsv("exposomics/results/20240724_mean_concentration_chemicals_matrix.txt")

# Run an ANOVA to look for difference in levels of chemicals by matrix
chem_test <- chem_plot %>% 
    mutate_at("conc.ppb", log2)

chem_test$matrix <- factor(chem_plot$matrix, levels = c("Brain \ntissue", "nEV", "TotalEV", "Serum"),
                           labels = c("Tissue", "nEV", "TotalEV", "Serum"))

unique_chem <- chem_test %>% 
    distinct(Chemical)

results <- vector(mode = "list", length = dim(unique_chem)[1])
names(results) <- unique_chem$Chemical

for(i in unique_chem$Chemical){
    sub <- chem_test %>% 
        filter(Chemical == i)
    
    results[[i]] <- aov(conc.ppb ~ matrix, data = sub)
    
}

summary(results[[2]])

results_tidy <- results %>% 
    map_df(broom::tidy) %>% 
    filter(term == "matrix") %>% 
    mutate(name = unique_chem$Chemical) %>% 
    mutate(fdr = p.adjust(p.value, method = "BH"))

posthoc_results <- results %>% 
    map(~TukeyHSD(.x))

posthoc_results_tidy <- posthoc_results %>% 
    map_df(broom::tidy)

# Contrast: Brain tissue - others
posthoc_results_nEV_brain <- posthoc_results_tidy %>% 
    filter(contrast == "nEV-Tissue") %>% 
    mutate(name = unique_chem$Chemical) %>% 
    mutate(qvalue = p.adjust(adj.p.value, method = "BH"))

posthoc_results_TotalEV_brain <- posthoc_results_tidy %>% 
    filter(contrast == "TotalEV-Tissue") %>% 
    mutate(name = unique_chem$Chemical) %>% 
    mutate(qvalue = p.adjust(adj.p.value, method = "BH"))

posthoc_results_serum_brain <- posthoc_results_tidy %>% 
    filter(contrast == "Serum-Tissue") %>% 
    mutate(name = unique_chem$Chemical) %>% 
    mutate(qvalue = p.adjust(adj.p.value, method = "BH"))

# Contrast: TotalEV - others
posthoc_results_serum_totalEV <- posthoc_results_tidy %>% 
    filter(contrast == "Serum-TotalEV") %>% 
    mutate(name = unique_chem$Chemical) %>% 
    mutate(qvalue = p.adjust(adj.p.value, method = "BH"))

posthoc_results_nEV_totalEV <- posthoc_results_tidy %>% 
    filter(contrast == "TotalEV-nEV") %>% 
    mutate(name = unique_chem$Chemical) %>% 
    mutate(qvalue = p.adjust(adj.p.value, method = "BH"))

# Contrast: nEV - serum
posthoc_results_serum_nEV <- posthoc_results_tidy %>% 
    filter(contrast == "Serum-nEV") %>% 
    mutate(name = unique_chem$Chemical) %>% 
    mutate(qvalue = p.adjust(adj.p.value, method = "BH"))

# Combine the posthoc results
all_posthoc <- rbind(posthoc_results_nEV_brain, posthoc_results_TotalEV_brain, posthoc_results_serum_brain,
                     posthoc_results_nEV_totalEV, posthoc_results_serum_nEV, posthoc_results_serum_totalEV)

all_posthoc %>% write_tsv("exposomics/results/20240725_chemicals_matrix_ANOVA_posthoc.txt")
