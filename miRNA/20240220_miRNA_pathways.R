#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# miRNA pathways across matrices
# 2024 February 20
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(tidyverse)
source("code/20230509_miRNA_brain_high.R")

# Call in results from miRpath
brain <- read_csv("HTG/results/mirpath/brain_miRPath-v4-miRNA-centric_analysis_2024-02-20_union.csv")
ev <- read_csv("HTG/results/mirpath/EVs_miRPath-v4-miRNA-centric_analysis_2024-02-20_union.csv")
nev <- read_csv("HTG/results/mirpath/NDEV_miRPath-v4-miRNA-centric_analysis_2024-02-20_union.csv")
serum <- read_csv("HTG/results/mirpath/serum_miRPath-v4-miRNA-centric_analysis_2024-02-20_union.csv")

source("code/20240823_nEV_miRs_corr_plot.R")
brain_30 <- brain %>% 
    arrange(FDR) %>% 
    slice(1:30) %>% 
    mutate(matrix = "Tissue") %>%
    mutate(n = 100 - 9)

ev_brain <- ev %>% 
    filter(`Term Name` %in% brain_30$`Term Name`) %>% 
    mutate(matrix = "TotalEV") %>%
    mutate(n = 100 - 46)

nev_brain <- nev %>% 
    filter(`Term Name` %in% brain_30$`Term Name`) %>% 
    mutate(matrix = "nEV") %>%
    mutate(n = 100 - 57)

serum_brain <- serum %>% 
    filter(`Term Name` %in% brain_30$`Term Name`) %>% 
    mutate(matrix = "Serum") %>%
    mutate(n = 100 - 58)

all_brain <- rbind(brain_30, ev_brain, nev_brain, serum_brain)

all_brain$matrix <- factor(all_brain$matrix, levels = c("Tissue", "Serum", "TotalEV", "nEV"),
                           labels = c("Brain \ntissue", "Serum", "TotalEV", "nEV"))

neuro_paths <- c("Axon guidance",
                 "Amyotrophic lateral sclerosis", 
                 "Pathways of neurodegeneration - multiple diseases",
                 "Huntington disease",
                 "Parkinson disease",
                 "Alzheimer disease",
                 "Spinocerebellar ataxia",
                 "Neurotrophin signaling pathway")

path_table <- all_brain %>% 
    mutate(prop =  `miRNAs (n)`/n) 

path_table_wide <- path_table %>% 
    select(`Term Name`, matrix, prop) %>% 
    spread(key = "matrix", value = "prop")

path_table_wide %>% write_tsv("HTG/results/miRNA_pathways_proportions_matrix.txt")

png("/Users/vrindakalia/Documents/post-doc/Grants/Howie_Wellcome_Trust//20250310_pathways_top30_brain_miRNA_neuro.png", res = 300, units = "in", h = 7, w = 6.5)
all_brain %>% 
    filter(`Term Name` %in% neuro_paths) %>% 
    mutate(prop =  `miRNAs (n)`/n) %>% 
    ggplot(aes(x = matrix, y = reorder(`Term Name`, prop))) +
    geom_point(aes(size = -log10(FDR), color = `miRNAs (n)`/n)) +
    #geom_text(aes(label = round(`miRNAs (n)`/n, 2)), size = 2.5, color = "white") +
    scale_color_gradient(low = "grey96", high = "midnightblue") +
    theme_bw() +
    labs(x = "", y = "", color = "Proportion\n of miRNA") +
    theme(axis.text.x = element_text(size = 10))
dev.off()               

mir_paths <- all_brain %>% 
    mutate(prop =  `miRNAs (n)`/n) %>% 
    ggplot(aes(x = matrix, y = reorder(`Term Name`, prop))) +
    geom_point(aes(size = -log10(FDR), color = `miRNAs (n)`/n)) +
    #geom_text(aes(label = round(`miRNAs (n)`/n, 2)), size = 2.5, color = "white") +
    scale_color_gradient(low = "grey96", high = "midnightblue") +
    theme_bw() +
    labs(x = "", y = "", color = "Proportion\n of miRNA", size = "-Log10(q)") +
    theme(axis.text.x = element_text(size = 10))

png("figures/figure2_miRs_grey.png", res = 300, units = "in",
    h = 7, w = 11)
cowplot::plot_grid(mir_paths, mir_cor_coef, nrow = 1, labels = c("A", "B"),
                   rel_widths = c(1, 0.6))
dev.off()

png("figures/figure2_updated_nEVs.png", res = 300, units = "in",
    h = 8, w = 16)
cowplot::plot_grid(mir_paths, mir_cor_coef, nEV_corrs, nrow = 1, labels = c("A", "B", "C"),
                   rel_widths = c(1, 0.65, 0.9))
dev.off()

test_brain <- all_brain %>% 
    mutate(prop =  `miRNAs (n)`/n)

aov.test <- aov(prop~matrix, data = test_brain)
aov.test %>% 
    broom::tidy() %>% 
    write_tsv("HTG/results/miR_pathways_prop_anova.txt")
TukeyHSD(aov.test) %>% 
    broom::tidy() %>% 
    write_tsv("HTG/results/miR_pathways_prop_tukey.txt")
