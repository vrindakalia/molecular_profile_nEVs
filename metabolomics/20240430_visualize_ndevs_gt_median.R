#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Call in the results from MetaboAnalyst
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(tidyverse)
library(RColorBrewer)
# Function for plotting colors side-by-side
pal <- function(col, border = "light gray", ...){
    n <- length(col)
    plot(0, 0, type="n", xlim = c(0, 1), ylim = c(0, 1),
         axes = FALSE, xlab = "", ylab = "", ...)
    rect(0:(n-1)/n, 0, 1:n/n, 1, col = col, border = border)
}
# Super class ----
ndev_super <- read_csv("results/NDEV/enrichment_gt_median/ndev_super_class/msea_ora_result.csv") %>% 
    mutate(matrix = "NDEV") %>% 
    dplyr::rename(Name = 1)
TotalEV_super <- read_csv("results/NDEV/enrichment_gt_median/Totalev_super_class/msea_ora_result.csv") %>% 
    mutate(matrix = "TotalEV") %>% 
    dplyr::rename(Name = 1)
serum_super <- read_csv("results/NDEV/enrichment_gt_median/serum_super_class/msea_ora_result.csv") %>% 
    mutate(matrix = "Serum") %>% 
    dplyr::rename(Name = 1)

super_class_all <- rbind(ndev_super, TotalEV_super, serum_super) 

super_class_all %>% 
    group_by(matrix) %>% 
    mutate(total = sum(hits)) %>% 
    ungroup() %>% 
    mutate(prop = hits/total*100) %>% 
    ggplot(aes(x = prop)) +
    geom_histogram(binwidth = 0.5)

tol21rainbow= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")
pal(tol21rainbow[1:16])

names <- super_class_all %>% 
    group_by(matrix) %>% 
    mutate(total = sum(hits)) %>% 
    ungroup() %>% 
    mutate(prop = hits/total*100) %>% 
    mutate(class = case_when(prop > 5 ~ Name,
                             prop <= 5 ~ "z_Other")) %>% 
    distinct(class) %>% 
    mutate(color = tol21rainbow[1:7]) %>% 
    column_to_rownames(var = "class")

png("results/figures/super_classes_totalEVs_blood_atleast5perc.png", res = 300, units = "in", height = 8, width = 6)
super_class_all %>% 
    group_by(matrix) %>% 
    mutate(total = sum(hits)) %>% 
    ungroup() %>% 
    mutate(prop = hits/total*100) %>% 
    mutate(class = case_when(prop > 5 ~ Name,
                             prop <= 5 ~ "z_Other")) %>% 
    group_by(matrix, class) %>% 
    mutate(hits.e = sum(hits)) %>% 
    ungroup() %>% 
    group_by(matrix) %>% 
    distinct(class, .keep_all = TRUE) %>% 
    ungroup() %>% 
    ggplot() +
    geom_col(aes(x = matrix, y = hits.e, fill = class), color = "grey20", width = 0.6, position = "fill") +
    scale_fill_manual(values = names$color) +
    theme_bw() +
    theme(legend.text = element_text(size=7),
          legend.key.size = unit(0.5, 'cm')) +
    guides(fill = guide_legend(ncol = 1)) +
    labs(x = "", y = "Proportion of hits", fill = "Class")
dev.off()

names_other <- super_class_all %>% 
    group_by(matrix) %>% 
    mutate(total = sum(hits)) %>% 
    ungroup() %>% 
    mutate(prop = hits/total*100) %>% 
    mutate(class = case_when(prop > 5 ~ Name,
                             prop <= 5 ~ "z_Other")) %>% 
    filter(class == "z_Other") %>% 
    distinct(Name) %>% 
    mutate(color = tol21rainbow[1:13]) %>% 
    column_to_rownames(var = "Name")

super_class_all %>% 
    group_by(matrix) %>% 
    mutate(total = sum(hits)) %>% 
    ungroup() %>% 
    mutate(prop = hits/total*100) %>% 
    mutate(class = case_when(prop > 5 ~ Name,
                             prop <= 5 ~ "z_Other")) %>% 
    filter(class == "z_Other") %>% 
    ggplot() +
    geom_bar(aes(x = matrix, y = prop, fill = Name), color = "grey20", width = 0.6, position = "stack", stat = "identity") +
    scale_fill_manual(values = names_other$color) +
    theme_bw() +
    theme(legend.text = element_text(size=7),
          legend.key.size = unit(0.5, 'cm')) +
    guides(fill = guide_legend(ncol = 1)) +
    labs(x = "", y = "Proportion of hits", fill = "Class")




# Main class ----
ndev_main <- read_csv("results/NDEV/enrichment_gt_median/ndev_main_class/msea_ora_result.csv") %>% 
    mutate(matrix = "NDEV") %>% 
    dplyr::rename(Name = 1)
Totalev_main <- read_csv("results/NDEV/enrichment_gt_median/Totalev_main_class/msea_ora_result.csv") %>% 
    mutate(matrix = "TotalEV") %>% 
    dplyr::rename(Name = 1)
serum_main <- read_csv("results/NDEV/enrichment_gt_median/serum_main_class/msea_ora_result.csv") %>% 
    mutate(matrix = "Serum") %>% 
    dplyr::rename(Name = 1)

main_class_all <- rbind(ndev_main, Totalev_main, serum_main)

main_class_all %>% 
    group_by(matrix) %>% 
    mutate(total = sum(hits)) %>% 
    ungroup() %>% 
    mutate(prop = hits/total*100) %>% 
    ggplot(aes(x = prop)) +
    geom_histogram(binwidth = 0.5)

tol21rainbow= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")
pal(tol21rainbow)

names_main <- main_class_all %>% 
    group_by(matrix) %>% 
    mutate(total = sum(hits)) %>% 
    ungroup() %>% 
    mutate(prop = hits/total*100) %>% 
    mutate(class = case_when(prop > 2.5 ~ Name,
                             prop <= 2.5 ~ "z_Other")) %>% 
    distinct(class) %>% 
    mutate(color = tol21rainbow[c(1:8)]) %>% 
    column_to_rownames(var = "class")


png("results/figures/main_classes_totalEVs_blood_atleast2.5perc.png", res = 300, units = "in", height = 8, width = 6)
main_class_all %>% 
    group_by(matrix) %>% 
    mutate(total = sum(hits)) %>% 
    ungroup() %>% 
    mutate(prop = hits/total*100) %>% 
    mutate(class = case_when(prop > 2.5 ~ Name,
                             prop <= 2.5 ~ "z_Other")) %>% 
    group_by(matrix, class) %>% 
    mutate(hits.e = sum(hits)) %>% 
    ungroup() %>% 
    group_by(matrix) %>% 
    distinct(class, .keep_all = TRUE) %>% 
    ungroup() %>% 
    ggplot() +
    geom_col(aes(x = matrix, y = hits.e, fill = class), color = "grey20", width = 0.6, position = "fill") +
    scale_fill_manual(values = names_main$color) +
    theme_bw() +
    theme(legend.text = element_text(size=7),
          legend.key.size = unit(0.5, 'cm')) +
    guides(fill = guide_legend(ncol = 1)) +
    labs(x = "", y = "Proportion of hits", fill = "Class")
dev.off()



# Sub class ----
ndev_sub <- read_csv("results/NDEV/enrichment_gt_median/ndev_sub_class/msea_ora_result.csv") %>% 
    mutate(matrix = "NDEV") %>% 
    dplyr::rename(Name = 1)
TotalEV_sub <- read_csv("results/NDEV/enrichment_gt_median/Totalev_sub_class/msea_ora_result.csv") %>% 
    mutate(matrix = "TotalEV") %>% 
    dplyr::rename(Name = 1)
serum_sub <- read_csv("results/NDEV/enrichment_gt_median/serum_sub_class/msea_ora_result.csv") %>% 
    mutate(matrix = "Serum") %>% 
    dplyr::rename(Name = 1)

sub_class_all <- rbind(ndev_sub, TotalEV_sub, serum_sub) 

sub_class_all %>% 
    group_by(matrix) %>% 
    mutate(total = sum(hits)) %>% 
    ungroup() %>% 
    mutate(prop = hits/total*100) %>% 
    ggplot(aes(x = prop)) +
    geom_histogram(binwidth = 0.5)

tol21rainbow= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")
pal(tol21rainbow)

names_sub <- sub_class_all %>% 
    group_by(matrix) %>% 
    mutate(total = sum(hits)) %>% 
    ungroup() %>% 
    mutate(prop = hits/total*100) %>% 
    mutate(class = case_when(prop > 2.5 ~ Name,
                             prop <= 2.5 ~ "z_Other")) %>% 
    distinct(class) %>% 
    mutate(color = tol21rainbow[c(1:11)]) %>% 
    column_to_rownames(var = "class")


png("results/figures/sub_classes_totalEVs_blood_atleast2.5perc.png", res = 300, units = "in", height = 8, width = 6)
sub_class_all %>% 
    group_by(matrix) %>% 
    mutate(total = sum(hits)) %>% 
    ungroup() %>% 
    mutate(prop = hits/total*100) %>% 
    mutate(class = case_when(prop > 2.5 ~ Name,
                             prop <= 2.5 ~ "z_Other")) %>% 
    group_by(matrix, class) %>% 
    mutate(hits.e = sum(hits)) %>% 
    ungroup() %>% 
    group_by(matrix) %>% 
    distinct(class, .keep_all = TRUE) %>% 
    ungroup() %>% 
    ggplot() +
    geom_col(aes(x = matrix, y = hits.e, fill = class), color = "grey20", width = 0.6, position = "fill") +
    scale_fill_manual(values = names_sub$color) +
    theme_bw() +
    theme(legend.text = element_text(size=7),
          legend.key.size = unit(0.5, 'cm')) +
    guides(fill = guide_legend(ncol = 1)) +
    labs(x = "", y = "Proportion of hits", fill = "Class")
dev.off()
