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

# Main class ----
ndev_main <- read_csv("results/NDEV/enrichment_gt_90perc/ndev_main_no_level5/msea_ora_result.csv") %>% 
    mutate(matrix = "NDEV") %>% 
    dplyr::rename(Name = 1)
Totalev_main <- read_csv("results/NDEV/enrichment_gt_90perc/Totalev_main_no_level5/msea_ora_result.csv") %>% 
    mutate(matrix = "TotalEV") %>% 
    dplyr::rename(Name = 1)
serum_main <- read_csv("results/NDEV/enrichment_gt_90perc/serum_main_no_level5/msea_ora_result.csv") %>% 
    mutate(matrix = "Serum") %>% 
    dplyr::rename(Name = 1)

main_class_all <- rbind(ndev_main, Totalev_main, serum_main)

hist_classes <- main_class_all %>% 
    group_by(matrix) %>% 
    mutate(total = sum(hits)) %>% 
    ungroup() %>% 
    mutate(prop = hits/total*100) %>% 
    ggplot(aes(x = prop)) +
    geom_histogram(binwidth = 0.5)

tol21rainbow= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")

names_main <- main_class_all %>% 
    filter(FDR<0.2) %>% 
    distinct(Name) %>% 
    mutate(color = c(tol21rainbow[c(1:8, 16:20)], "grey86")) %>% 
    column_to_rownames(var = "Name")

main_class_all$matrix <- factor(main_class_all$matrix, levels = c("Serum", "TotalEV", "NDEV"),
                                labels = c("Serum", "TotalEV", "nEV"))

main_class_plot <- main_class_all %>% 
    group_by(matrix) %>% 
    mutate(total = sum(hits)) %>% 
    ungroup() %>% 
    mutate(prop = hits/total*100) %>% 
    mutate(class = case_when(prop > 2.5 ~ Name,
                             prop <= 2.5 ~ "Other")) %>% 
    group_by(matrix, class) %>% 
    mutate(hits.e = sum(hits)) %>% 
    ungroup() %>% 
    group_by(matrix) %>% 
    distinct(class, .keep_all = TRUE) %>% 
    ungroup()

main_class_plot$class <- factor(main_class_plot$class, levels = c("Benzene and substituted derivatives",
                                                                  "Carboxylic acids and derivatives",
                                                                  "Fatty Acyls",
                                                                  "Hydroxy acids and derivatives",
                                                                  "Imidazopyrimidines",
                                                                  "Indoles and derivatives",
                                                                  "Organonitrogen compounds",
                                                                  "Organooxygen compounds",
                                                                  "Prenol lipids",
                                                                  "Pyrrolidines",
                                                                  "Steroids and steroid derivatives",
                                                                  "Unsaturated hydrocarbons",
                                                                  "Other"),
                                labels = c("Benzene and substituted \nderivatives",
                                           "Carboxylic acids and derivatives",
                                           "Fatty Acyls",
                                           "Hydroxy acids and derivatives",
                                           "Imidazopyrimidines",
                                           "Indoles and derivatives",
                                           "Organonitrogen compounds",
                                           "Organooxygen compounds",
                                           "Prenol lipids",
                                           "Pyrrolidines",
                                           "Steroids and steroid derivatives",
                                           "Unsaturated hydrocarbons",
                                           "Other"))

#png("results/NDEV/to10perc_main_classes_totalEVs_blood_atleast2.5perc.png", 
    #res = 300, units = "in", height = 7, width = 2.5)
main_class_plot %>% 
    ggplot() +
    geom_col(aes(x = matrix, y = hits.e, fill = Name), color = "grey20", width = 0.6, position = "fill") +
    #scale_fill_manual(values = names_main$color) +
    theme_bw() +
    theme(legend.text = element_text(size=7),
          legend.key.size = unit(0.5, 'cm'),
          legend.position = "bottom") +
    guides(fill = guide_legend(ncol = 2)) +
    labs(x = "", y = "Proportion of hits", fill = "")
#dev.off()

classes_prop <- main_class_plot %>% 
    ggplot() +
    geom_col(aes(x = matrix, y = hits.e, fill = class), color = "grey20", width = 0.6, position = "fill") +
    scale_fill_manual(values = names_main$color) +
    theme_bw() +
    theme(legend.text = element_text(size=7),
          legend.key.size = unit(0.5, 'cm'),
          legend.position = "bottom") +
    labs(x = "", y = "Proportion of hits", fill = "")

other_classes <- main_class_all %>% 
    mutate(prop = hits/total*100) %>%
    filter(prop <= 2.5)

main_classes_all_fdr <- main_class_all %>% 
    filter(FDR < 0.05)

main_classes_all_holm <- main_class_all %>% 
    filter(`Holm p` < 0.05)

main_classes_all_rawp <- main_class_all %>% 
    filter(`Raw p` < 0.05)

main_classes_all_rawp %>% 
    distinct(Name)

# Add indication for FDR < 0.05?
plot_data <- main_class_all %>% 
    mutate(label.sig = case_when(FDR < 0.05 ~ "*",
                                 FDR >= 0.05 ~ "")) %>% 
    filter(FDR < 0.2) %>%
    group_by(matrix) %>% 
    mutate(sum.hits = sum(hits)) %>% 
    ungroup() %>% 
    mutate(prop.hits = hits/sum.hits)

png("results/NDEV/main_classes_totalEVs_blood_FDR020.png", 
    res = 300, units = "in", height = 4, width = 4.5)
plot_data %>% 
    ggplot(aes(x = matrix, y = prop.hits, fill = Name, label = label.sig)) +
    geom_bar(stat = "identity", color = "grey20", width = 0.6, position = "fill") +
    geom_text(size = 3, position = position_stack(vjust = 0.5), color = "grey86") +
    scale_fill_manual(values = names_main$color) +
    theme_bw() +
    theme(legend.text = element_text(size=7),
          legend.key.size = unit(0.5, 'cm'),
          legend.position = "right") +
    guides(fill = guide_legend(ncol = 1)) +
    labs(x = "", y = "Proportion of hits", fill = "")

dev.off()

classes_prop_2 <- plot_data %>% 
    ggplot(aes(x = matrix, y = prop.hits, fill = Name, label = label.sig)) +
    geom_bar(stat = "identity", color = "grey20", width = 0.6, position = "fill") +
    geom_text(size = 5, position = position_stack(vjust = 0.45), color = "grey90") +
    scale_fill_manual(values = names_main$color) +
    theme_bw() +
    theme(legend.text = element_text(size=7),
          legend.key.size = unit(0.5, 'cm'),
          legend.position = "bottom") +
    guides(fill = guide_legend(ncol = 1)) +
    labs(x = "", y = "Proportion of hits", fill = "")

