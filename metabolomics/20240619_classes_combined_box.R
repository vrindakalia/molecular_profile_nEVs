#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Combine metabolite class information across the columns
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(tidyverse)
source("code/20240619_visualize_gt_90perc.R")
# Call in the files
hilpos_classes_box <- read_tsv("results/NDEV/20240619_hilpos_nev_tev_serum_classes_for_boxes") %>% 
    mutate(column = "HILIC pos")
hilneg_classes_box <- read_tsv("results/NDEV/20240619_hilneg_nev_tev_serum_classes_for_boxes") %>% 
    mutate(column = "HILIC neg")
c18pos_classes_box <- read_tsv("results/NDEV/20240619_c18pos_nev_tev_serum_classes_for_boxes") %>% 
    mutate(column = "C18 pos")
c18neg_classes_box <- read_tsv("results/NDEV/20240619_c18neg_nev_tev_serum_classes_for_boxes") %>% 
    mutate(column = "C18 neg")

all_classes_box <- rbind(hilpos_classes_box, hilneg_classes_box,
                         c18pos_classes_box, c18neg_classes_box)

all_classes_box$type_compare <- factor(all_classes_box$type_compare,
                                       levels = c("Serum", "Total EV", "nEV"),
                                       labels = c("Serum", "TotalEV", "nEV"))

png("results/NDEV/2024_all_classes_matrices_boxes.png",
    res = 300, units = "in", h = 7, w = 11)
all_classes_box %>% 
    filter(`Main class` != "Alkyl fluorides") %>% 
    ggplot(aes(x = type_compare, y = log2(intensity), fill = type_compare)) +
    geom_boxplot(width = 0.5, alpha = 0.8) +
    facet_wrap(~`Main class`, scales = "free_y", ncol = 6) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(x = "") +
    scale_fill_manual(values = c("orange", "#B66D99", "#2E71AF"))
dev.off()

boxes <- all_classes_box %>% 
    filter(`Main class` != "Alkyl fluorides") %>% 
    ggplot(aes(x = type_compare, y = log2(intensity), fill = type_compare)) +
    geom_boxplot(width = 0.5, alpha = 0.8) +
    facet_wrap(~`Main class`, scales = "free_y", ncol = 6) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 8)) +
    labs(x = "", y = "Log2(Intensity)") +
    scale_fill_manual(values = c("orange", "#B66D99", "#2E71AF"))

png("results/figures/figure4_classes_boxes.png", res = 300, units = "in",
    h = 7, w = 13.5)
cowplot::plot_grid(classes_prop_2, boxes, nrow = 1, labels = c("A", "B"), rel_widths = c(0.2, 1))
dev.off()

png("results/NDEV/2024_all_classes_matrices_bars.png",
    res = 300, units = "in", h = 8, w = 12)
all_classes_box %>% 
    ggplot(aes(x = type_compare, fill = type_compare)) +
    geom_bar(stat = "count", width = 0.6) +
    facet_wrap(~`Main class`, scales = "free_y", ncol = 6) +
    theme_bw()+
    theme(legend.position = "none") +
    labs(x = "") +
    scale_fill_manual(values = c("#87321F", "#B66D99", "#2E71AF"))
dev.off()

png("results/NDEV/2024_bile_acids_matrices_boxes.png",
    res = 300, units = "in", h = 5, w = 6)
all_classes_box %>% 
    filter(`Main class` == "Bile acids") %>% 
    ggplot(aes(x = type_compare, y = log2(intensity), fill = type_compare)) +
    geom_boxplot(width = 0.6) +
    facet_wrap(~`Name`, scales = "free_y", ncol = 3) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(x = "", title = "Bile acids") +
    scale_fill_manual(values = c("#87321F", "#B66D99", "#2E71AF"))
dev.off()

png("results/NDEV/2024_fatty_amides_matrices_boxes.png",
    res = 300, units = "in", h = 3, w = 4)
all_classes_box %>% 
    filter(`Main class` == "Fatty amides") %>% 
    ggplot(aes(x = type_compare, y = log2(intensity), fill = type_compare)) +
    geom_boxplot(width = 0.6) +
    facet_wrap(~`Name`, scales = "free_y", ncol = 3) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(x = "", title = "Fatty amides") +
    scale_fill_manual(values = c("#87321F", "#B66D99", "#2E71AF"))
dev.off()

png("results/NDEV/2024_porphyrins_matrices_boxes.png",
    res = 300, units = "in", h = 3, w = 4)
all_classes_box %>% 
    filter(`Main class` == "Porphyrins") %>% 
    ggplot(aes(x = type_compare, y = log2(intensity), fill = type_compare)) +
    geom_boxplot(width = 0.6) +
    facet_wrap(~`Name`, scales = "free_y", ncol = 3) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(x = "", title = "Porphyrins") +
    scale_fill_manual(values = c("#87321F", "#B66D99", "#2E71AF"))
dev.off()

png("results/NDEV/2024_sphingoid_bases_matrices_boxes.png",
    res = 300, units = "in", h = 3, w = 4)
all_classes_box %>% 
    filter(`Main class` == "Sphingoid bases") %>% 
    ggplot(aes(x = type_compare, y = log2(intensity), fill = type_compare)) +
    geom_boxplot(width = 0.6) +
    facet_wrap(~`Name`, scales = "free_y", ncol = 3) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(x = "", title = "Sphingoid bases") +
    scale_fill_manual(values = c("#87321F", "#B66D99", "#2E71AF"))
dev.off()

png("results/NDEV/2024_aa_peptides_matrices_boxes.png",
    res = 300, units = "in", h = 8, w = 12)
all_classes_box %>% 
    filter(`Main class` == "Amino acids and peptides") %>% 
    ggplot(aes(x = type_compare, y = log2(intensity), fill = type_compare)) +
    geom_boxplot(width = 0.6) +
    facet_wrap(~`Name`, scales = "free_y", ncol = 5) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(x = "", title = "Amino acids and peptides") +
    scale_fill_manual(values = c("#87321F", "#B66D99", "#2E71AF"))
dev.off()

png("results/NDEV/2024_fatty_acids_matrices_boxes.png",
    res = 300, units = "in", h = 6, w = 8)
all_classes_box %>% 
    filter(`Main class` == "Fatty acids") %>% 
    ggplot(aes(x = type_compare, y = log2(intensity), fill = type_compare)) +
    geom_boxplot(width = 0.6) +
    facet_wrap(~`Name`, scales = "free_y", ncol = 3) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(x = "", title = "Fatty acids") +
    scale_fill_manual(values = c("#87321F", "#B66D99", "#2E71AF"))
dev.off()

png("results/NDEV/2024_PFOS_matrices_boxes.png",
    res = 300, units = "in", h = 3, w = 4)
all_classes_box %>% 
    filter(`Main class` == "Alkyl fluorides") %>% 
    ggplot(aes(x = type_compare, y = log2(intensity), fill = type_compare)) +
    geom_boxplot(width = 0.6) +
    facet_wrap(~`Name`, scales = "free_y", ncol = 3) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(x = "", title = "Alkyl fluorides") +
    scale_fill_manual(values = c("#87321F", "#B66D99", "#2E71AF"))
dev.off()

pfos <- all_classes_box %>% 
    filter(`Main class` == "Alkyl fluorides") %>% 
    select(type_compare, intensity, File.Name) %>%
    mutate(sample = stringr::str_sub(File.Name, start = -3)) %>% 
    mutate(sample_num = case_when(sample == "012" ~ "1",
                                  sample == "014" ~ "1",
                                  sample == "026" ~ "2",
                                  sample == "036" ~ "3",
                                  sample == "038" ~ "2",
                                  sample == "057" ~ "3",
                                  sample == "064" ~ "1",
                                  sample == "071" ~ "2",
                                  sample == "076" ~ "4",
                                  sample == "078" ~ "3",
                                  sample == "083" ~ "4")) %>% 
    select(-File.Name, -sample) %>% 
    spread(key = "type_compare", value = "intensity")

pfos %>% 
    ggplot(aes(x = Serum, y = TotalEV)) +
    geom_point() +
    geom_smooth(method = "lm") +
    theme_bw()

all_classes_box %>% write_tsv("results/NDEV/20240619_classes_all_columns.txt")

table(all_classes_box$type_compare, all_classes_box$`Main class`)

png("results/NDEV/20240619_main_class_member_numbers.png", 
    res = 300, units = "in", h = 5.5, w = 3)
all_classes_box %>% 
    group_by(`Main class`) %>% 
    summarise(num.mem = n()/11) %>% 
    ggplot(aes(x = reorder(`Main class`, num.mem), y = num.mem)) +
    geom_bar(stat = "identity") +
    theme_bw() +
    coord_flip() +
    labs(y = "Number of members", x = "")
dev.off()
