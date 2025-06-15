#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Venn overlap 
# 2024 February 21
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(tidyverse)
library(VennDetail)
library(VennDiagram)
library(ggvenn)
library(ggVennDiagram)
library(venn)
library(openxlsx)
library(RColorBrewer)

# Call in the summary data
all_chem <- read.xlsx("exposomics/results/chemicals_detection.xlsx", sheet = 1)
none_chem <- read.xlsx("exposomics/results/chemicals_detection.xlsx", sheet = 3)
brain <- read.xlsx("exposomics/results/chemicals_detection.xlsx", sheet = 6)
ndev <- read.xlsx("exposomics/results/chemicals_detection.xlsx", sheet = 7)
ev <- read.xlsx("exposomics/results/chemicals_detection.xlsx", sheet = 9)
serum <- read.xlsx("exposomics/results/chemicals_detection.xlsx", sheet = 10)
dep <- read.xlsx("exposomics/results/chemicals_detection.xlsx", sheet = 8)

groups<- read_tsv("exposomics/summary_chemicals_group_prop_updated.txt") %>% 
    rename(Detail = chemical)

brain_chem <- brain %>% 
    filter(Number.of.individuals > 2) %>% 
    .$Chemical
ndev_chem <- ndev %>% 
    filter(Number.of.individuals > 2) %>% 
    .$Chemical
ev_chem <- ev %>% 
    filter(Number.of.individuals > 2) %>% 
    .$Chemical
serum_chem <- serum %>% 
    filter(Number.of.individuals > 2) %>% 
    .$Chemical
dep_chem <- dep %>% 
    filter(Number.of.individuals > 2) %>% 
    .$Chemical

chem <- list(Brain = brain_chem,
             NDEV = ndev_chem,
             TotalEV = ev_chem,
             Serum = serum_chem,
             Depl = dep_chem)

ven <- venndetail(chem)

ggvenn(chem, 
       fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF", "orange"),
       stroke_size = 0.5, set_name_size = 4)

ggVennDiagram(chem, edge_size = 2, label = "count", label_alpha = 0) + scale_fill_distiller(palette = "RdBu")

venn(chem, ilab=TRUE, zcolor = "style", ilcs = 1)

plot(ven, type = "upset")

ven_result <- result(ven, wide = TRUE)

dplot(ven, order = TRUE, textsize = 4)

# Merge the venn result file with the detection rates

plot_venn_1 <- brain %>%     
    filter(Number.of.individuals > 2) %>%
    dplyr::rename(Detail = Chemical, brain_num = Number.of.individuals) %>% 
    left_join(ven_result, ., by = "Detail")

plot_venn_2 <- ndev %>%     
    filter(Number.of.individuals > 2) %>%
    dplyr::rename(Detail = Chemical, ndev_num = Number.of.individuals) %>% 
    left_join(plot_venn_1, ., by = "Detail")

plot_venn_3 <- ev %>%     
    filter(Number.of.individuals > 2) %>%
    dplyr::rename(Detail = Chemical, ev_num = Number.of.individuals) %>% 
    left_join(plot_venn_2, ., by = "Detail")

plot_venn_4 <- serum %>%     
    filter(Number.of.individuals > 2) %>%
    dplyr::rename(Detail = Chemical, serum_num = Number.of.individuals) %>% 
    left_join(plot_venn_3, ., by = "Detail")

plot_venn_5 <- dep %>%     
    filter(Number.of.individuals > 2) %>%
    dplyr::rename(Detail = Chemical, dep_num = Number.of.individuals) %>% 
    left_join(plot_venn_4, ., by = "Detail")

plot_venn_long <- plot_venn_5 %>% 
    gather(key = "matrix", value = "num", brain_num, ndev_num, ev_num, serum_num) %>% 
    filter(!is.na(num))

plot_venn_long$matrix <- factor(plot_venn_long$matrix, levels = c("brain_num", "serum_num", "ev_num", "ndev_num"),
                                labels = c("Brain \nTissue", "Serum", "TotalEV", "nEV"))

plot_groups <- left_join(plot_venn_long, groups, by = "Detail")

plot_groups$chem.group.2 <- factor(plot_groups$chem.group.2, levels = c("dioxin", "PCB", "PAH", 
                                                                        "organochlorine", "organophosphate",
                                                                        "pyrethroid", "triaza", "other"),
                                   labels = c("Dioxins", "PCBs", "PAHs", "Organochlorines", "Organophosphates",
                                              "Pyrethroids", "Triaza- compounds", "Others"))

png("figures/20240221_chemical_detection_51.png", res = 300, units = "in", h = 8.5, w = 5.5)
plot_groups %>% 
    ggplot(aes(x = matrix, y = reorder(Detail, SharedSets), fill = factor(num))) +
    geom_tile(color = "white") +
    theme_bw() +
    scale_fill_manual(values = c("#FFD580", "coral", "indianred4")) +
    labs(x = "", y = "", fill = "Detected in \nsamples") +
    ggforce::facet_col(facets = vars(`chem.group.2`),
                       scales = "free_y", space = "free") +
    theme(strip.text.y = element_text(angle=0))
dev.off()

detection_plot <- plot_groups %>% 
    ggplot(aes(x = matrix, y = reorder(Detail, SharedSets), fill = factor(num))) +
    geom_tile(color = "white") +
    theme_bw() +
    scale_fill_manual(values = c("#FFD580", "coral", "indianred4")) +
    labs(x = "", y = "", fill = "Detected in \nsamples") +
    ggforce::facet_col(facets = vars(`chem.group.2`),
                       scales = "free_y", space = "free") +
    theme(strip.text.y = element_text(angle=0))

png("figures/figure3_GC_detection_correlation.png", res = 300, units = "in", h = 10, w = 12)
cowplot::plot_grid(detection_plot, all_chem_corr, labels = c("A", "B"),
                   nrow = 1)
dev.off()

