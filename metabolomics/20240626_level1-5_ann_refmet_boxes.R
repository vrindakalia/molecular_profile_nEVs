#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Boxplots all main classes with level 1 - 5 annotations
# 2024.06.26
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(tidyverse)

ann <- read_tsv("results/NDEV/all_matrix_hmdbids.txt")
ann_refmet <- read_tsv("results/NDEV/all_matrix_refmet_results.txt")

ann_classes <- merge(ann, ann_refmet, by.x = "use_name", by.y = "Input name")

# Call in the data
hilpos <- read_tsv("nev_metabolomics/hilicpos_nev_serum_tev_replaced.txt")
labels <- select(hilpos, File.Name, type_compare)
hilpos_t <- hilpos %>%
    select(-type_compare) %>% 
    t() %>% as.data.frame() %>% 
    janitor::row_to_names(1) %>% 
    rownames_to_column(var = "feature")

# Fix variable class type
hilpos_cmp <- hilpos_t %>% 
    select(-feature) %>% 
    map_df(~as.character(.x)) %>% 
    map_df(~as.numeric(.x)) %>% 
    mutate(feature = hilpos_t$feature) %>% 
    select(feature, everything())

# Turn into long format
hilpos_long <- hilpos_cmp %>% 
    gather(key = "File.Name", value = "intensity", -feature) %>% 
    merge(., labels, by = "File.Name")

# Call in the data
c18pos <- read_tsv("nev_metabolomics/c18pos_nev_serum_tev_replaced.txt")
labels <- select(c18pos, File.Name, type_compare)
c18pos_t <- c18pos %>%
    select(-type_compare) %>% 
    t() %>% as.data.frame() %>% 
    janitor::row_to_names(1) %>% 
    rownames_to_column(var = "feature")

# Fix variable class type
c18pos_cmp <- c18pos_t %>% 
    select(-feature) %>% 
    map_df(~as.character(.x)) %>% 
    map_df(~as.numeric(.x)) %>% 
    mutate(feature = c18pos_t$feature) %>% 
    select(feature, everything())

# Turn into long format
c18pos_long <- c18pos_cmp %>% 
    gather(key = "File.Name", value = "intensity", -feature) %>% 
    merge(., labels, by = "File.Name")

# Call in the data
c18neg <- read_tsv("nev_metabolomics/c18neg_nev_serum_tev_replaced.txt")
labels <- select(c18neg, File.Name, type_compare)
c18neg_t <- c18neg %>%
    select(-type_compare) %>% 
    t() %>% as.data.frame() %>% 
    janitor::row_to_names(1) %>% 
    rownames_to_column(var = "feature")

# Fix variable class type
c18neg_cmp <- c18neg_t %>% 
    select(-feature) %>% 
    map_df(~as.character(.x)) %>% 
    map_df(~as.numeric(.x)) %>% 
    mutate(feature = c18neg_t$feature) %>% 
    select(feature, everything())

# Turn into long format
c18neg_long <- c18neg_cmp %>% 
    gather(key = "File.Name", value = "intensity", -feature) %>% 
    merge(., labels, by = "File.Name") 

# Call in the data
hilneg <- read_tsv("nev_metabolomics/hilicneg_nev_serum_tev_replaced.txt")
labels <- select(hilneg, File.Name, type_compare)
hilneg_t <- hilneg %>%
    select(-type_compare) %>% 
    t() %>% as.data.frame() %>% 
    janitor::row_to_names(1) %>% 
    rownames_to_column(var = "feature")

# Fix variable class type
hilneg_cmp <- hilneg_t %>% 
    select(-feature) %>% 
    map_df(~as.character(.x)) %>% 
    map_df(~as.numeric(.x)) %>% 
    mutate(feature = hilneg_t$feature) %>% 
    select(feature, everything())

# Turn into long format
hilneg_long <- hilneg_cmp %>% 
    gather(key = "File.Name", value = "intensity", -feature) %>% 
    merge(., labels, by = "File.Name") 

# Combine data from four columns
all_long <- rbind(hilpos_long, hilneg_long, c18pos_long, c18neg_long) %>% 
    rename(mz_time = feature)

# Merge with annotation and class information
all_int_classes <- left_join(all_long, ann_classes, by = "mz_time")

all_classes_box <- all_int_classes %>% 
    filter(!is.na(`Main class`))

#png("results/NDEV/20240618_c18neg_MainClass_boxplot.png", res = 300, units = "in",
#    h = 8, w = 8)
all_classes_box %>% 
    ggplot(aes(x = type_compare, y = log2(intensity))) +
    geom_boxplot() +
    facet_wrap(~`Main class`, scales = "free_y") +
    theme_bw()
#dev.off()

all_classes_box %>% 
    filter(confidence != "Level 5") %>% 
    ggplot(aes(x = type_compare, y = log2(intensity))) +
    geom_boxplot() +
    facet_wrap(~`Main class`, scales = "free_y") +
    theme_bw()

boxes <- all_classes_box %>% 
    filter(confidence != "Level 5") 

boxes$type_compare <- factor(boxes$type_compare,
                                       levels = c("Serum", "Total EV", "nEV"),
                                       labels = c("Serum", "TotalEV", "nEV"))

png("results/NDEV/20240626_main_classes_level1-3_boxes_all_columns.png", res = 300,
    units = "in", h = 16, w = 16)
boxes %>% 
    ggplot(aes(x = type_compare, y = log2(intensity), fill = type_compare)) +
    geom_boxplot(width = 0.5, alpha = 0.8) +
    facet_wrap(~`Main class`, scales = "free_y", ncol = 9) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 8)) +
    labs(x = "", y = "Log2(Intensity)") +
    scale_fill_manual(values = c("orange", "#B66D99", "#2E71AF"))
dev.off()

png("results/NDEV/20240626_super_classes_level1-3_boxes_all_columns.png", res = 300,
    units = "in", h = 8, w = 6.25)
boxes_fig <- boxes %>% 
    ggplot(aes(x = type_compare, y = log2(intensity), fill = type_compare)) +
    geom_boxplot(width = 0.5, alpha = 0.8) +
    facet_wrap(~`Super class`, scales = "free_y", ncol = 3) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 8)) +
    labs(x = "", y = "Log2(Intensity)") +
    scale_fill_manual(values = c("orange", "#B66D99", "#2E71AF"))
dev.off()

boxes %>% 
    filter(`Main class` == "Bile acids") %>% 
    ggplot(aes(x = type_compare, y = log2(intensity), fill = type_compare)) +
    geom_boxplot(width = 0.6) +
    facet_wrap(~`Name`, scales = "free_y", ncol = 3) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(x = "", title = "Bile acids") +
    scale_fill_manual(values = c("orange", "#B66D99", "#2E71AF"))

# Total number of metabolites in each sub group
subgroup_nums <- boxes %>% 
    distinct(use_name, .keep_all = TRUE) %>%
    mutate(sub_class_clean=case_when(`Sub class` == "NAE" ~ "N-acyl ethanolamines",
                                     `Sub class` == "PS" ~" Phosphoserine",
                                     `Sub class` == "Unsaturated FA" ~ "Unsaturated fatty acid",
                                     `Sub class` == "Saturated FA" ~ "Saturated fatty acid",
                                     `Sub class` != "NAE" & `Sub class` != "PS" & `Sub class` != "Unsaturated FA" & `Sub class` != "Saturated FA" ~ `Sub class`)) %>% 
    group_by(sub_class_clean) %>% 
    summarise(total_num = n()) %>% 
    ungroup()

# ANOVAs and posthoc Tukey test for super class membership    
boxes_test <- boxes %>% 
    mutate_at("intensity", log2)

unique_super_class <- boxes_test %>% 
    distinct(`Super class`)

results_sup <- vector(mode = "list", length = dim(unique_super_class)[1])
names(results_sup) <- unique_super_class$`Super class`

for(i in unique_super_class$`Super class`){
    sub <- boxes_test %>% 
        filter(`Super class` == i)
    
    results_sup[[i]] <- aov(intensity ~ type_compare, data = sub)
    
}

summary(results_sup[[2]])

results_sup_tidy <- results_sup %>% 
    map_df(broom::tidy) %>% 
    filter(term == "type_compare") %>% 
    mutate(name = unique_super_class$`Super class`) %>% 
    mutate(fdr = p.adjust(p.value, method = "BH"))

posthoc_results_sup <- results_sup %>% 
    map(~TukeyHSD(.x))

posthoc_results_sup_tidy <- posthoc_results_sup %>% 
    map_df(broom::tidy)

posthoc_results_sup_nEV_serum <- posthoc_results_sup_tidy %>% 
    filter(contrast == "nEV-Serum") %>% 
    mutate(name = unique_super_class$`Super class`) %>% 
    mutate(qvalue = p.adjust(adj.p.value, method = "BH"))

posthoc_results_sup_nEV_serum_sig <- posthoc_results_sup_nEV_serum %>% 
    filter(qvalue < 0.05)

ndev_higher_sup <- posthoc_results_sup_nEV_serum_sig %>% 
    filter(estimate >0)

ndev_lower_sup <- posthoc_results_sup_nEV_serum_sig %>% 
    filter(estimate < 0)

posthoc_results_sup_nEV_totalEV <- posthoc_results_sup_tidy %>% 
    filter(contrast == "nEV-TotalEV") %>% 
    mutate(name = unique_super_class$`Super class`) %>% 
    mutate(qvalue = p.adjust(adj.p.value, method = "BH"))

posthoc_results_sup_nEV_totalEV_sig <- posthoc_results_sup_nEV_totalEV %>% 
    filter(qvalue < 0.05)

posthoc_results_sup_nEV_totalEV_notsig <- posthoc_results_sup_nEV_totalEV %>% 
    filter(qvalue > 0.05)

ndev_higher_tev_sup <- posthoc_results_sup_nEV_totalEV_sig %>% 
    filter(estimate >0)

ndev_lower_tev_sup <- posthoc_results_sup_nEV_totalEV_sig %>% 
    filter(estimate < 0)

posthoc_results_sup_totalEV_serum <- posthoc_results_sup_tidy %>% 
    filter(contrast == "TotalEV-Serum") %>% 
    mutate(name = unique_super_class$`Super class`) %>% 
    mutate(qvalue = p.adjust(adj.p.value, method = "BH"))

posthoc_results_sup_totalEV_serum_sig <- posthoc_results_sup_totalEV_serum %>% 
    filter(qvalue < 0.05)

posthoc_results_sup_totalEV_serum_notsig <- posthoc_results_sup_totalEV_serum %>% 
    filter(qvalue > 0.05)

overlap_sup <- VennDetail::venndetail(list("nEV-serum" = posthoc_results_sup_nEV_serum_sig$name, 
                                           "nEV-tEV" = posthoc_results_sup_nEV_totalEV_sig$name,
                                           "tEV-serum" = posthoc_results_sup_totalEV_serum_sig$name))

plot(overlap_sup)

posthoc_results_sup_all <- rbind(posthoc_results_sup_nEV_serum,posthoc_results_sup_nEV_totalEV,posthoc_results_sup_totalEV_serum)
posthoc_results_sup_all %>% write_tsv("results/NDEV/20240726_superclass_ANOVA_posthoc.txt")

# ANOVAs and posthoc Tukey test
boxes_test <- boxes %>% 
    mutate_at("intensity", log2)

# Number of metabolites
unique_annotations <- boxes_test %>% 
    distinct(use_name)

results <- vector(mode = "list", length = dim(unique_annotations)[1])
names(results) <- unique_annotations$use_name

for(i in unique_annotations$use_name){
    sub <- boxes_test %>% 
        filter(use_name == i)
    
    results[[i]] <- aov(intensity ~ type_compare, data = sub)
    
}

summary(results[[2]])

results_tidy <- results %>% 
    map_df(broom::tidy) %>% 
    filter(term == "type_compare") %>% 
    mutate(name = unique_annotations$use_name)

posthoc_results <- results %>% 
    map(~TukeyHSD(.x))

posthoc_results_tidy <- posthoc_results %>% 
    map_df(broom::tidy)

posthoc_results_nEV_serum <- posthoc_results_tidy %>% 
    filter(contrast == "nEV-Serum") %>% 
    mutate(name = unique_annotations$use_name) %>% 
    mutate(qvalue = p.adjust(adj.p.value, method = "BH"))

posthoc_results_nEV_serum_sig <- posthoc_results_nEV_serum %>% 
    filter(qvalue < 0.05)

ndev_higher <- posthoc_results_nEV_serum_sig %>% 
    filter(estimate >0)

ndev_lower <- posthoc_results_nEV_serum_sig %>% 
    filter(estimate < 0)


ndev_higher %>% write_tsv("results/NDEV/comparisons/20240628_nEV-serum_higher.txt")
ndev_lower %>% write_tsv("results/NDEV/comparisons/20240628_nEV-serum_lower.txt")

posthoc_results_nEV_totalEV <- posthoc_results_tidy %>% 
    filter(contrast == "nEV-TotalEV") %>% 
    mutate(name = unique_annotations$use_name) %>% 
    mutate(qvalue = p.adjust(adj.p.value, method = "BH"))

posthoc_results_nEV_totalEV_sig <- posthoc_results_nEV_totalEV %>% 
    filter(qvalue < 0.05)

ndev_higher_tev <- posthoc_results_nEV_totalEV_sig %>% 
    filter(estimate >0)

ndev_lower_tev <- posthoc_results_nEV_totalEV_sig %>% 
    filter(estimate < 0)

posthoc_results_totalEV_serum <- posthoc_results_tidy %>% 
    filter(contrast == "TotalEV-Serum") %>% 
    mutate(name = unique_annotations$use_name) %>% 
    mutate(qvalue = p.adjust(adj.p.value, method = "BH"))

posthoc_results_totalEV_serum_sig <- posthoc_results_totalEV_serum %>% 
    filter(qvalue < 0.05)

plot(overlap_sup)

posthoc_results_metab_all <- rbind(posthoc_results_nEV_serum, posthoc_results_nEV_totalEV, posthoc_results_totalEV_serum)
posthoc_results_metab_all %>% write_tsv("results/NDEV/20240726_individual_metab_ANOVA_posthoc.txt")

# donut plot with higher lower information
# Create test data.
data_1 <- data.frame(category=c("nEV", "Serum"),
                     count=c(31, 414))

# Compute percentages
data_1$fraction = data_1$count / sum(data_1$count)
# Compute the cumulative percentages (top of each rectangle)
data_1$ymax = cumsum(data_1$fraction)
# Compute the bottom of each rectangle
data_1$ymin = c(0, head(data_1$ymax, n=-1))
# Compute label position
data_1$labelPosition <- (data_1$ymax + data_1$ymin) / 2
# Compute a good label
data_1$label <- paste0(data_1$category, ": ", data_1$count)
# Make the plot
png("results/NDEV/20240627_donut_nEV_Serum.png", res = 300, units = "in",
    h = 4, w= 4)
do_1 <- ggplot(data_1, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
    geom_rect() +
    geom_label( x=4.25, aes(y=labelPosition, label=label), size=3, color = "white") +
    scale_fill_manual(values = c("#2E71AF", "orange")) +
    coord_polar(theta="y") +
    xlim(c(2, 4)) +
    theme_void() +
    theme(legend.position = "none") +
    annotate(geom = 'text', x = 2, y = 0, label = "Total: \n445", size = 3)
dev.off()

# nEV - TotalEV

data_2 <- data.frame(category=c("nEV", "TotalEV"),
                     count=c(39, 237))

# Compute percentages
data_2$fraction = data_2$count / sum(data_2$count)
# Compute the cumulative percentages (top of each rectangle)
data_2$ymax = cumsum(data_2$fraction)
# Compute the bottom of each rectangle
data_2$ymin = c(0, head(data_2$ymax, n=-1))
# Compute label position
data_2$labelPosition <- (data_2$ymax + data_2$ymin) / 2
# Compute a good label
data_2$label <- paste0(data_2$category, ": ", data_2$count)
# Make the plot
png("results/NDEV/20240627_donut_nEV_TotalEV.png", res = 300, units = "in",
    h = 4, w= 4)
do_2 <- ggplot(data_2, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
    geom_rect() +
    geom_label(x = 4.25, aes(y=labelPosition, label=label), size=3, color = "white") +
    scale_fill_manual(values = c("#2E71AF", "#B66D99")) +
    coord_polar(theta="y") +
    xlim(c(2, 4)) +
    theme_void() +
    theme(legend.position = "none") +
    annotate(geom = 'text', x = 2, y = 0, label = "Total: \n276", size = 3)
dev.off()

# Compare nEV-Serum higher to nEV-tEV higher
ndev_higher$name
ndev_higher_tev$name

overlap <- VennDetail::venndetail(list("Higher than \nSerum" = ndev_higher$name, 
                                       "Higher than \nTotalEV" = ndev_higher_tev$name))

plot(overlap)

overlap_res <- overlap@result

common <- filter(overlap_res, Subset == "Shared")

common %>% write_tsv("results/NDEV/20240628_higher_in_ndev_shared.txt")

boxes_common <- boxes %>% 
    filter(use_name %in% overlap_res$Detail) %>% 
    mutate(sub_class_clean=case_when(`Sub class` == "NAE" ~ "N-acyl ethanolamines",
                                     `Sub class` == "PS" ~" Phosphoserine",
                                     `Sub class` == "Unsaturated FA" ~ "Unsaturated fatty acid",
                                     `Sub class` == "Saturated FA" ~ "Saturated fatty acid",
                                     `Sub class` != "NAE" & `Sub class` != "PS" & `Sub class` != "Unsaturated FA" & `Sub class` != "Saturated FA" ~ `Sub class`)) 

png("results/NDEV/20240628_boxes_all_nEV_higher_subclass.png", res = 300, units = "in",
    h = 10, w = 14)
boxes_common %>% 
    ggplot(aes(x = type_compare, y = log2(intensity), fill = type_compare)) +
    geom_boxplot(width = 0.5, alpha = 0.8) +
    facet_wrap(~`sub_class_clean` + use_name, nrow = 8, scales = "free_y") +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 8)) +
    labs(x = "", y = "Log2(Intensity)") +
    scale_fill_manual(values = c("orange", "#B66D99", "#2E71AF"))
dev.off()

# HMDBIDs of metabolites higher in nEV compared to serum and tEV
boxes %>% 
    distinct(use_name, .keep_all = TRUE) %>% 
    filter(use_name %in% ndev_higher$name) %>% 
    write_tsv("results/NDEV/comparisons/20240628_nEV_higher_than_serum_IDs.txt")

boxes %>% 
    distinct(use_name, .keep_all = TRUE) %>% 
    filter(use_name %in% ndev_higher_tev$name) %>% 
    write_tsv("results/NDEV/comparisons/20240628_nEV_higher_than_totalEVs_IDs.txt")

boxes %>% 
    distinct(use_name, .keep_all = TRUE) %>% 
    filter(use_name %in% ndev_lower$name) %>% 
    write_tsv("results/NDEV/comparisons/20240628_nEV_lower_than_serum_IDs.txt")

boxes %>% 
    distinct(use_name, .keep_all = TRUE) %>% 
    filter(use_name %in% ndev_lower_tev$name) %>% 
    write_tsv("results/NDEV/comparisons/20240628_nEV_lower_than_totalEVs_IDs.txt")

# Create bar plot with data from enrichment analysis of metabolites higher in nEVs
nev_higher_serum <- read_csv("results/NDEV/comparisons/nEV_higher_than_serum_main_class/msea_pie_data.csv") %>% 
    mutate(compare = "nEV-Serum")
nev_higher_totalev <- read_csv("results/NDEV/comparisons/nEV_higher_than_TotalEV_main_class/msea_pie_data.csv") %>% 
    mutate(compare = "nEV-TotalEV")

both_comb <- rbind(nev_higher_serum, nev_higher_totalev)

both_comb %>% 
    ggplot(aes(x = Hits, y = reorder(X1, Hits))) +
    geom_bar(stat = "identity", fill = "#2E71AF", width = 0.6) +
    theme_bw() +
    facet_wrap(~compare) +
    labs(x = "Number of metabolites", y = "")

both_comb %>% 
    ggplot(aes(x = Hits, y = reorder(X1, Hits))) +
    geom_segment(aes(xend = 0, yend = X1)) +
    geom_point(size = 2, color = "#2E71AF") +
    theme_bw() +
    facet_wrap(~compare) +
    labs(x = "Number of metabolites", y = "")

# Check against refmet classification
nev_higher_serum_ref <- read_tsv("results/NDEV/comparisons/20240628_nEV_higher_than_serum_IDs.txt") %>% 
    mutate(compare = "nEV-Serum")
nev_higher_totalev_ref <- read_tsv("results/NDEV/comparisons/20240628_nEV_higher_than_totalEVs_IDs.txt") %>% 
    mutate(compare = "nEV-TotalEV")

both_comb_ref <- rbind(nev_higher_serum_ref, nev_higher_totalev_ref) %>% 
    mutate(sub_class_clean=case_when(`Sub class` == "NAE" ~ "N-acyl ethanolamines",
                                     `Sub class` == "PS" ~" Phosphoserine",
                                     `Sub class` == "Unsaturated FA" ~ "Unsaturated fatty acid",
                                     `Sub class` == "Saturated FA" ~ "Saturated fatty acid",
                                     `Sub class` != "NAE" & `Sub class` != "PS" & `Sub class` != "Unsaturated FA" & `Sub class` != "Saturated FA" ~ `Sub class`)) %>% 
    group_by(compare, sub_class_clean) %>% 
    summarise(num_metab = n()) %>% 
    ungroup() %>% 
    merge(., subgroup_nums, by = "sub_class_clean") %>% 
    mutate(ratio = paste0(num_metab, "/", total_num)) %>% 
    mutate(ratio_num = num_metab/total_num)

both_comb_ref %>% 
    ggplot(aes(x = num_metab, y = reorder(sub_class_clean, num_metab))) +
    geom_bar(stat = "identity", fill = "#2E71AF", width = 0.6) +
    geom_text(aes(label = ratio), hjust = -0.1, size = 3) +
    theme_bw() +
    facet_wrap(~compare) +
    labs(x = "Number of metabolites", y = "") +
    xlim(0,4)

png("results/NDEV/20240628_nev_higher_sub_classes.png", res = 300,
    units = "in", h = 5, w = 5)
lo <- both_comb_ref %>% 
    ggplot(aes(x = num_metab, y = reorder(sub_class_clean, num_metab))) +
    geom_segment(aes(xend = 0, yend = sub_class_clean)) +
    geom_point() + #aes(color = ratio_num)) +
    #geom_text(aes(label = ratio), hjust = -0.5, size = 3) +
    theme_bw() + 
    #scale_color_gradient(high = "#2E71AF", low = "grey86")+
    facet_wrap(~compare) +
    labs(x = "Number of metabolites", y = "", color = "Proportion of total detected") +
    theme(legend.position = "bottom") +
    guides(color = guide_legend(title.position="top", override.aes = list(size=4))) 
dev.off()

left_top <- cowplot::plot_grid(do_1, do_2, labels = c("C", "D"), nrow = 1)
left <- cowplot::plot_grid(left_top, lo, labels = c("","E"), rel_heights = c(0.3, 1), nrow = 2)

source("code/20240619_visualize_gt_90perc.R")
classes_prop_leg <- classes_prop_2 +
    labs(fill = "Metabolite main class") +
    guides(fill = guide_legend(title.position="top", title.hjust = 0.5, ncol = 1))

png("results/figures/figure4_classes_boxes_compare.png", res = 300, units = "in",
    h = 8, w = 13.5)
cowplot::plot_grid(boxes_fig, classes_prop_leg, left, nrow = 1, labels = c("A", "B", ""), 
                   rel_widths = c(0.475, 0.175, 0.35))
dev.off()


# Annotations summary
all_classes_box_uni <- all_classes_box %>% 
    filter(confidence == "Level 3") %>% 
    distinct(use_name, .keep_all = TRUE)

tot_feat <- 6956+2280+6257+2366
1673/tot_feat

ann_level3 <- ann %>% 
    filter(confidence == "Level 3")
