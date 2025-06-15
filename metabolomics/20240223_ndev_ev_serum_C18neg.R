#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# NDEV, EV (exoQ), and serum 
# Compare levels across the three matrices
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(tidyverse)
library(gplots)
library(RColorBrewer)
replacezero <- function(x) "[<-"(x, x == 0, min(x[x>0])/2)


# Call in filtration information ----
filter.val <- read_tsv("results/blank_filtration/C18_neg_exo_sec.txt")

# Call in C18 + data ----
C18_neg <- read_tsv("C18_column/RAW_mzcalibrated_untargeted_featuretable_C18neg.txt")
new_names <- gsub(".mzXML", "", names(C18_neg))
names(C18_neg) <- new_names
map <- read_tsv("C18_column/map_edited_C18_neg.txt") %>% 
    mutate(type_source = case_when(Sample_type == "Water" ~ "Water",
                                   Sample_type == "NIST" ~ "NIST",
                                   Sample_type == "HRE_Pool" ~ "HRE_Pool",
                                   type_2 == "SEC blank" ~ "Blank",
                                   type_2 == "Exo blank" ~ "Blank",
                                   grepl("P", Sample.ID) ~ "Plasma",
                                   grepl("S", Sample.ID) ~ "Serum",
                                   grepl("Plasma", Sample.ID) ~ "Plasma (total)",
                                   grepl("Serum", Sample.ID) ~ "Serum (total)",
                                   grepl("NDX 7", Sample.ID) ~ "Serum",
                                   grepl("NDX 8", Sample.ID) ~ "Serum",
                                   grepl("NDX 9", Sample.ID) ~ "Serum",
                                   grepl("NDX 10", Sample.ID) ~ "Serum",
                                   grepl("NDX 11", Sample.ID) ~ "Serum",
                                   grepl("NDX 12", Sample.ID) ~ "Plasma",
                                   grepl("NDX 13", Sample.ID) ~ "Plasma",
                                   grepl("NDX14", Sample.ID) ~ "Plasma",
                                   grepl("NDX15", Sample.ID) ~ "Plasma",
                                   grepl("NDX16", Sample.ID) ~ "Plasma")) %>% 
    mutate(type_method = case_when(Sample_type == "Water" ~ "Water",
                                   Sample_type == "NIST" ~ "NIST",
                                   Sample_type == "HRE_Pool" ~ "HRE_Pool",
                                   type_2 == "SEC blank" ~ "Blank",
                                   type_2 == "Exo blank" ~ "Blank",
                                   grepl("ExoQ", Sample.ID) ~ "ExoQ",
                                   grepl("SEC", Sample.ID) ~ "SEC",
                                   grepl("NDX 7|NDX 8|NDX 9", Sample.ID) ~ "NDEV dry beads",
                                   grepl("NDX 10", Sample.ID) ~ "NDEV glycine",
                                   grepl("NDX 11", Sample.ID) ~ "NDEV elution",
                                   grepl("NDX 12|NDX 13|NDX14", Sample.ID) ~ "NDEV dry beads",
                                   grepl("NDX15", Sample.ID) ~ "NDEV glycine",
                                   grepl("NDX16", Sample.ID) ~ "NDEV elution",
                                   grepl("Plasma", Sample.ID) ~ "Blood",
                                   grepl("Serum", Sample.ID) ~ "Blood"))

# Filter metabolomic data ----
filter_exo <- filter.val %>% 
    filter(gt_2_sd_exo == 1)

C18_neg_exo <- C18_neg %>% 
    unite("mz_time", c(mz, time), sep = "_", remove = FALSE) %>% 
    filter(mz_time %in% filter_exo$mz_time)

# File names for comparing blood, total EVs, and NDEVs
C18_neg_map_compare <- map %>% 
    filter(type_source == "Serum") %>% 
    filter(!grepl("100", type_2)) %>% 
    filter(!grepl("200", type_2)) %>% 
    filter(!grepl("SEC", type_2)) %>% 
    mutate(type_compare = case_when(grepl("NDEV", type_broad) ~ "nEV",
                                    grepl("ExoQ", type_broad) ~ "Total EV",
                                    grepl("Matrix", type_broad) ~ "Serum")) %>% 
    filter(type_method != "NDEV glycine") %>% 
    filter(type_method != "NDEV elution") %>% 
    dplyr::select(File.Name, type_compare)

C18_neg_compare <- C18_neg_exo %>% 
    dplyr::select(mz_time, C18_neg_map_compare$File.Name)

# Find features with differential abundance in the comparison groups
C18_neg_compare_labels <- C18_neg_compare %>% 
    t() %>% 
    as.data.frame() %>% 
    janitor::row_to_names(1) %>% 
    rownames_to_column(var = "File.Name") %>% 
    right_join(C18_neg_map_compare, ., by = "File.Name")

# Features most abundant in NDEVs
C18_neg_ndev <- C18_neg_compare_labels %>% 
    filter(type_compare == "nEV") %>% 
    dplyr::select(contains("_"), -type_compare) %>% 
    map_df(~as.character(.x)) %>% 
    map_df(~as.numeric(.x)) %>% 
    map_df(~mean(.x)) %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "feature") %>% 
    dplyr::rename(mean.int = V1) %>% 
    arrange(-mean.int) %>% 
    filter(mean.int != 0)

C18_neg_ndev %>% write_tsv("results/NDEV/C18_neg_serum_features_NDEV.txt")

# Features most abundant in total ev
C18_neg_tev <- C18_neg_compare_labels %>% 
    filter(type_compare == "Total EV") %>% 
    dplyr::select(contains("_"), -type_compare) %>% 
    map_df(~as.character(.x)) %>% 
    map_df(~as.numeric(.x)) %>% 
    map_df(~mean(.x)) %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "feature") %>% 
    dplyr::rename(mean.int = V1) %>% 
    arrange(-mean.int) %>% 
    filter(mean.int != 0)

C18_neg_tev %>% write_tsv("results/NDEV/C18_neg_serum_features_TotalEV.txt")

# Features most abundant in serum
C18_neg_serum <- C18_neg_compare_labels %>% 
    filter(type_compare == "Serum") %>% 
    dplyr::select(contains("_"), -type_compare) %>% 
    map_df(~as.character(.x)) %>% 
    map_df(~as.numeric(.x)) %>% 
    map_df(~mean(.x)) %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "feature") %>% 
    dplyr::rename(mean.int = V1) %>% 
    arrange(-mean.int) %>% 
    filter(mean.int != 0)

C18_neg_serum %>% write_tsv("results/NDEV/C18_neg_serum_features_serum.txt")

C18_neg_mean_all <- C18_neg_compare_labels %>% 
    gather(key = "feature", value = "intensity", -type_compare, -File.Name) %>% 
    mutate_at("intensity", as.character) %>% 
    mutate_at("intensity", as.numeric) %>% 
    group_by(type_compare, feature) %>% 
    summarise(mean.int = mean(intensity)) %>% 
    ungroup() %>% 
    arrange(-mean.int)

C18_neg_mean <- C18_neg_compare_labels %>% 
    gather(key = "feature", value = "intensity", -type_compare, -File.Name) %>% 
    mutate_at("intensity", as.character) %>% 
    mutate_at("intensity", as.numeric) %>% 
    group_by(feature) %>%
    summarise(mean.int = mean(intensity))  %>% 
    ungroup()

zero_feature <- C18_neg_mean %>% 
    filter(mean.int == 0)

# Select features most abundant in nEVs ----
all_top_100 <- C18_neg_mean_all %>% 
    arrange(-mean.int) %>% 
    slice(1:100)

# Plot mean intensity 
C18_neg_mean_all %>% 
    filter(feature %in% all_top_100$feature) %>% 
    separate(feature, into = c("mz", "time"), sep = "_", remove = FALSE) %>% 
    ggplot(aes(x = type_compare, y = time, fill = log2(mean.int))) +
    geom_tile(color = "white") +
    theme_bw()

# Plot heatmap
rows <- C18_neg_compare_labels %>% 
    group_by(type_compare) %>% 
    mutate(row_name = paste0(type_compare, "_", seq_along(type_compare))) %>% 
    ungroup() %>% 
    .$row_name

ndev_heatmap <- C18_neg_compare_labels %>% 
    group_by(type_compare) %>% 
    mutate(row_name = paste0(type_compare, "_", seq_along(type_compare))) %>% 
    ungroup() %>% 
    dplyr::select(row_name, all_top_100$feature) %>% 
    column_to_rownames(var = "row_name") %>% 
    map_df(~as.character(.x)) %>% 
    map_df(~as.numeric(.x)) %>% 
    map_df(~scale(.x))

rownames(ndev_heatmap) <- rows

ndev_scaled_heat <- ndev_heatmap %>% 
    t() %>% 
    as.data.frame()

colMain <- colorRampPalette(brewer.pal(8, "RdYlBu"))(7)

heatmap.2(as.matrix(ndev_scaled_heat),
          Rowv=TRUE,
          Colv=TRUE,
          trace='none',
          density.info="none",
          dendrogram = "col",
          #ColSideColors = colSide,
          #labRow = "",
          #labCol = "",
          #xlab = "", ylab =  "",
          col = colMain,
          keysize = 1,
          key.title = "",
          #margins=c(8,12),
          cexRow=1.3,
          cexCol=1.3) 

ndev_heatmap$group <- rows

ndev_scaled_tile <- ndev_heatmap %>% 
    gather(key = "feature", value = "intensity", -group)

ndev_scaled_tile %>% 
    ggplot(aes(x = group, y = reorder(feature, intensity), fill = intensity)) +
    geom_tile(color = "white") +
    theme_bw()

# PCA
zero_feature <- C18_neg_mean %>% 
    filter(mean.int == 0)

met_pca <- C18_neg_compare_labels %>% 
    dplyr::select(-File.Name, -type_compare, -all_of(zero_feature$feature)) %>% 
    map_df(~as.character(.x)) %>% 
    map_df(~as.numeric(.x)) %>% 
    map_df(~replacezero(.x)) %>% 
    prcomp(.)

png("results/NDEV/C18neg_PCA.png", res = 300, units = "in", h = 3, w = 4)
factoextra::fviz_pca_ind(met_pca, axes = c(1,2), habillage = hilic_pos_compare_labels$type_compare, 
                         geom = c("point"), title = "C18 -") 
dev.off()


C18_neg_compare_labels %>% 
    dplyr::select(-File.Name, -type_compare, -all_of(zero_feature$feature)) %>% 
    map_df(~as.character(.x)) %>% 
    map_df(~as.numeric(.x)) %>% 
    map_df(~replacezero(.x)) %>% 
    mutate(File.Name = C18_neg_compare_labels$File.Name) %>% 
    mutate(type_compare = C18_neg_compare_labels$type_compare) %>% 
    select(File.Name, type_compare, everything()) %>% 
    write_tsv("nev_metabolomics/C18neg_nev_serum_tev_replaced.txt")


