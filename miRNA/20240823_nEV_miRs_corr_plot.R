#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Visualize pathways for nEVs miR correlated with brain
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(tidyverse)

# Call in the data
high <- read_csv("HTG/results/mirpath/nEVs_high_corr_miRPath-v4-miRNA-centric_analysis_2024-02-20_union.csv")
low <- read_csv("HTG/results/mirpath/nEVs_low_corr_miRPath-v4-miRNA-centric_analysis_2024-02-20_union.csv")

high_arrange_30 <- high %>% 
    arrange(-`miRNAs (n)`) %>% 
    slice(1:30) %>% 
    mutate(corr = "Positive") %>% 
    mutate(tot = 88) %>% 
    mutate(prop = `miRNAs (n)`/tot)

low_arrange_30 <- low %>% 
    arrange(-`miRNAs (n)`) %>% 
    slice(1:30) %>% 
    mutate(corr = "Negative") %>% 
    mutate(tot = 86) %>% 
    mutate(prop = `miRNAs (n)`/tot)

both <- rbind(high_arrange_30, low_arrange_30)
both$corr <- factor(both$corr, levels = c("Positive", "Negative"))
png("figures/figure2_nEVs_pathways.png", res = 300, units = "in",
    h = 11, w = 6)
both %>% 
    ggplot(aes(x = corr, y = reorder(`Term Name`, `miRNAs (n)`))) +
    geom_point(aes(color = prop, size = -log10(FDR))) +
    theme_bw() +
    scale_size_binned() +
    scale_color_gradient(high = "maroon", low = "orange") +
    labs(x = "Correlation direction", y = "", size = "-Log10(q)", color = "Proportion \n of miRNA")
dev.off()

nEV_corrs <- both %>% 
    ggplot(aes(x = corr, y = reorder(`Term Name`, `miRNAs (n)`))) +
    geom_point(aes(color = prop, size = -log10(FDR))) +
    theme_bw() +
    scale_size_binned() +
    scale_color_gradient(high = "maroon", low = "orange") +
    labs(x = "Correlation direction", y = "", size = "-Log10(q)", color = "Proportion \n of miRNA")

# For supplemental table
high_all <- high %>% 
    mutate(correlation = "positive") %>% 
    mutate(total = 88) %>% 
    mutate(prop = `miRNAs (n)`/total) %>% 
    arrange(-`miRNAs (n)`)

low_all <- low %>% 
    mutate(correlation = "negative") %>% 
    mutate(total = 86) %>% 
    mutate(prop = `miRNAs (n)`/total)%>% 
    arrange(-`miRNAs (n)`)

both_all <- rbind(high_all, low_all)

both_all %>% write_tsv("HTG/results/miR_pathways_nEV_corrs.txt")
