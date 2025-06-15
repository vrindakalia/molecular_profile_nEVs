#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# correlations for chemicals in all
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(tidyverse)
library(corrplot)
library(RColorBrewer)

# data
dat <- read_tsv("exposomics/clean_data/in_four_replaced_atl3.txt")
chem_key <- read_tsv("exposomics/clean_data/chemical_key.txt")

c25 <- c(
    "dodgerblue2", "#E31A1C", # red
    "green4",
    "#6A3D9A", # purple
    "#FF7F00", # orange
    "black", 
    "gold1",
    "skyblue2", "#FB9A99", # lt pink
    "palegreen2",
    "#CAB2D6", # lt purple
    "#FDBF6F", # lt orange
    "gray70", "khaki2",
    "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
    "darkturquoise", "green1", "yellow4", "yellow3",
    "darkorange4", "brown", "forestgreen", "royalblue3", "indianred4"
)


dat %>% 
    separate(ID, into = c("source", "sample"), sep = "_", remove = FALSE) %>% 
    gather(key = "index", value = "level",  -source, -sample, -ID) %>% 
    mutate_at("level", as.character) %>% 
    mutate_at("level", as.numeric) %>% 
    ggplot(aes(x = index, y = log2(level), color = source)) +
    geom_boxplot() +
    geom_line(aes(group = source))

# Consider correlation with all chemicals, grouped by chemicals and sample
# Paired across chemicals and people
dat.all <- dat %>% 
    separate(ID, into = c("source", "sample"), sep = "_", remove = FALSE) %>% 
    select(-ID) %>% 
    gather(key = "index", value = "level",  -source, -sample) %>% 
    mutate_at("level", as.character) %>% 
    mutate_at("level", as.numeric)

# NDEV

dat.all.ndev <- dat.all %>% 
    filter(source %in% c("Tissue", "NDEV")) %>% 
    spread(key = "source", value = "level") #%>% 
    #filter(log2(NDEV) > -12)
cor.test( ~ NDEV + Tissue, data = dat.all.ndev,  method = "spearm", exact = FALSE)
    
dat.all.ndev %>% 
    ggplot(aes(x = log2(NDEV), y = log2(Tissue), color = index)) +
    geom_point() +
    geom_smooth(method = "lm") + 
    theme_bw()

# Total EV
dat.all.ev <- dat.all %>% 
    filter(source %in% c("Tissue", "TotalEV")) %>% 
    spread(key = "source", value = "level") #%>% 
    #filter(log2(TotalEV) > -12)
cor.test( ~ TotalEV + Tissue, data = dat.all.ev,  method = "spearm", exact = FALSE)

dat.all.ev %>% 
    ggplot(aes(x = log2(TotalEV), y = log2(Tissue), color = index)) +
    geom_point() +
    geom_smooth(method = "lm") +
    theme_bw()

# Serum
dat.all.serum <- dat.all %>% 
    filter(source %in% c("Tissue", "Serum")) %>% 
    spread(key = "source", value = "level") %>% 
    filter(log2(Serum) > -12)
cor.test( ~ Serum + Tissue, data = dat.all.serum,  method = "spearm", exact = FALSE)

dat.all.serum %>% 
    ggplot(aes(x = log2(Serum), y = log2(Tissue), color = index)) +
    geom_point() +
    geom_smooth(method = "lm") +
    theme_bw()

# Depleted
dat.all.depl <- dat.all %>% 
    filter(source %in% c("Tissue", "Depleted")) %>% 
    spread(key = "source", value = "level")
cor.test( ~ Depleted + Tissue, data = dat.all.depl,  method = "spearm", exact = FALSE)

dat.all.depl %>% 
    ggplot(aes(y = log2(Depleted), x = log2(Tissue), color = index)) +
    geom_point() +
    geom_smooth(method = "lm") +
    theme_bw()

tissue_all <- left_join(dat.all.serum, dat.all.ev, by = c("index", "sample", "Tissue")) %>% 
    left_join(., dat.all.ndev, by = c("index", "sample", "Tissue")) %>% 
    gather(key = "matrix", value = "level", NDEV, Serum, TotalEV) %>% 
    merge(., chem_key, by = "index")

tissue_all$matrix <- factor(tissue_all$matrix, levels = c("Serum", "TotalEV", "NDEV"), 
                            labels = c("Serum", "TotalEV", "nEV"))

ann_text <- data.frame(
    label = c("corr = 0.64", "corr = 0.68", "corr = 0.64"),
    matrix =  c("Serum", "nEV", "TotalEV"),
    x     = c(rep(-20, 3)),
    y     = c(rep(5, 3)))

png("figures/gc_total_corr.png", res = 300, units = "in", h = 7, w = 5)
tissue_all %>% 
    filter(log2(level) > -12) %>% 
    ggplot(aes(x = log2(Tissue), y = log2(level))) +
    geom_smooth(method = "lm", color = "grey30", linetype = "dotted", fill = "grey86", size = 0.7) +
    geom_point(aes(color = reorder(Chemical, -level))) +
    geom_smooth(aes(color = reorder(Chemical, -level)), method = "lm", fill = "grey76") +
    theme_bw() +
    facet_wrap(~matrix, nrow = 3, scales = "free_y") +
    theme(legend.position = "right") +
    scale_color_manual(values = c25) +
    #geom_text(data = ann_text,  aes(x = -Inf, y = -Inf, label = label),
    #          hjust   = -0.1,
    #          vjust   = -1) +
    labs(x = "Log2(Tissue concentration)", y = "Log2(Matrix concentration)", color = "")
dev.off()

all_chem_corr <- tissue_all %>% 
    filter(log2(level) > -12) %>% 
    ggplot(aes(x = log2(Tissue), y = log2(level))) +
    geom_smooth(method = "lm", color = "grey30", linetype = "dotted", fill = "grey86", size = 0.7) +
    geom_point(aes(color = reorder(Chemical, -level))) +
    geom_smooth(aes(color = reorder(Chemical, -level)), method = "lm", fill = "grey76") +
    theme_bw() +
    facet_wrap(~matrix, nrow = 3, scales = "free_y") +
    theme(legend.position = "right") +
    scale_color_manual(values = c25) +
    #geom_text(data = ann_text,  aes(x = -Inf, y = -Inf, label = label),
    #          hjust   = -0.1,
    #          vjust   = -1) +
    labs(x = "Log2(Brain Tissue Concentration)", y = "Log2(Compartment Concentration)", color = "Chemical")

tissue_all_depl <-  left_join(dat.all.serum, dat.all.ev, by = c("index", "sample", "Tissue")) %>% 
    left_join(., dat.all.ndev, by = c("index", "sample", "Tissue")) %>% 
    left_join(., dat.all.depl, by = c("index", "sample", "Tissue")) %>% 
    gather(key = "matrix", value = "level", NDEV, Serum, TotalEV, Depleted) %>% 
    merge(., chem_key, by = "index")

tissue_all_depl$matrix <- factor(tissue_all_depl$matrix, levels = c("Serum", "TotalEV", "NDEV", "Depleted"))

ann_text_depl <- data.frame(
    label = c("corr = 0.70", "corr = 0.72", "corr = 0.57", "corr = 0.33"),
    matrix =  c("Serum", "NDEV", "TotalEV", "Depleted"),
    x     = c(rep(-20, 4)),
    y     = c(rep(5, 4)))

png("figures/gc_total_corr_depleted.png", res = 300, units = "in", h = 4, w = 8)

tissue_all_depl %>% 
    #filter(log2(level) > -36) %>% 
    ggplot(aes(x = log2(Tissue), y = log2(level))) +
    geom_point(aes(color = reorder(Chemical, -level))) +
    geom_smooth(aes(color = reorder(Chemical, -level)), method = "lm") +
    theme_bw() +
    facet_wrap(~matrix, nrow = 1, scales = "free_y") +
    theme(legend.position = "bottom") +
    scale_color_manual(values = c25) +
    #geom_text(data = ann_text_depl,  aes(x = -Inf, y = -Inf, label = label),
    #          hjust   = -0.1,
    #          vjust   = -1) +
    labs(x = "Log2(Tissue concentration)", y = "Log2(Matrix concentration)", color = "")
dev.off()

# Brain and NDEV
dat_brain_ndev <- dat %>% 
    separate(ID, into = c("source", "sample"), sep = "_", remove = FALSE) %>% 
    filter(source %in% c("Tissue", "NDEV")) %>% 
    select(-ID) %>% 
    gather(key = "index", value = "level",  -source, -sample) %>% 
    mutate_at("level", as.character) %>% 
    mutate_at("level", as.numeric) 
    #%>% 
    #spread(key = "source", value = "level")

chem_indices = unique(dat_brain_ndev$index)

brain_ndev_result_list <- NULL

for(i in 1:length(chem_indices)){
    
    brain_ndev_result_list[[i]] <- wilcox.test(level ~ source, data = filter(dat_brain_ndev, index == chem_indices[i]), paired = TRUE, alternative = "two.sided")
    
}

brain_ndev_result_list_tidy <- brain_ndev_result_list %>% 
    map_df(broom::tidy) %>% 
    mutate(index = chem_indices)

dat_brain_ndev %>% 
    ggplot(aes(x = index, y = log2(level), color = source)) +
    geom_boxplot() +
    geom_line(aes(group = source)) #+
    #facet_wrap(~sample, scales = "free_y")

# Brain and serum
dat_brain_serum <- dat %>% 
    separate(ID, into = c("source", "sample"), sep = "_", remove = FALSE) %>% 
    filter(source %in% c("Tissue", "Serum")) %>% 
    select(-ID) %>% 
    gather(key = "index", value = "level",  -source, -sample) %>% 
    mutate_at("level", as.character) %>% 
    mutate_at("level", as.numeric) #%>% 
    #spread(key = "source", value = "level")


brain_serum_result_list <- NULL

for(i in 1:length(chem_indices)){
    
    brain_serum_result_list[[i]] <- wilcox.test(level ~ source, data = filter(dat_brain_serum, index == chem_indices[i]), paired = TRUE, alternative = "two.sided")
    
}

brain_serum_result_list_tidy <- brain_serum_result_list %>% 
    map_df(broom::tidy) %>% 
    mutate(index = chem_indices)

dat_brain_serum %>% 
    ggplot(aes(x = source, y = level, group = sample, color = sample)) +
    geom_point() +
    geom_line() +
    facet_wrap(~index, scales = "free_y")

# Brain and total EVs
dat_brain_EV <- dat %>% 
    separate(ID, into = c("source", "sample"), sep = "_", remove = FALSE) %>% 
    filter(source %in% c("Tissue", "TotalEV")) %>% 
    select(-ID) %>% 
    gather(key = "index", value = "level",  -source, -sample) %>% 
    mutate_at("level", as.character) %>% 
    mutate_at("level", as.numeric) #%>% 
#spread(key = "source", value = "level")

brain_EV_result_list <- NULL

for(i in 1:length(chem_indices)){
    
    brain_EV_result_list[[i]] <- wilcox.test(level ~ source, data = filter(dat_brain_EV, index == chem_indices[i]), paired = TRUE, alternative = "two.sided")
    
}

brain_EV_result_list_tidy <- brain_EV_result_list %>% 
    map_df(broom::tidy) %>% 
    mutate(index = chem_indices)

dat_brain_EV %>% 
    ggplot(aes(x = source, y = level, group = sample, color = sample)) +
    geom_point() +
    geom_line() +
    facet_wrap(~index, scales = "free_y")

# Brain and depleted
dat_brain_deplete <- dat %>% 
    separate(ID, into = c("source", "sample"), sep = "_", remove = FALSE) %>% 
    filter(source %in% c("Tissue", "Depleted")) %>% 
    select(-ID) %>% 
    gather(key = "index", value = "level",  -source, -sample)  %>% 
    mutate_at("level", as.character) %>% 
    mutate_at("level", as.numeric) #%>% 
#spread(key = "source", value = "level")


brain_deplete_result_list <- NULL

for(i in 1:length(chem_indices)){
    
    brain_deplete_result_list[[i]] <- wilcox.test(level ~ source, data = filter(dat_brain_deplete, index == chem_indices[i]), paired = TRUE, alternative = "two.sided")
    
}

brain_deplete_result_list_tidy <- brain_deplete_result_list %>% 
    map_df(broom::tidy) %>% 
    mutate(index = chem_indices)

all_brain_cor <- merge(ndev_brain_correlations, deplete_brain_correlations, by = "index") %>% 
    merge(., EV_brain_correlations, by = "index") %>% 
    merge(., serum_brain_correlations, by = "index")

all_brain_cor_plot <- all_brain_cor %>% 
    select(-contains("p.value")) %>% 
    gather(key = "type", value = "estimate", -index) %>% 
    separate(type, into = c("matrix", "parameter")) %>% 
    merge(., chem_key, by = "index")

all_brain_pval_plot <- all_brain_cor %>% 
    select(-contains("cor")) %>% 
    gather(key = "type", value = "pvalue", -index) %>% 
    separate(type, into = c("matrix", "parameter")) %>% 
    merge(., chem_key, by = "index")
