#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal library annotations
# 2024, Februrary 20
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(tidyverse)
library(xMSanalyzer)

# C18 NEG
# Call in the annotation files
#c18_neg_lib <- read_csv("xMSannotator/C18 NEG/Stage5.csv")
c18_neg_lib_2 <- read_tsv("results/C18neg_annotations_level1-5.txt")

# Call in the feature list of interest
ev_C18neg <- read_tsv("results/NDEV/C18_neg_serum_features_TotalEV.txt")

dataA <- ev_C18neg %>% 
    separate(feature, into = c("mz", "time"), sep = "_", remove = FALSE) %>% 
    dplyr::select(mz, time, everything()) %>% 
    mutate_at("mz", as.numeric) %>% 
    mutate_at("time", as.numeric) %>% 
    mutate(index.A = seq(mz)) 

dataB <- c18_neg_lib_2 %>% 
    dplyr::select(mz, time, everything()) %>% 
    mutate_at("mz", as.numeric) %>% 
    mutate_at("time", as.numeric) %>% 
    mutate(index.B = seq(mz)) 

overlap <- find.Overlapping.mzs(dataA, dataB, mz.thresh = 10, time.thresh = 10,
                                alignment.tool=NA)

lib_overlap <- dataB[overlap$index.B,]
data_overlap <- dataA[overlap$index.A,]

data_lib <- cbind(data_overlap, lib_overlap, overlap)

data_lib %>% write_tsv("results/NDEV/level1-5_c18_neg_features_totalEV.txt")


# C18 POS
# Call in the annotation files
#c18_neg_lib <- read_csv("xMSannotator/C18 NEG/Stage5.csv")
c18_pos_lib_2 <- read_tsv("results/C18pos_annotations_level1-5.txt")

# Call in the feature list of interest
ev_C18pos <- read_tsv("results/NDEV/C18_pos_serum_features_TotalEV.txt")

dataA <- ev_C18pos %>% 
    separate(feature, into = c("mz", "time"), sep = "_", remove = FALSE) %>% 
    dplyr::select(mz, time, everything()) %>% 
    mutate_at("mz", as.numeric) %>% 
    mutate_at("time", as.numeric) %>% 
    mutate(index.A = seq(mz)) 

dataB <- c18_pos_lib_2 %>% 
    dplyr::select(mz, time, everything()) %>% 
    mutate_at("mz", as.numeric) %>% 
    mutate_at("time", as.numeric) %>% 
    mutate(index.B = seq(mz)) 

overlap <- find.Overlapping.mzs(dataA, dataB, mz.thresh = 10, time.thresh = 10,
                                alignment.tool=NA)

lib_overlap <- dataB[overlap$index.B,]
data_overlap <- dataA[overlap$index.A,]

data_lib <- cbind(data_overlap, lib_overlap, overlap)

data_lib %>% write_tsv("results/NDEV/level1-5_c18_pos_features_TotalEV.txt")

# HILIC POS
# Call in the annotation files
#c18_neg_lib <- read_csv("xMSannotator/C18 NEG/Stage5.csv")
hilic_pos_lib_2 <- read_tsv("results/hilpos_annotations_level1-5.txt")

# Call in the feature list of interest
ev_hilicpos <- read_tsv("results/NDEV/hilic_pos_serum_features_TotalEV.txt")

dataA <- ev_hilicpos %>% 
    separate(feature, into = c("mz", "time"), sep = "_", remove = FALSE) %>% 
    dplyr::select(mz, time, everything()) %>% 
    mutate_at("mz", as.numeric) %>% 
    mutate_at("time", as.numeric) %>% 
    mutate(index.A = seq(mz)) 

dataB <- hilic_pos_lib_2 %>% 
    dplyr::select(mz, time, everything()) %>% 
    mutate_at("mz", as.numeric) %>% 
    mutate_at("time", as.numeric) %>% 
    mutate(index.B = seq(mz)) 

overlap <- find.Overlapping.mzs(dataA, dataB, mz.thresh = 10, time.thresh = 10,
                                alignment.tool=NA)

lib_overlap <- dataB[overlap$index.B,]
data_overlap <- dataA[overlap$index.A,]

data_lib <- cbind(data_overlap, lib_overlap, overlap)

data_lib %>% write_tsv("results/NDEV/level1-5_hilic_pos_features_TotalEV.txt")

# HILIC NEG
# Call in the annotation files
#c18_neg_lib <- read_csv("xMSannotator/C18 NEG/Stage5.csv")
hilic_neg_lib_2 <- read_tsv("results/hilneg_annotations_level1-5.txt")

# Call in the feature list of interest
ev_hilicneg <- read_tsv("results/NDEV/hilic_neg_serum_features_TotalEV.txt")

dataA <- ev_hilicneg %>% 
    separate(feature, into = c("mz", "time"), sep = "_", remove = FALSE) %>% 
    dplyr::select(mz, time, everything()) %>% 
    mutate_at("mz", as.numeric) %>% 
    mutate_at("time", as.numeric) %>% 
    mutate(index.A = seq(mz)) 

dataB <- hilic_neg_lib_2 %>% 
    dplyr::select(mz, time, everything()) %>% 
    mutate_at("mz", as.numeric) %>% 
    mutate_at("time", as.numeric) %>% 
    mutate(index.B = seq(mz)) 

overlap <- find.Overlapping.mzs(dataA, dataB, mz.thresh = 10, time.thresh = 10,
                                alignment.tool=NA)

lib_overlap <- dataB[overlap$index.B,]
data_overlap <- dataA[overlap$index.A,]

data_lib <- cbind(data_overlap, lib_overlap, overlap)

data_lib %>% write_tsv("results/NDEV/level1-5_hilic_neg_features_TotalEV.txt")
