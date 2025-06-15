#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Chemical detection rates
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(tidyverse)
library(xlsx)
source("code/functions.R")

# call in the files
chem_key <- read_tsv("exposomics/clean_data/chemical_key.txt")
chem_levels <- read_tsv("exposomics/clean_data/chemical_data.txt")
chem_LOD <- read_tsv("exposomics/clean_data/chemical_lod.txt")
id_key <- read_tsv("exposomics/clean_data/gc_ids.txt") 

# Recode entries <LOD
# 1: NF/ NA
# 2: <LOD
# 3: detected

chem_recode <- chem_levels %>% 
    select(-ID) %>% 
    map_df(~codeLODNF(.x)) %>% 
    mutate(ID = chem_levels$ID) %>% 
    select(ID, everything())

# Calculate detection count by chemical
chem_detect <- chem_recode %>% 
    separate(ID, into = c("matrix", "person"), sep = "_", remove = FALSE) %>% 
    filter(matrix != "Depleted") %>% 
    gather(key = "index", value = "code", - ID, -matrix, -person) %>% 
    group_by(index) %>% 
    summarise(num.detect = sum(code)) %>% 
    ungroup()

# Calculate detection count by chemical and matrix
chem_matrix_detect <- chem_recode %>% 
    separate(ID, into = c("matrix", "person"), sep = "_", remove = FALSE) %>% 
    filter(matrix != "Depleted") %>% 
    gather(key = "index", value = "code", - ID, -matrix, -person) %>% 
    group_by(index, matrix) %>% 
    summarise(num.detect = sum(code)) %>% 
    ungroup()

chem_matrix_detect_namess <- merge(chem_matrix_detect, chem_key, by = "index" )

# Calculate detection count by individual
chem_person_detect <- chem_recode %>% 
    separate(ID, into = c("matrix", "person"), sep = "_", remove = FALSE) %>% 
    filter(matrix != "Depleted") %>% 
    gather(key = "index", value = "code", - ID, -matrix, -person) %>% 
    group_by(index, person) %>% 
    summarise(num.detect = sum(code)) %>% 
    ungroup()

# Chemicals not detected in any matrix
chemical_none <- chem_detect %>% 
    filter(num.detect == 0) %>% 
    left_join(., chem_key, by = "index")

# Chemicals detected in all four matrices in a person
chemical_all <- chem_person_detect %>% 
    filter(num.detect == 4) %>% 
    left_join(., chem_key, by = "index")

chemical_all_unique <- chem_person_detect %>% 
    filter(num.detect == 4) %>% 
    group_by(index) %>% 
    summarise(num.person = n()) %>% 
    ungroup() %>% 
    left_join(., chem_key, by = "index")

# Chemicals detected in four matrices in a person, no depleted fraction
chem_person_detect_no_depl <- chem_recode %>% 
    separate(ID, into = c("matrix", "person"), sep = "_", remove = FALSE) %>% 
    filter(matrix != "Depleted") %>% 
    gather(key = "index", value = "code", - ID, -matrix, -person) %>% 
    group_by(index, person) %>% 
    summarise(num.detect = sum(code)) %>% 
    ungroup()

chemical_all_no_depl <- chem_person_detect_no_depl %>% 
    filter(num.detect > 2) %>% 
    left_join(., chem_key, by = "index")

chemical_all_unique_no_depl <- chem_person_detect_no_depl %>% 
    filter(num.detect > 2) %>% 
    group_by(index) %>% 
    summarise(num.person = n()) %>% 
    ungroup() %>% 
    left_join(., chem_key, by = "index")

# Detection count in tissue and NDEV
chem_neuro_detect <- chem_recode %>% 
    separate(ID, into = c("matrix", "person"), sep = "_", remove = FALSE) %>% 
    gather(key = "index", value = "code", - ID, -matrix, -person) %>% 
    filter(matrix %in% c("Tissue", "NDEV")) %>% 
    group_by(index, person) %>% 
    summarise(num.detect = sum(code)) %>% 
    ungroup()

# Chemicals detected in neuro matrices in a person
chemical_neuro <- chem_neuro_detect %>% 
    filter(num.detect == 2) %>% 
    left_join(., chem_key, by = "index")

chemical_neuro_unique <- chem_neuro_detect %>% 
    filter(num.detect == 2)  %>% 
    group_by(index) %>% 
    summarise(num.person = n()) %>% 
    ungroup() %>% 
    left_join(., chem_key, by = "index")

# Chemicals detected in tissue
chem_tissue_detect <- chem_recode %>% 
    separate(ID, into = c("matrix", "person"), sep = "_", remove = FALSE) %>% 
    gather(key = "index", value = "code", - ID, -matrix, -person) %>% 
    filter(matrix %in% c("Tissue")) %>% 
    group_by(index, person) %>% 
    summarise(num.detect = sum(code)) %>% 
    ungroup()

# Chemicals detected in tissue in a person
chemical_tissue <- chem_tissue_detect %>% 
    filter(num.detect == 1) %>% 
    left_join(., chem_key, by = "index")

chemical_tissue_unique <- chem_tissue_detect %>% 
    filter(num.detect == 1)  %>% 
    group_by(index) %>% 
    summarise(num.person = n()) %>% 
    ungroup() %>% 
    left_join(., chem_key, by = "index")

# Chemicals detected in NDEV
chem_ndev_detect <- chem_recode %>% 
    separate(ID, into = c("matrix", "person"), sep = "_", remove = FALSE) %>% 
    gather(key = "index", value = "code", - ID, -matrix, -person) %>% 
    filter(matrix %in% c("NDEV")) %>% 
    group_by(index, person) %>% 
    summarise(num.detect = sum(code)) %>% 
    ungroup()

# Chemicals detected in NDEVs in a person
chemical_ndev <- chem_ndev_detect %>% 
    filter(num.detect == 1) %>% 
    left_join(., chem_key, by = "index")

chemical_ndev_unique <- chem_ndev_detect %>% 
    filter(num.detect == 1)  %>% 
    group_by(index) %>% 
    summarise(num.person = n()) %>% 
    ungroup() %>% 
    left_join(., chem_key, by = "index")

# Chemicals detected in depleted fraction
chem_deplete_detect <- chem_recode %>% 
    separate(ID, into = c("matrix", "person"), sep = "_", remove = FALSE) %>% 
    gather(key = "index", value = "code", - ID, -matrix, -person) %>% 
    filter(matrix %in% c("Depleted")) %>% 
    group_by(index, person) %>% 
    summarise(num.detect = sum(code)) %>% 
    ungroup()

# Chemicals detected in depleted fraction in a person
chemical_deplete <- chem_deplete_detect %>% 
    filter(num.detect == 1) %>% 
    left_join(., chem_key, by = "index")

chemical_deplete_unique <- chem_deplete_detect %>% 
    filter(num.detect == 1)  %>% 
    group_by(index) %>% 
    summarise(num.person = n()) %>% 
    ungroup() %>% 
    left_join(., chem_key, by = "index")

# Chemicals detected in total EVs
chem_tev_detect <- chem_recode %>% 
    separate(ID, into = c("matrix", "person"), sep = "_", remove = FALSE) %>% 
    gather(key = "index", value = "code", - ID, -matrix, -person) %>% 
    filter(matrix %in% c("TotalEV")) %>% 
    group_by(index, person) %>% 
    summarise(num.detect = sum(code)) %>% 
    ungroup()

# Chemicals detected in total EVs in a person
chemical_tev <- chem_tev_detect %>% 
    filter(num.detect == 1) %>% 
    left_join(., chem_key, by = "index")

chemical_tev_unique <- chem_tev_detect %>% 
    filter(num.detect == 1)  %>% 
    group_by(index) %>% 
    summarise(num.person = n()) %>% 
    ungroup() %>% 
    left_join(., chem_key, by = "index")

# Chemicals detected in Serum
chem_serum_detect <- chem_recode %>% 
    separate(ID, into = c("matrix", "person"), sep = "_", remove = FALSE) %>% 
    gather(key = "index", value = "code", - ID, -matrix, -person) %>% 
    filter(matrix %in% c("Serum")) %>% 
    group_by(index, person) %>% 
    summarise(num.detect = sum(code)) %>% 
    ungroup()

# Chemicals detected in Serum in a person
chemical_serum <- chem_serum_detect %>% 
    filter(num.detect == 1) %>% 
    left_join(., chem_key, by = "index")

chemical_serum_unique <- chem_serum_detect %>% 
    filter(num.detect == 1)  %>% 
    group_by(index) %>% 
    summarise(num.person = n()) %>% 
    ungroup() %>% 
    left_join(., chem_key, by = "index")

# Chemicals detected in any matrix
chemical_allM <- chem_detect %>% 
    #filter(num.detect == 1) %>% 
    left_join(., chem_key, by = "index")

chemical_allM_unique <- chem_detect %>% 
    filter(num.detect > 0)  %>% 
    group_by(index) %>% 
    summarise(num.person = n()) %>% 
    ungroup() %>% 
    left_join(., chem_key, by = "index")

# Saving data in excel
dataframe1 = chemical_none %>%  # Chemicals not detected in any matrix (n = 12) 
    select(-index) %>% 
    select(Chemical, `Number of samples` = num.detect) %>% 
    as.data.frame()

dataframe2 = chemical_all_unique %>%  # Chemicals detected in all 5 matrices (n = 20)
    select(-index) %>% 
    select(Chemical, `Number of individuals` = num.person) %>% 
    as.data.frame()

dataframe3 = chemical_neuro_unique %>% # Chemicals detected in the two neuro matrices: brain and NDEVs (n = 32)
    select(-index) %>% 
    select(Chemical, `Number of individuals` = num.person) %>% 
    as.data.frame()

dataframe4 = chemical_tissue_unique %>% # Chemicals detected in brain tissue (n = 60)
    select(-index) %>% 
    select(Chemical, `Number of individuals` = num.person) %>% 
    as.data.frame()

dataframe5 = chemical_ndev_unique %>% # Chemicals detected in NDEVs (n = 64)
    select(-index) %>% 
    select(Chemical, `Number of individuals` = num.person) %>% 
    as.data.frame()

dataframe6 = chemical_deplete_unique %>% # Chemicals detected in NDEV depleted fraction (n = 63)
    select(-index) %>% 
    select(Chemical, `Number of individuals` = num.person) %>% 
    as.data.frame()

dataframe7 = chemical_tev_unique %>% # Chemicals detected in total EVs (n = 74)
    select(-index) %>% 
    select(Chemical, `Number of individuals` = num.person) %>% 
    as.data.frame()

dataframe8 = chemical_serum_unique %>% # Chemicals detected in serum (n = 65)
    select(-index) %>% 
    select(Chemical, `Number of individuals` = num.person) %>% 
    as.data.frame()

dataframe9 = chemical_allM %>% # Chemicals detected in any matrix (n = 109)
    select(-index) %>% 
    filter(num.detect > 0) %>% 
    select(Chemical, `Number of samples` = num.detect) %>% 
    as.data.frame()

dataframe10 = chemical_all_unique_no_depl %>%  # Chemicals detected 4 matrices (n = 27, no depleted)
    select(-index) %>% 
    select(Chemical, `Number of individuals` = num.person) %>% 
    as.data.frame()

dataframe11 = chem_LOD %>% 
    select(-index) %>% 
    as.data.frame()

write.xlsx(dataframe11, file="exposomics/results/chemicals_detection.xlsx", sheetName="chemicals", row.names = FALSE)
write.xlsx(dataframe9, file="exposomics/results/chemicals_detection.xlsx", sheetName="chemicals_in_any", append=TRUE, row.names=FALSE)
write.xlsx(dataframe1, file="exposomics/results/chemicals_detection.xlsx", sheetName="chemicals_in_none", append=TRUE, row.names=FALSE)
write.xlsx(dataframe2, file="exposomics/results/chemicals_detection.xlsx", sheetName="chemicals_in_all", append=TRUE, row.names=FALSE)
write.xlsx(dataframe3, file="exposomics/results/chemicals_detection.xlsx", sheetName="chemicals_in_neuro", append=TRUE, row.names=FALSE)
write.xlsx(dataframe4, file="exposomics/results/chemicals_detection.xlsx", sheetName="chemicals_in_brain", append=TRUE, row.names=FALSE)
write.xlsx(dataframe5, file="exposomics/results/chemicals_detection.xlsx", sheetName="chemicals_in_ndev", append=TRUE, row.names=FALSE)
write.xlsx(dataframe6, file="exposomics/results/chemicals_detection.xlsx", sheetName="chemicals_in_depleted", append=TRUE, row.names=FALSE)
write.xlsx(dataframe7, file="exposomics/results/chemicals_detection.xlsx", sheetName="chemicals_in_totalEV", append=TRUE, row.names=FALSE)
write.xlsx(dataframe8, file="exposomics/results/chemicals_detection.xlsx", sheetName="chemicals_in_serum", append=TRUE, row.names=FALSE)
write.xlsx(dataframe10, file="exposomics/results/chemicals_detection.xlsx", sheetName="chemicals_in_4_mat_updated", append=TRUE, row.names=FALSE)

# Plots
chemical_tissue_unique %>% 
    ggplot(aes(x = num.person, y = reorder(Chemical, num.person))) +
    geom_bar(width = 0.6, stat = "identity") +
    theme_bw() +
    labs(x = "Count", y = "")

chemical_ndev_unique %>% 
    ggplot(aes(x = num.person, y = reorder(Chemical, num.person))) +
    geom_bar(width = 0.6, stat = "identity") +
    theme_bw() +
    labs(x = "Count", y = "")

chemical_deplete_unique %>% 
    ggplot(aes(x = num.person, y = reorder(Chemical, num.person))) +
    geom_bar(width = 0.6, stat = "identity") +
    theme_bw() +
    labs(x = "Count", y = "")

chemical_tev_unique %>% 
    ggplot(aes(x = num.person, y = reorder(Chemical, num.person))) +
    geom_bar(width = 0.6, stat = "identity") +
    theme_bw() +
    labs(x = "Count", y = "")

chemical_serum_unique %>% 
    ggplot(aes(x = num.person, y = reorder(Chemical, num.person))) +
    geom_bar(width = 0.6, stat = "identity") +
    theme_bw() +
    labs(x = "Count", y = "")

chemical_allM %>% 
    ggplot(aes(x = num.detect, y = reorder(Chemical, num.detect))) +
    geom_bar(width = 0.6, stat = "identity") +
    theme_bw() +
    labs(x = "Count", y = "")
