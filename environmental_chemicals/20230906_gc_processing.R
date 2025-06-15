#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Targeted GC data
# Replace LOD
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(tidyverse)
library(xlsx)
source("code/functions.R")

# call in the files
chem_key <- read_tsv("exposomics/clean_data/chemical_key.txt")
chem_levels <- read_tsv("exposomics/clean_data/chemical_data.txt")
chem_LOD <- read_tsv("exposomics/clean_data/chemical_lod.txt")
id_key <- read_tsv("exposomics/clean_data/gc_ids.txt") 

chemical_all <- read.xlsx(file="exposomics/results/chemicals_detection.xlsx", sheetName="chemicals_in_4_mat_updated")

chem_all_filter <- chemical_all %>% 
    merge(., chem_key, by = "Chemical") %>% 
    filter(Number.of.individuals > 2)

# Restrict data to those chemicals detected in all matrices in at least 3 people
chem_levels_filter <- chem_levels %>% 
    select(chem_all_filter$index) %>% 
    mutate_at("X15", as.character) %>% 
    replace_na(list(X15 = "<LOD")) # need to check with Kate why value is NA, treated as <LOD for initial analysis

# Find the minimum value detected for each chemical
chem_Dmin <- chem_levels_filter %>% 
    map_df(~as.character(.x)) %>% 
    map_df(~replaceLODzero(.x)) %>%
    map_df(~as.numeric(.x)) %>% 
    map_df(~FindMin(.x)) %>% 
    gather(key = "index", value = "min.detection")

chemicals_summary <- chem_LOD %>% 
    merge(chem_Dmin,., by = "index") %>% 
    mutate(replace.use = ifelse(Serum_LOD > min.detection, "min.detection", "serum_LOD")) %>% 
    mutate(replace.value = case_when(replace.use == "serum_LOD" ~ Serum_LOD/2,
                                     replace.use == "min.detection" ~ min.detection/2)) %>% 
    mutate_at("replace.value", as.character)

# Replace <LOD 
chemicals_summary$index <- chemicals_summary$index[match(names(chem_levels_filter), chemicals_summary$index)]

for(j in 1:dim(chem_levels_filter)[2]){
    for(i in 1:dim(chem_levels_filter)[1]){
        if(chem_levels_filter[i,j] == "<LOD") chem_levels_filter[i,j] = chemicals_summary$replace.value[j]
        else chem_levels_filter[i,j] = chem_levels_filter[i,j]
    }
}

chem_levels_filter$ID <- chem_levels$ID

chem_levels_filter %>% 
    select(ID, everything()) %>% 
    write_tsv("exposomics/clean_data/in_four_replaced_atl3_updated.txt")

