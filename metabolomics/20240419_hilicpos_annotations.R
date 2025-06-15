#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Annotations for HILIC + column
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(tidyverse)
library(xMSanalyzer)
library(openxlsx)

# Call in the metabolomic data
hilpos <- read_tsv("HILIC_column/RAW_mzcalibrated_untargeted_featuretable_HILICpos.txt") %>%
    unite("mz_time", c(mz, time), sep = "_", remove = FALSE) %>% 
    dplyr::select(mz, time, everything())
#dplyr::select(-`id`) %>% 
#t() %>% 
#as.data.frame() %>% 
#janitor::row_to_names(1) %>% 
#rownames_to_column(var = "feature") %>% 
#separate(feature, into = c("mz", "time"), sep = "_") %>% 
#dplyr::select(mz, time, everything()) %>% 
#mutate_at("mz", as.numeric) %>% 
#mutate_at("time", as.numeric)

hilpos_meanint <- hilpos %>% 
    dplyr::select(contains("_"), -median_CV) %>% 
    column_to_rownames(var = "mz_time") %>% 
    mutate(mean = rowMeans(.)) %>% 
    rownames_to_column(var = "mz_time") %>% 
    dplyr::select(mz_time, mean) 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#              LEVEL 1 Confidence           #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# HILIC pos
lib_pos <- read.xlsx("/Users/vk2316/Documents/Baccarelli_lab/EVs/with_DW/HRE_Lab/library/230411-HILICpos_Merged_InchKey_Filtered_Confirmed_ChemicalDB.xlsx") %>% 
    dplyr::rename(mz.std = Confirmed.mz, time.std = Confirmed.rt) %>% 
    dplyr::select(mz.std, time.std, everything()) %>% 
    dplyr::rename(Name = Metabolite_Name)

# Find overlapping features: positive
dataA = hilpos
dataB = lib_pos
pos_overlap <- find.Overlapping.mzs(dataA, dataB, mz.thresh = 20, time.thresh = 30, alignment.tool = "apLCMS")

std_pos_ann <- lib_pos %>% 
    slice(pos_overlap$index.B) %>%
    mutate(confidence = "Level 1") %>% 
    unite("mz_time_std", c(mz.std, time.std), sep = "_", remove = FALSE)

std_pos_ann_filter_mz <- pos_overlap %>% 
    unite("mz_time_std", c(mz.data.B, time.data.B), sep = "_", remove = FALSE) %>% 
    merge(., std_pos_ann, by = c("mz_time_std")) %>% 
    dplyr::rename(mz = mz.data.A, time = time.data.A) %>% 
    mutate(delta.mz = abs(((mz - mz.std)/mz.std)*10^6)) %>% 
    dplyr::select(mz, time, mz_time_std, delta.mz, time.difference, Name, confidence, PubChem_CID, `Chemical Formula` = MolecularFormula) %>% 
    unite("mz_time", c(mz, time), sep = "_", remove = FALSE) %>%
    merge(., hilpos_meanint, by = "mz_time") %>% 
    group_by(mz_time) %>% 
    mutate(min.delta.mz = min(delta.mz)) %>% 
    ungroup() %>% 
    mutate(filter.mz = case_when(delta.mz == min.delta.mz ~ 0,
                                 delta.mz > min.delta.mz ~ 1)) %>% 
    filter(filter.mz == 0) 

std_pos_ann_filter_mz_time <- std_pos_ann_filter_mz %>% 
    group_by(mz_time) %>% 
    mutate(min.delta.time = min(time.difference)) %>% 
    ungroup() %>% 
    mutate(filter.time = case_when(time.difference == min.delta.time ~ 0,
                                   time.difference > min.delta.time ~ 1)) %>% 
    filter(filter.time == 0) %>% 
    group_by(mz_time_std) %>% 
    mutate(min.delta.time = min(time.difference)) %>% 
    ungroup() %>% 
    mutate(filter.time.std = case_when(time.difference == min.delta.time ~ 0,
                                       time.difference > min.delta.time ~ 1)) %>% 
    filter(filter.time.std == 0) 

std_pos_ann_filter_mz_time_duplicate <- std_pos_ann_filter_mz_time %>% 
    group_by(mz_time) %>% 
    mutate(rep = seq_along(mz_time)) %>% 
    ungroup() %>% 
    filter(rep > 1)

std_pos_ann_filter_mz_time_unique <- std_pos_ann_filter_mz_time %>% 
    distinct(mz_time, .keep_all = TRUE)  %>% 
    mutate(Annotation.confidence.score = NA)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#              LEVEL 3 or greater           #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# HILIC POS
ann_pos <- read_csv("xMSannotator/HILIC POS/Stage5.csv")

ann_pos_isotopes <- ann_pos %>% 
    filter(is.na(delta_ppm)) %>% 
    distinct(., mz, .keep_all = TRUE) %>% 
    dplyr::select(mz, time, chemical_ID, Annotation.confidence.score, delta_ppm, Name, Formula, Adduct,mz.theor = theoretical.mz) 

ann_pos_sig <- ann_pos %>% 
    filter(!is.na(delta_ppm)) %>% 
    group_by(mz) %>% 
    mutate(min.delta = min(delta_ppm, na.rm = TRUE)) %>% 
    ungroup() %>% 
    filter(delta_ppm <= min.delta) %>% 
    group_by(mz) %>% 
    mutate(num = seq_along(Name)) %>% 
    mutate(tot = max(num)) %>% 
    ungroup() %>% 
    filter(!mz %in% std_pos_ann_filter_mz_time_unique$mz)

ann_pos_sig_tot1 <- ann_pos_sig %>% 
    filter(tot == 1) %>% 
    dplyr::select(mz, time, chemical_ID, Annotation.confidence.score, delta_ppm, Name, Formula, Adduct,mz.theor = theoretical.mz) 

ann_pos_sig_totgt1 <- ann_pos_sig %>% 
    filter(tot > 1) %>% 
    mutate(name.semi = paste0(Name, ";")) %>% 
    unite("mz_time", c(mz, time), sep = "_", remove = FALSE)

ann_pos_sig_totgt1_distinct <- ann_pos_sig_totgt1 %>%
    unite("mz_time", c(mz, time), sep = "_", remove = FALSE) %>% 
    distinct(., mz_time, .keep_all = TRUE)

names.concat = as.data.frame(matrix(NA, nrow = dim(ann_pos_sig_totgt1_distinct)[1], ncol = 3))
names(names.concat) = c("mz", "time", "name.concat")

for(i in 1:dim(ann_pos_sig_totgt1_distinct)[1]){
    
    ann_pos_sub = ann_pos_sig_totgt1 %>% 
        filter(mz_time == ann_pos_sig_totgt1_distinct$mz_time[i]) 
    
    mz_name = paste(ann_pos_sub$name.semi, collapse = "")
    
    names.concat$mz[i] = ann_pos_sig_totgt1_distinct$mz[i]
    names.concat$time[i] = ann_pos_sig_totgt1_distinct$time[i]
    names.concat$name.concat[i] = mz_name
}

ann_pos_sig_totgt1_names <- ann_pos_sig_totgt1 %>% 
    left_join(., names.concat, by = c("mz", "time"))

ann_pos_sig_clean <- ann_pos_sig_totgt1_names %>% 
    distinct(mz_time, .keep_all = TRUE) %>% 
    dplyr::select(mz, time, chemical_ID, Annotation.confidence.score, delta_ppm, Name = name.concat, Formula, Adduct, mz.theor = theoretical.mz) %>% 
    rbind(ann_pos_sig_tot1, ann_pos_isotopes, .) %>% 
    mutate(confidence = case_when(Annotation.confidence.score > 1 ~ "Level 3",
                                  Annotation.confidence.score <= 1 ~ "Level 5")) %>% 
    dplyr::rename(HMDB_ID = chemical_ID, `Chemical Formula` = Formula, delta.mz = delta_ppm) %>% 
    mutate(time.theor = time) %>% 
    dplyr::select(mz, time, mz.theor, time.theor, Name, `Chemical Formula`, HMDB_ID, confidence, delta.mz, Adduct, Annotation.confidence.score) %>% 
    unite("mz_time", c(mz, time), sep = "_", remove = FALSE) 


# Combining annotation data from all confidence levels
names(ann_pos_sig_clean)
names(std_pos_ann_filter_mz_time_unique)
ann_pos_sig_clean <- ann_pos_sig_clean %>% 
    dplyr::rename(ID = HMDB_ID)

pos_ann_level135 <- std_pos_ann_filter_mz_time_unique %>% 
    mutate(Adduct = NA) %>% 
    separate(mz_time_std, into = c("mz.theor", "time.theor"), sep = "_") %>% 
    mutate(ID = paste0("PubChem_", PubChem_CID)) %>% 
    dplyr::select(mz, time, mz_time, mz.theor, delta.mz, time.theor, Name, `Chemical Formula`, ID, confidence, Adduct, Annotation.confidence.score) %>% 
    rbind(., ann_pos_sig_clean)

pos_ann_level135 %>% write_tsv("results/hilpos_annotations_level1-5.txt")

