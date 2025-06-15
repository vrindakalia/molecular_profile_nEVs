#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Chemical characteristics and correlations
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(tidyverse)
library(openxlsx)

all_brain_cor <- read_tsv("exposomics/results/brain_matrix_spearman_corr.txt")
chem_char <- read.xlsx("exposomics/CCD-Batch-Search_2024-06-24_03_27_28.xlsx", sheet = "Main Data")

cor_chem_char <- merge(all_brain_cor, chem_char, by.x = "Chemical", by.y = "INPUT") %>% 
    mutate(LOGP.med = case_when(OCTANOL_WATER_PARTITION_LOGP_OPERA_PRED <= median(OCTANOL_WATER_PARTITION_LOGP_OPERA_PRED, na.rm = TRUE) ~ 0,
                                OCTANOL_WATER_PARTITION_LOGP_OPERA_PRED > median(OCTANOL_WATER_PARTITION_LOGP_OPERA_PRED, na.rm = TRUE) ~ 1))

cor_chem_char_sub <- cor_chem_char %>% 
    select(Chemical, contains("cor"), contains("LOGP"), AVERAGE_MASS, BIODEGRADATION_HALF_LIFE_DAYS_DAYS_OPERA_PRED, OPERA_PKAA_OPERA_PRED, `DENSITY_G/CM^3_TEST_PRED`)

cor_chem_char_sub %>% 
    filter(!is.na(LOGP.med)) %>% 
    ggplot(aes(x = factor(LOGP.med), y = ev.cor)) +
    geom_boxplot() +
    geom_smooth(method = "loess")

cor.test(cor_chem_char_sub$ndev.cor, cor_chem_char_sub$OCTANOL_WATER_PARTITION_LOGP_OPERA_PRED, method = "spear")
cor.test(cor_chem_char_sub$serum.cor, cor_chem_char_sub$OCTANOL_WATER_PARTITION_LOGP_OPERA_PRED, method = "spear")
cor.test(cor_chem_char_sub$ev.cor, cor_chem_char_sub$OCTANOL_WATER_PARTITION_LOGP_OPERA_PRED, method = "spear")

cor.test(cor_chem_char_sub$ndev.cor, cor_chem_char_sub$BIODEGRADATION_HALF_LIFE_DAYS_DAYS_OPERA_PRED, method = "spear")
cor.test(cor_chem_char_sub$serum.cor, cor_chem_char_sub$BIODEGRADATION_HALF_LIFE_DAYS_DAYS_OPERA_PRED, method = "spear")
cor.test(cor_chem_char_sub$ev.cor, cor_chem_char_sub$BIODEGRADATION_HALF_LIFE_DAYS_DAYS_OPERA_PRED, method = "spear")

cor.test(cor_chem_char_sub$ndev.cor, cor_chem_char_sub$AVERAGE_MASS, method = "spear")
cor.test(cor_chem_char_sub$serum.cor, cor_chem_char_sub$AVERAGE_MASS, method = "spear")
cor.test(cor_chem_char_sub$ev.cor, cor_chem_char_sub$AVERAGE_MASS, method = "spear")

cor.test(cor_chem_char_sub$ndev.cor, cor_chem_char_sub$`DENSITY_G/CM^3_TEST_PRED`, method = "spear")
cor.test(cor_chem_char_sub$serum.cor, cor_chem_char_sub$`DENSITY_G/CM^3_TEST_PRED`, method = "spear")
cor.test(cor_chem_char_sub$ev.cor, cor_chem_char_sub$`DENSITY_G/CM^3_TEST_PRED`, method = "spear")

cor_chem_char_sub %>% 
    ggplot(aes(x = `DENSITY_G/CM^3_TEST_PRED`, y = serum.cor)) +
    geom_point() +
    geom_smooth(method = "lm")
