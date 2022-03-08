library(tidyverse)

efo <- read.csv("data/full_disease_list.csv")
map <- read.csv("data/UK_Biobank_master_file.tsv", sep = "\t")

mapped <- efo %>% left_join(map, by = c("specificDiseaseId" = "MAPPED_TERM_URI"))
tmp <- mapped %>% select(diseaseId, diseaseName, specificDiseaseId, specificDiseaseName, description, ZOOMA.QUERY, ICD10_CODE.SELF_REPORTED_TRAIT_FIELD_CODE, MAPPED_TERM_LABEL, MAPPING_TYPE ) 

ukb <- read.csv("data/ukb_multimorbidity.csv")

mm <- tmp %>% left_join(ukb, by = c("ICD10_CODE.SELF_REPORTED_TRAIT_FIELD_CODE" = "Disease1"))
tmp2 <- mm %>% filter(!is.na(Disease2))
