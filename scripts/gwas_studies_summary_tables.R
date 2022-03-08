library(tidyverse)
library(googlesheets4)

save_dir = "data/analysis/"
load(paste0(save_dir, "ard_leads_filtered.Rda"))
load(paste0(save_dir, "individual_disease_n_communities.Rda"))


## Get GWAS study tables
## N GWAS studies, sample size, n variants, n genetic signals per trait and per summarised phenotype
ard_studies <- ard_leads %>% 
  select(studyId, 
         diseaseName, 
         n_cases, 
         n_initial, 
         ancestry_initial, 
         ancestry_replication)
study_size <- ard_studies %>% 
  group_by(diseaseName) %>% 
  bind_rows(ard_studies %>% 
              select(-diseaseName) %>% 
              unique() %>% 
              mutate(diseaseName = "Total")) %>% 
  summarise_at(.vars = "n_initial", 
               .funs = c("min", "max", "median"))

total_gwas_tbl <- ard_leads %>% 
  select(studyId, lead_variantId, has_sumstats) %>%
  unique() %>%
  group_by(studyId) %>% 
  summarise(has_sumstats = max(has_sumstats), 
            n_variants = length(studyId)) %>% 
  ungroup() %>% 
  summarise(n_gwas_studies = length(studyId), 
            n_variants = sum(n_variants), 
            has_sumstats = sum(has_sumstats)) %>%
  mutate(diseaseName = "Total")

gwas_tbl <- 
  ard_leads %>% 
  group_by(diseaseId, 
           diseaseName, 
           specificDiseaseId, 
           specificDiseaseName, 
           studyId)  %>% 
  summarise(has_sumstats = max(has_sumstats), 
            n_variants = length(studyId)) %>% 
  ungroup() %>% 
  group_by(diseaseId, diseaseName, specificDiseaseId, specificDiseaseName) %>% 
  summarise(n_gwas_studies = length(specificDiseaseName), 
            has_sumstats = max(has_sumstats), 
            n_variants = sum(n_variants))

summarised_gwas_tbl <- 
  gwas_tbl %>% 
  ungroup() %>% 
  group_by(diseaseName) %>% 
  summarise(n_gwas_studies = sum(n_gwas_studies), 
            has_sumstats = sum(has_sumstats), 
            n_variants = sum(n_variants)) %>% 
  bind_rows(total_gwas_tbl) %>% 
  mutate(has_sumstats = paste0(has_sumstats, 
                               " (", 
                               round((has_sumstats/n_gwas_studies)*100), 
                               "%)")) %>% 
  left_join(study_size) %>% 
  full_join(ard_leads %>% 
              select(diseaseName, morbidity) %>% 
              unique()) %>%
  full_join(individual_communities) %>%
  select(-morbidity) %>%
  select(diseaseName, n_gwas_studies, has_sumstats, min, max, median, n_variants, everything())

full_gwas_tbl <- gwas_tbl %>% 
  full_join(diseases) %>% 
  select(-therapeuticAreas) %>%
  filter(!is.na(n_gwas_studies) | diseaseId == specificDiseaseId) %>%
  arrange(diseaseName) 

gs4_create("gwas_table", sheets = summarised_gwas_tbl)
gs4_create("gwas_specific_table", sheets = full_gwas_tbl)

ancestries <- ard_leads %>% select(studyId, ancestry_initial) %>% unique() %>% unnest(ancestry_initial) %>% separate(ancestry_initial, into = c("ancestry", "n"), sep = "=") %>% group_by(studyId) %>% mutate(n = as.numeric(n), total = sum(n), proportion = n/total)

majority_population <- ancestries %>% group_by(studyId) %>% top_n(1, proportion) %>%
  mutate(single_population = ifelse(proportion == 1, TRUE, FALSE)) 

majority_population %>% ungroup() %>% count(single_population, ancestry) %>% mutate(p = n/sum(n)*100) %>% arrange(-p)
