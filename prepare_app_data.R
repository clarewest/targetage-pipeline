library(tidyverse)
load("data/analysis/target_annotations.Rda")             # Target annotations

### Functions
add_go_hallmarks <- function(df){

  ## Get hallmarks for each target
  go_terms <- target_annotations %>% 
    select(targetId, ageingHallmarks) %>% 
    unnest(cols = c(ageingHallmarks)) %>%
    left_join(go_icons, by = "goHallmark") %>%
    select(targetId, goIcon, goHallmark, goHallmarkTerm) %>%
    unique() %>% 
    group_by(targetId) %>%
    summarise(
      n_hallmarks = length(unique(goHallmark)),
      goHallmarks = toString(unique(goHallmark)),
      goHallmarkTerms = toString(unique(goHallmarkTerm)),
      icon.goHallmarks = paste(unique(goIcon), collapse = " "),
      .groups = "drop_last"
    ) %>%
    mutate(icon.goHallmarks = str_replace_all(
      icon.goHallmarks,
      "></i>",
      paste0(
        " data-toggle=\"tooltip\" data-placement=\"right\" title=\"",
        goHallmarks,
        "\"></i>"
      )
    )) 
  
  df %>% 
    left_join(go_terms, by = "targetId") %>%
    unite(col = icon.goHallmarks, icon.goHallmarks, CellAge.senescence_icon, sep = "<br>", na.rm = TRUE)
  
}
  

add_genage_icons <- function(df){
  gen <- df %>% 
    select(targetSymbol, GenAge) %>% 
    unnest(cols = c(GenAge)) %>% 
    left_join(model_imgs, by = "GenAge.model.organism") %>%
    mutate(img = ifelse(is.na(model.organism.img.path), 
                        NA, 
                        paste0("<img src=\"", model.organism.img.path, "\" height=\"20\" data-toggle=\"tooltip\" data-placement=\"right\" title=\"", GenAge.longevity.influence, " in ", GenAge.model.organism, "\"></img>"))) %>%
    left_join(longevity_imgs, by = "GenAge.longevity.influence") %>%
    mutate(GenAge.longevity.icon = ifelse(is.na(GenAge.longevity.influence), 
                                          NA, 
                                          paste0("<i class=\"fa fa-", GenAge.longevity.iconcode,"\" data-toggle=\"tooltip\" data-placement=\"right\" title=\"", GenAge.longevity.influence, " in ", GenAge.model.organism, "\"></i>"))) %>% 
    unite(col = GenAge.longevity.icon, c(GenAge.longevity.icon, img), na.rm = TRUE, sep = "") %>%
    select(targetSymbol, GenAge.ID, GenAge.longevity.icon) %>% 
    add_genage_links() %>%
    ungroup()
  df %>% select(-GenAge) %>% 
    left_join(gen)
}

add_genage_links <- function(df){
  df %>%
    group_by(targetSymbol, GenAge.ID) %>% 
    summarise(icon.GenAge = paste(GenAge.longevity.icon, collapse = "<br>")) %>% 
    mutate(GenAge.URL = paste0("https://genomics.senescence.info/genes/entry.php?hgnc=", targetSymbol),
           icon.GenAge = ifelse(icon.GenAge=="NA", NA, icon.GenAge),
           GenAge.link = paste0("<a target=\"_blank\" href=\"", GenAge.URL, "\">", GenAge.ID, "</a>")) %>%
    unite(col = icon.GenAge, c(GenAge.link, icon.GenAge), sep = "<br>", na.rm = TRUE)
}

add_cellage_icons <- function(df){
  cell <- df %>% 
    select(targetId, CellAge) %>% 
    unnest(cols = c(CellAge)) %>% 
    select(targetId, CellAge.organism, CellAge.senescence_effect) %>% 
    left_join(data.frame(
      CellAge.senescence_effect = c("Induces", "Inhibits", "Unclear", NA),
      senescence_effect_iconcode = c("arrow-up", "arrow-down", "question", NA)),
      by = "CellAge.senescence_effect") %>% 
    mutate(CellAge.senescence_icon = paste0("<div data-toggle=\"tooltip\" data-placement=\"right\" title=\"", CellAge.senescence_effect, " Senescence\">",
                                            "<i class=\"fa fa-", senescence_effect_iconcode,"\"</i>",
                                            "<i class=\"fab fa-canadian-maple-leaf\"</i></div>")) %>%
    select(-senescence_effect_iconcode, -CellAge.organism) 
  df %>% left_join(cell)
}

add_details_column <- function(df){
  description <- df %>% select(targetId, functionDescriptions) %>%
    unnest(functionDescriptions) %>% 
    group_by(targetId) %>% 
    summarise(details.description = paste(functionDescriptions, collapse = " ")) %>% 
    ungroup() %>% select(targetId, details.description)
  
  
  df %>% left_join(description) %>% 
    rowwise() %>% 
    mutate(details_output = data.frame(
      "ID" = targetId,
      "Gene symbol" = targetSymbol,
      "Gene name" = targetName,
      "Description" = details.description,
      "Biotype" = biotype, 
      "PDB structures" = PDB,
      "Synonyms" = symbolSynonyms,
      "GO Hallmarks" = goHallmarkTerms,
      "Max literature evidence count" = ifelse(is.na(maxLiteratureEvidencePhenotype), 
                                               NA, 
                                               paste0(maxLiteratureEvidencePhenotype," (", maxLiteratureEvidenceCount, ")")
                                               ),
      "Druggable Genome" = details.druggable_genome,
   #   "ChEMBL compounds" = tractability.smallmolecule.high_quality_compounds,
      "Chemical probes" = details.probes,
      "Target Enabling Package" = details.tep,
      "Safety Warnings" = details.safety_warnings,
      check.names = FALSE)) %>% 
    ungroup() %>%
      select(-starts_with("details."))
    #    "Druggable Genome" = details$Small_Molecule_Druggable_Genome_Member,
    #   "DrugEBIlity" = ifelse(
    #        details$DrugEBIlity_score == -1,
    #        NA,
    #        details$DrugEBIlity_score
    #    ),
    #    "PDBs with ligand" = details$PDB_Known_Ligand %>% str_remove_all("\\[|\\]|\\'"),
    #    "Safety Warning" = safety_warning,
    #    "Drugs (Drugbank)" = length(ot_details$data.drugs.drugbank[[1]]),
    #    "Drugs (ChEMBL)" = length(ot_details$data.drugs.chembl_drugs[[1]]$synonyms),
    #    "Drugs (OT evidence)" = length(ot_details$data.drugs.evidence_data[[1]]),
    #    "ChEMBL compounds" = details$High_Quality_ChEMBL_compounds,
    #    "Age of Onset Cluster" = details$age_of_onset.cluster,
    #    "GenAge ID" = details$GenAge.ID,
    #    "Description" = ot_details$data.description,

}

add_overall_mm <- function(df, diseases, minscore = 0.05){
  df %>%
    mutate(overallMultimorbidityCount =
             rowSums(select(., any_of(diseases)) >= minscore, na.rm = TRUE))
}

add_teps_probes <- function(df){
# summarise probes and links to resources
probes <- df %>% select(targetId, chemicalProbes) %>% unnest(chemicalProbes)  %>% 
  select(targetId, drugId, id, urls) %>% 
  unnest(urls) %>%
  rename(chemicalprobe = id, 
         source = niceName) %>% 
  mutate(details.probe = paste0("<a target=\"_blank\" href=\"", url, "\">", source, "</a>")) %>%
  group_by(targetId, chemicalprobe ) %>%
  summarise(details.probe_link = paste0("(",paste(details.probe, collapse = ", "),")")) %>%
  unite(details.probe_link, c(chemicalprobe, details.probe_link), sep = " ", remove = TRUE) %>% 
  summarise(details.probes = paste(details.probe_link, collapse = "; \n")) %>%
  mutate(probe = ifelse(!is.na(details.probes), 
                 "<span class=\"badge badge-dark\">Probe</span>"))

# summarise teps and links to resources
teps <-  df %>% 
  select(targetId, tep.description, tep.url) %>%
  mutate(tep = ifelse(!is.na(tep.url), "<span class=\"badge badge-dark\">TEP</span>", NA)) %>%
  mutate(details.tep = ifelse(
    !is.na(tep.url), 
    paste0("<a target=\"_blank\" href=\"", tep.url, "\"><span class=\"badge badge-dark\">", tep.description, "</span></a>"), 
    NA)) %>%
  select(targetId, tep, details.tep)

# combine
df %>% 
  left_join(probes) %>% 
  left_join(teps) %>%  
  unite("icon.tep_probe", c("tep", "probe"), sep = " ", remove = TRUE, na.rm = TRUE)
}
  




add_tractability <- function(df){
  tractability_icons <- tractability_classifications_top %>% 
    ungroup() %>% 
    select(targetId, modality, name) %>% 
    pivot_wider(names_from = modality, values_from = name) %>% 
    mutate(icon.tractability = case_when(
      stringr::str_detect(paste(PR), "Discovery_Opportunity") ~ "Discovery Opportunity",
      stringr::str_detect(paste(AB, SM), "Predicted_Tractable") ~ "Predicted Tractable",
      stringr::str_detect(paste(PR), "Literature_Precedence") ~ "Literature Precedence",
      stringr::str_detect(paste(AB, SM), "Discovery_Opportunity") ~ "Discovery Precedence",
      stringr::str_detect(paste(AB, OC, PR, SM), "Clinical_Precedence") ~ "Clinical Precedence",
      NA ~ NA_character_)
    ) %>%
    mutate(icon.tractability = ifelse(is.na(icon.tractability), NA, paste0("<span class=\"badge badge-pill badge-dark\">", icon.tractability, "</span>"))) 
  
  # druggable genome
  druggable <- df %>% 
    select(targetId, tractability) %>% 
    unnest(tractability) %>% 
    filter(id == "Druggable Family" & value == "TRUE") %>% 
    select(targetId) %>% 
    mutate(details.druggable_genome = TRUE)
  
  df %>% left_join(tractability_icons) %>% left_join(druggable)
}

add_safety <- function(df){
  safety <- df %>% select(targetId, safetyLiabilities) %>% unnest(safetyLiabilities) %>%
    group_by(targetId) %>% 
    summarise(details.safety_warnings = paste(unique(event), collapse = "; "),
              details.safety_sources = paste(unique(datasource), collapse = "; ")) %>%
    mutate(icon.safety_warning = as.character(
      shiny::icon("exclamation-triangle", lib = "font-awesome")
    )) 
  
  df %>% left_join(safety, on = targetId)
}


old_add_safety <- function(df){
  safety_liability <- df %>% 
    filter(!is.na(safety.safety_risk_info)) %>% 
    select(targetId, safety.safety_risk_info) %>% 
    unnest(safety.safety_risk_info) %>%
    unnest(references) %>%
    select(targetId, safety_liability, ref_label) %>%
    group_by(targetId, safety_liability) %>% 
    summarise(ref_label = paste(ref_label, collapse = "/")) %>%
    mutate(safety_liability = paste0(safety_liability, " [", ref_label, "].")) %>% 
    unique() %>% 
    group_by(targetId) %>% 
    summarise(safety_liability = paste(safety_liability, collapse = " "))
  
  experimental_toxicity <- df %>% 
    filter(!is.na(safety.experimental_toxicity)) %>% 
    select(targetId, safety.experimental_toxicity) %>% 
    unnest(safety.experimental_toxicity) %>% 
    reduce(data.frame) %>%
    mutate(summary = paste0("Experimental toxicity in ", assay_format_type, " ", cell_short_name, tissue, " assay [", elt, "].")) %>% 
    select(targetId = out, summary) %>% 
    unique() %>% 
    group_by(targetId) %>% 
    summarise(experimental_toxicity = paste0(summary, collapse = " "))
  
  adverse_effects <- df %>% 
    filter(!is.na(safety.adverse_effects)) %>% 
    select(targetId, safety.adverse_effects) %>% 
    unnest(safety.adverse_effects) %>% 
    unnest(organs_systems_affected) %>% 
    unnest(references) %>% 
    select(targetId, mapped_term, ref_label) %>% 
    unique() %>% 
    mutate(effect = paste0(mapped_term, " [", ref_label, "]")) %>% 
    unique() %>% 
    group_by(targetId) %>% 
    summarise(adverse_effects = paste0(effect, collapse = " and ")) %>% 
    mutate(adverse_effects = paste0("Reported adverse effects affecting the ", adverse_effects, ".")) 
  
  safety <- safety_liability %>% 
    full_join(experimental_toxicity, by = "targetId") %>% 
    full_join(adverse_effects, by = "targetId") %>% 
    mutate(icon.safety_warning = as.character(
      shiny::icon("exclamation-triangle", lib = "font-awesome")
    )) %>% 
    unite(details.safety_warnings, c(safety_liability, experimental_toxicity, adverse_effects), sep = " ", na.rm = TRUE)
  
  df %>% left_join(safety, by = "targetId")
  
}

make_legends <- function(){
  ### CellAge icons 
  cellage_legend <<- 
    data.frame(Icon = c("<i class=\"fab fa-canadian-maple-leaf\"</i><i class=\"fa fa-arrow-up\"</i>", 
                        "<i class=\"fab fa-canadian-maple-leaf\"</i><i class=\"fa fa-arrow-down\"</i>"),
               Meaning = c("Induces senescence", "Inhibits senescence"))
  
  ## Genage icons
  model_imgs <<- data.frame(
    GenAge.model.organism = c("Caenorhabditis elegans",
                              "Drosophila melanogaster",
                              "Mus musculus",
                              "Saccharomyces cerevisiae"),
    model.organism.img.path = c(
      #   "noun_worm_1817517.png",
      "noun_worm_1817517.png",
      #  "noun_Fish_657747.png",
      "noun_fly_2278070.png",
      #  "noun_Hamster_3511141.png",
      "noun_Mouse_3551360.png",
      #  "noun_mold_fungus_3292550.png",
      #  "noun_yeast_275970.png",
      "noun_yeast_275970.png"
    )
  )
  
  longevity_imgs <<- data.frame(
    GenAge.longevity.influence = c("Pro-Longevity", "Anti-Longevity", "Unclear", NA),
    GenAge.longevity.iconcode = c("arrow-up", "arrow-down", "question", NA))
  
  ## CellAge and GenAge legend
  genage_legend <<- bind_rows(
    model_imgs %>%
      filter(!is.na(model.organism.img.path)) %>%
      mutate(Icon = paste0("<img src=\"", model.organism.img.path, "\" height=\"20\"></img>")) %>%
      select(Icon, Meaning = GenAge.model.organism),
    longevity_imgs %>%
      filter(!is.na(GenAge.longevity.iconcode)) %>%
      mutate(Icon = paste0("<i class=\"fa fa-", GenAge.longevity.iconcode, "\"></i>")) %>%
      select(Icon, Meaning = GenAge.longevity.influence)
  )
  
  ## GO icons
  go_icons <<- read.csv("data/goterm_icons.csv", header = FALSE, col.names = c("goHallmark", "goIconCode", "goIcon"),quote="" )
  go_legend <- go_icons %>% select(goIcon, goHallmark)
  
  ## Save legends for app
  save(go_legend, file = "TargetAgeApp/data/hallmarks_legend.Rda")
  save(cellage_legend, file = "TargetAgeApp/data/cellage_legend.Rda")
  save(genage_legend, file = "TargetAgeApp/data/genage_legend.Rda")
  
}

# Format specific disease names 
format_specific_diseases <- function(df){
  df %>% 
    mutate(specificDiseaseName = recode(specificDiseaseName, 
                                        `low density lipoprotein cholesterol measurement` = "LDL cholesterol measurement",
                                        `high density lipoprotein cholesterol measurement` = "HDL cholesterol measurement",
                                        `very low density lipoprotein cholesterol measurement` = "VLDL cholesterol measurement",
                                        `type II diabetes mellitus` = "type 2 diabetes")) %>%
    mutate(specificDiseaseName = ifelse(specificDiseaseName != morbidity, specificDiseaseName, NA))
}


#### Data ####
load("data/analysis/target_annotations.Rda")             # Target annotations
load("data/analysis/target_tractability.Rda")
load("data/analysis/diseases_with_associations.Rda")

ard_leads %>% 
  format_specific_diseases()




# top_genes <- top_l2g %>% 
#   select(studyId, lead_variantId, targetSymbol, gene.id, L2G = yProbaModel, hasColoc, distanceToLocus) %>% 
#   mutate(targetSymbol = ifelse(hasColoc == TRUE, paste0('atop(bold("',targetSymbol,'")'), targetSymbol)) 

tbl <- g_all$cn %>%
  select(id, cluster, n_morbidities, morbidity, has_sumstats, everything(), -color, -label, -c_size) %>% 
  group_by(cluster) %>% 
  mutate(communities = length(unique(community))) %>% 
  ungroup() %>%
  filter(!is.na(cluster)) %>% 
  filter(n_morbidities > 1) %>% 
  left_join(top_l2g) %>%
  left_join(ard_leads %>% select(studyId, lead_variantId, morbidity, specificDiseaseName)) %>%
  arrange(communities, -n_morbidities, cluster) 


tbl_leads <- tbl %>%
  left_join(ard_leads %>% select(-c("morbidity", "diseaseId", "diseaseName", "specificDiseaseId",   "specificDiseaseName")) %>% unique())

tbl_genes <- tbl_leads  %>% 
  left_join(top_l2g) 


################################################

columns <- c(
  "targetId",
  "targetSymbol",
  "bioType",
  "symbolSynonyms",
  "chemicalProbes"
)

## add cellage, genage, and hallmarks of ageing 
make_legends()
tbl_targets <- target_annotations %>% 
  rowwise() %>% 
  mutate(symbolSynonyms = paste(symbolSynonyms, collapse='; '),
         biotype = str_replace_all(biotype, "_", " ")) %>% 
  ungroup() %>%
  add_genage_icons() %>%
  add_cellage_icons() %>%
  add_go_hallmarks() %>% 
  add_tractability() %>% select(-tractability) %>% 
  add_teps_probes() %>% select(-chemicalProbes) %>% 
  add_safety() %>% 
  add_overall_mm(diseases = d, minscore = 0.05) %>% 
  add_details_column()

save(tbl_targets, file = "TargetAgeApp/data/prepared_table.Rda")
load(file = "data/analysis/targetage_geneids.Rda")
save(targetage, file = "TargetAgeApp/data/targetage_geneids.Rda")
save(d, file = "TargetAgeApp/data/diseases_with_associations.Rda")



## Genetics data 
load("data/analysis/ltg_all.Rda")
load("data/analysis/graph_all_morbidities.Rda")
load("data/analysis/ard_leads_filtered.Rda")

load("data/analysis/top_l2g.Rda")

ard_leads <- ard_leads %>% 
  mutate(specificDiseaseName = recode(specificDiseaseName, 
                                      `low density lipoprotein cholesterol measurement` = "LDL cholesterol measurement",
                                      `high density lipoprotein cholesterol measurement` = "HDL cholesterol measurement",
                                      `very low density lipoprotein cholesterol measurement` = "VLDL cholesterol measurement",
                                      `type II diabetes mellitus` = "type 2 diabetes")) %>%
  mutate(specificDiseaseName = ifelse(specificDiseaseName != morbidity, specificDiseaseName, NA))

top_l2g <- l2g_all_joined %>%
  group_by(studyId, lead_variantId) %>% 
  top_n(1,yProbaModel ) %>% ## need to remove duplicates for 5 variants 
  select(studyId, lead_variantId, gene.symbol, gene.id, L2G = yProbaModel, hasColoc, distanceToLocus)

top_genes <- top_l2g %>%
  mutate(gene.symbol = ifelse(hasColoc == TRUE, paste0('atop(bold("',gene.symbol,'")'), gene.symbol)) %>%
  top_n(1, L2G) %>%
  summarise(gene.symbol = paste0(gene.symbol, collapse = ", "),
            L2G = max(L2G))

cluster_tbl <- g_all$cn %>%
  select(id, cluster, n_morbidities, morbidity, has_sumstats, everything(), -color, -label, -c_size) %>% 
  group_by(cluster) %>% 
  mutate(communities = length(unique(community))) %>% 
  ungroup() %>%
  filter(!is.na(cluster)) %>% 
  filter(n_morbidities > 1) %>% 
  left_join(top_l2g) %>%
  left_join(ard_leads %>% select(studyId, lead_variantId, morbidity, specificDiseaseName)) %>%
  group_by(cluster) %>% 
  mutate(targetIds = paste(unique(gene.id), collapse = ";"),
         targetSymbols = paste(unique(gene.symbol), collapse = ";")) %>% 
  arrange(communities, -n_morbidities, cluster) 

tbl_leads <- cluster_tbl %>%
  left_join(ard_leads %>% 
              select(-c("morbidity", "diseaseId", "diseaseName", "specificDiseaseId",   "specificDiseaseName")) %>% 
              unique())

tbl_genes <- tbl_leads  %>% 
  left_join(top_l2g) 

save(tbl_genes, file = "TargetAgeApp/data/genetics_table.Rda")
save(g_all, file = "TargetAgeApp/data/graph_all_morbidities.Rda")


