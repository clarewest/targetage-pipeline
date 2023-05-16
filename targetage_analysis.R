library(tidyverse)
source("targetage_analysis_functions.R")

## Global settings
default_save_dir = "data/analysis/"  # save location for intermediate Rda files
default_overwrite = FALSE            # if FALSE will load previously saved Rda files

###########################################################################
###########################################################################
###                                                                     ###
###                       SECTION 1:                                    ###
###                     Read in data                                    ###
###                                                                     ###
###########################################################################
###########################################################################

### Age-related diseases/traits
diseases <- read.csv("data/full_disease_list.csv")

### Lead variants for genetic associations from all GWASs for age-related diseases/traits
ard_leads <- get_lead_variants(overwrite = default_overwrite, input_file = "data/targetage/ard_leads.parquet/part-00000-682df92e-1d00-4977-9cd7-1834de53b5c2-c000.snappy.parquet")

### GWAS studies
gwas_ids <- ard_leads %>% select(c('studyId', 'trait_reported')) %>% unique()
#write.csv(gwas_ids, paste0(default_save_dir,"SITable5.csv"))

### Diseases with GWAS evidence (some in `diseases` don't have any)
d <- ard_leads$morbidity %>% unique()
#save(d, file = "data/analysis/diseases_with_associations.Rda")

### Overlaps between genetic associations
overlaps <- get_overlaps(overwrite = default_overwrite, 
                         coloc_input_file = "data/targetage/coloc_ard_leads.parquet/part-00000-9479cf2c-da10-447a-ab81-7c1750b5bc78-c000.snappy.parquet",
                         overlap_input_file = "data/targetage/overlap_ard_leads.parquet/part-00000-5b27ad8f-95fb-4b3b-85a9-bc24c5e30ea1-c000.snappy.parquet")

### Target details and genetic associations with each phenotype
target_annotations <- get_associations_annotations(overwrite = default_overwrite, 
                                                   annotations_input_file = "data/targetage/ard_annotations.parquet/part-00000-16237769-f92d-48c3-b396-6d5f5d797719-c000.snappy.parquet",
                                                   associations_input_file = "data/targetage/ard_associations.parquet/part-00000-b1ee0a40-e1b1-4ee6-9b41-fc41115d0a05-c000.snappy.parquet",
                                                   go_input_file = "data/full_goterm_list.csv",
                                                   diseases = diseases)


###########################################################################
###########################################################################
###                                                                     ###
###                       SECTION 2:                                    ###
###                     Graph analysis                                  ###
###                                                                     ###
###########################################################################
###########################################################################

### Graph analysis
g_all <- create_targetage_graph(overwrite = default_overwrite, save_dir = default_save_dir)

## Get the number of clusters within a disease (i.e. number of individual signals)
individual_diseases = get_diseases_individually(ard_leads, overlaps, overwrite = TRUE, save_dir = default_save_dir)

## Get the number of clusters with just one community
n_multicommunities <- get_n_multicommunity_clusters(g_all)

### Calculate overlap between each combination of diseases 
ft_all <- get_pairwise_overlaps(g_all)
qq_plot <- get_qq_plot(ft_all, overwrite = TRUE)

## Make a heatmap
shared_heatmap <- plot_shared_heatmap(ft_all, detailed = FALSE, overwrite = FALSE)
shared_heatmap <- plot_shared_heatmap(ft_all, detailed = TRUE, overwrite = FALSE)

### Table of clusters with number of nodes and morbidities
tbl <- g_all$cn %>% 
  group_by(cluster, n_morbidities) %>% 
  summarise(nodes = length(id), morbidities = paste0(unique(morbidity), collapse = ", ")) %>% 
  arrange(-n_morbidities) %>% filter(!is.na(cluster))

###########################################################################
###########################################################################
###                                                                     ###
###                       SECTION 3:                                    ###
###           Get genes for each cluster via OT L2G                     ###
###                                                                     ###
###########################################################################
###########################################################################

### Get L2G from Open Targets Genetics
l2g_all_joined <- get_L2G(overwrite = default_overwrite)

### Get details of xQTL colocalisation 
l2g_qtl_coloc <- get_L2G_coloc(overwrite = default_overwrite)

## Top gene per lead variant 
top_l2g <- l2g_all_joined %>%
  group_by(studyId, lead_variantId) %>% 
  top_n(1,yProbaModel ) %>% 
  ## 4 variants have more than one gene with identical L2G 
  ## in this case choose the closest by distanceToLocus
  ## 3 have L2G < 0.5 and therefore aren't in TargetAge anyway
  ## 1 has L2G > 0.5 and maps to two genes
  ## putative ENSG00000285901/AC008012.1 and CCND2 instead 
  ## choose one (CCND2)
  top_n(-1, distanceToLocus) %>%
  slice(1)
#save(top_l2g, file = paste0(default_save_dir, "top_l2g.Rda"))

## Add top gene to each node
# NB this df contains a row for every node but L2G was only retrieved for nodes:
# - with >1 morbidity
# - in a cluster (i.e. !is.na(cluster))
# Remaining nodes without an L2G (~500) didn't have any genes with >0.05 L2G
top_cluster_genes <- g_all$cn %>% 
  left_join(top_l2g) %>% 
  left_join(tbl) 

## Top genes for nodes in clusters with more than one morbidity
top_mm_cluster_genes <- top_cluster_genes %>% 
  filter(n_morbidities>1)

## How many clusters is each gene associated with 
top_cluster_genes_summarised <- top_mm_cluster_genes %>% 
  filter(yProbaModel >= 0.5) %>% 
  ungroup() %>% 
  group_by(gene.id, targetSymbol, cluster) %>% 
  count() %>% 
  filter(!is.na(cluster)) %>% 
  ungroup() %>% 
  group_by(gene.id, targetSymbol) %>% 
  summarise(n_clusters = length(cluster)) %>% 
  arrange(-n_clusters)

top_cluster_genes_summarised %>% 
  ungroup() %>% 
  count(n_clusters) %>% 
  mutate(prop = n/sum(n)*100)


###########################################################################
###########################################################################
###                                                                     ###
###                       SECTION 4:                                    ###
###                   GWAS summary information                          ###
###                                                                     ###
###########################################################################
###########################################################################


### Look at ancestry of GWAS
populations <- get_ancestries(ard_leads)
populations %>% ungroup %>% count(single_population, ancestry) %>% mutate(p = n/sum(n)*100) %>% arrange(-p)

### Update SI table in Google Drive with GWAS and cluster details
# The list of ARDs, number of GWAS studies, and number of (lead) variants, proportion with summary statistics, number of communities
genetics_tables <- get_genetics_table(ard_leads, diseases, g_all, individual_diseases)
#update_gs_genetics_table(genetics_tables)


###########################################################################
###########################################################################
###                                                                     ###
###                       SECTION 5:                                    ###
###          TargetAge genes annotations and tractability               ###
###                                                                     ###
###########################################################################
###########################################################################


# all multimorbidity genes from clusters
targetage_all <- top_cluster_genes_summarised$gene.id 

## protein coding and has OT annotations
targetage <- targetage_all[targetage_all %in% target_annotations$targetId]
targetage_cluster_summarised = subset(top_cluster_genes_summarised, gene.id %in% targetage)

targetage_cluster_summarised %>% 
  ungroup() %>% 
  count(n_clusters) %>% 
  mutate(prop = n/sum(n)*100)

targetage_annotations <- target_annotations %>% filter(targetId %in% targetage)

#write.table(targetage_cluster_summarised, file = paste0(default_save_dir, "top_L2G_genes.txt"), sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
#save(targetage, file = paste0(default_save_dir, "targetage_geneids.Rda"))
load(file = paste0(default_save_dir, "targetage_geneids.Rda"))


###################################################
## Overlap with GenAge and CellAge and Hallmarks ##
###################################################

## Get gene sets
gene_sets = get_gene_sets(overwrite = default_overwrite)
gene_sets_entrez = get_entrez_gene_sets(overwrite = default_overwrite)

## Venn Diagram of overlaps
gplots::venn(gene_sets[1:3])

## Test overlaps
cellage_overlap = test_overlap(gene_sets, "TargetAge", "CellAge")
genage_overlap = test_overlap(gene_sets, "TargetAge", "GenAge")
hallmarks_overlap = test_overlap(gene_sets, "TargetAge", "Hallmarks")
hallmarks_a_overlap = test_overlap(gene_sets, "ARDs", "Hallmarks")
hallmarks_g_overlap = test_overlap(gene_sets, "GenAge", "Hallmarks")
hallmarks_c_overlap = test_overlap(gene_sets, "CellAge", "Hallmarks")



###################################################
##                  Enrichment                   ##
###################################################

all_enriched <- perform_enrichment(overwrite = default_overwrite)

top_enriched <- all_enriched %>% 
  filter(set != "Reactome") %>% 
  group_by(set, Cluster) %>% slice_min(p.adjust, n = 10)
plot_enriched <- all_enriched %>% filter(ID %in% top_enriched$ID)

ggplot(plot_enriched, aes(x=Cluster, y = reorder(Description, p.adjust))) + geom_point() + facet_wrap(~set, ncol = 1, strip.position = "right", scales = "free_y") + theme_bw()

hallmark_terms <- read.csv("data/full_goterm_list.csv")
clusterprofiler_enriched_hallmarks <- all_enriched %>% filter(ID %in% hallmark_terms$goId) %>% filter(p.adjust<0.05)

###################################################
##                  Enrichment                   ##
###################################################

## Hallmarks
targetage_hallmarks <- targetage_annotations %>% select(targetId, ageingHallmarks) %>% unnest(ageingHallmarks)

hallmark_counts <- get_hallmark_counts(targetage_hallmarks, overwrite = default_overwrite)

enriched_hallmarks <- enrich_hallmarks(gene_sets) 

hallmarks_barplot <- plot_hallmarks_barplot(enriched_hallmarks, overwrite = default_overwrite)
hallmarks_barplot_ta <- plot_mini_hallmarks_barplot(enriched_hallmarks, overwrite = default_overwrite)

## OT Tractability

#targetage_tract <- get_subset_annotations(targetage_annotations, targetage, "tractability")
targetage_tract <- targetage_annotations %>% 
  filter(targetId %in% targetage) %>% 
  dplyr::select(targetId, targetSymbol, all_of("tractability"))  %>% 
  unnest("tractability")

# N with approved drug:
tractability_wide <- targetage_tract %>% group_by(targetSymbol, modality, id) %>% pivot_wider(names_from = "id", values_from = "value")
## get categories
tractability_classifications <- 
tractability_wide %>% 
  mutate(
    `Clinical_Precedence` = any(`Approved Drug`, `Advanced Clinical`, `Phase 1 Clinical`) == TRUE,
    `Discovery_Opportunity` = any(`Structure with Ligand`, `High-Quality Ligand`, `UniProt Ubiquitination`, `Database Ubiquitination`, `Half-life Data`, `Small Molecule Binder`) == TRUE,
    `Literature_Precedence` = any(`Literature`) == TRUE,
    `Predicted_Tractable` = any(`High-Quality Pocket`, `Med-Quality Pocket`, `Druggable Family`) == TRUE,
    `Predicted_Tractable_(H)` = any(`GO CC high conf`, `UniProt loc high conf`) == TRUE,
    `Predicted_Tractable_(M/L)` = any(`UniProt loc med conf`, `UniProt SigP or TMHMM`, `GO CC med conf`) == TRUE,
    `No_Evidence` = TRUE
    ) %>% 
  select(targetId, targetSymbol, modality, contains("_"))

classes = c("Clinical_Precedence", "Literature_Precedence", "Discovery_Opportunity", "Predicted_Tractable", "Predicted_Tractable_(H)", "Predicted_Tractable_(M/L)", "No_Evidence" )

class_order = data.frame(
  name = factor(classes, levels = classes),
  priority = c(1:length(classes))
)

tractability_classifications_long <- tractability_classifications %>% pivot_longer(-c(targetId, targetSymbol, modality)) %>% right_join(class_order)
tractability_classifications_top <- tractability_classifications_long %>% group_by(targetId, targetSymbol, modality) %>% filter(value==TRUE) %>% slice(which.min(priority))
tractability_classifications_summary <- tractability_classifications_top %>% 
  ungroup() %>% group_by(modality) %>% 
  count(name) %>% 
  rename(top_category = "name") %>%
  mutate(modality = factor(modality, levels = c("SM", "AB", "PR", "OM")),
         top_category = factor(gsub("_", " ", top_category), levels = gsub("_", " ", classes)))
save(tractability_classifications_top, file = paste0(default_save_dir, "target_tractability.Rda"))


clinical_targets <- tractability_classifications_top %>% filter(name=="Clinical_Precedence") %>% group_by(targetSymbol) %>% count()
discovery_precedence <- tractability_classifications_top %>% filter(! targetSymbol %in% clinical_targets$targetSymbol) %>% filter(name == "Discovery_Opportunity") %>% filter(modality == "SM")
predicted_tractable <- tractability_classifications_top %>% filter(! targetSymbol %in% clinical_targets$targetSymbol) %>% filter(!targetSymbol %in% discovery_precedence$targetSymbol) %>% filter(name %in% c("Predicted_Tractable", "Predicted_Tractable_(H)", "Predicted_Tractable_(M/L)")) %>% filter(modality != "PR") %>% group_by(targetSymbol) %>% count()

tractability_plot <- 
  tractability_classifications_summary %>% 
  filter(modality != "OC") %>% 
  ggplot(., aes(x=modality, y = n, fill = top_category, label = ifelse(n>10, paste0(top_category, " (", n, ")"), NA))) + 
  geom_col() + 
  geom_text(colour = "white", position = position_stack(vjust = 0.5), size = 3) + 
  theme_classic() + 
  theme(legend.position = "bottom", panel.grid=element_blank(), axis.title.x = element_blank()) + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_discrete(limits = c("SM", "AB", "PR"), labels = c("Small Molecule", "Antibody", "PROTAC")) + 
  scale_fill_manual(values = c(viridis::viridis(7)[1:6], "grey")) + 
  guides(fill = "none") + 
  labs(y = "Number of targets", subtitle = "Open Targets Target Tractability Assessment")


blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

simple_tractability_plot <- 
  tractability_classifications_summary %>% 
  filter(modality %in% c("AB", "SM")) %>% 
    mutate(top_category = as.character(top_category),
           modality = factor(modality, levels = c("AB", "SM"))) %>%
  mutate(top_category = ifelse(stringr::str_detect(top_category, "Predicted Tractable"), "Predicted Tractable", top_category)) %>%
    mutate(top_category = factor(top_category, levels = c("Clinical Precedence", "Discovery Opportunity", "Predicted Tractable", "No Evidence"))) %>% 
    group_by(modality, top_category) %>% summarise(n= sum(n)) %>%
    ggplot(., aes(x=modality, y = n, fill = top_category, label = n)) + 
    geom_col() + 
    geom_text(colour = "white", position = position_stack(vjust = 0.5), size = 3) + 
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.border = element_blank(),
      panel.grid=element_blank(),
      axis.ticks = element_blank(),
      plot.title=element_text(size=14, face="bold"),
      legend.position = "bottom",
      legend.direction = "vertical", 
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      legend.title = element_blank(),
      legend.margin=margin(0,0,0,0),
      legend.box.margin=margin(-10,-10,-10,-10)
      
    ) + 
   # theme(legend.position = "bottom", panel.grid=element_blank(), axis.title.x = element_blank()) + 
    scale_y_continuous(expand = c(0,0)) + 
    scale_x_discrete(limits = c("AB","SM"), labels = c("Antibody", "Small Molecule")) + 
    scale_fill_manual(values = c(viridis::viridis(4)[1:3], "grey")) + 
    coord_polar("y", start = 0)
ggsave(simple_tractability_plot, width = 3, height = 2.5, dpi = 600, file = paste0(default_save_dir, "tractability_mini.png"))
ggsave(simple_tractability_plot, width = 3, height = 2.5, file = paste0(default_save_dir, "tractability_mini.pdf"))


library(patchwork)
tract_fig <- hallmarks_barplot + tractability_plot + plot_annotation(tag_levels = 'a') + plot_layout(widths = c(0.65, 2.35))
ggsave(tract_fig, width = 12.1, height = 4.5, dpi = 600, file = paste0(default_save_dir, "fig3.png"))
ggsave(tract_fig, width = 12.1, height = 4.5, file = paste0(default_save_dir, "fig3.pdf"))


## Chemical Probes 

#probes <- get_subset_annotations(target_annotations, targetage, "chemicalProbes")
probes <- targetage_annotations %>% 
  filter(targetId %in% targetage) %>% 
  dplyr::select(targetId, targetSymbol, all_of("chemicalProbes"))  %>% 
  unnest("chemicalProbes")

probes %>% group_by(targetId, targetSymbol) %>% count()

## 38 targets have high quality probes 
high_quality_probes <- probes %>% filter(isHighQuality==TRUE) %>% group_by(targetId, targetSymbol) %>% count()
top_any_modality <- tractability_classifications_long %>% ungroup() %>% group_by(targetId, targetSymbol) %>% filter(value==TRUE) %>% slice(which.min(priority)) %>% ungroup()

## 11 have no clinical precedence
high_quality_probes_tractability <- high_quality_probes %>% left_join(top_any_modality) 
high_quality_probes_tractability %>% group_by(name) %>% count()

## SI Probe Table 
probe_si_table <- probes %>% 
  filter(isHighQuality==TRUE) %>% 
  left_join(top_any_modality) %>% 
  select(targetId, targetSymbol, inchiKey, drugId, id, control, urls, tractability = name) %>% 
  unnest(urls) %>% 
  arrange(tractability, targetSymbol)

write.csv(probe_si_table, paste0(default_save_dir,"SITable4.csv"))


