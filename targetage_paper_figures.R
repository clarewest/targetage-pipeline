library(tidyverse)

## Global settings
default_save_dir = "data/analysis/"  # save location for intermediate Rda files
paper_save_dir = "data/paper_data/"  # save location for figures/tables for paper
default_overwrite = FALSE            # if FALSE will load previously saved Rda files

### Age-related diseases/traits
diseases <- read.csv("data/full_disease_list.csv")

### Lead variants for genetic associations from all GWASs for age-related diseases/traits
ard_leads <- get_lead_variants(overwrite = default_overwrite)

### Diseases with GWAS evidence (some in `diseases` don't have any)
d <- ard_leads$morbidity %>% unique()

############## SI Fig 1

### Overlaps between genetic associations
overlaps <- get_overlaps(overwrite = default_overwrite)

## Varying overlap threshold
jaccard_communities_results <- get_jaccard_cutoff_communities(overwrite = default_overwrite)
ji_figure <- plot_jaccard_figure(overwrite = default_overwrite)
get_n_communities_per_cluster()


############# Fig 2

## Required data 
target_annotations <- get_associations_annotations(overwrite = default_overwrite, diseases = d)
load(file = paste0(default_save_dir, "targetage_geneids.Rda"))
load(file = paste0(default_save_dir, "targets_all_go.Rda"))

## Prepare for enrichment analysis
## (Do in separate session to avoid painful interactions)
prepare_enrichment_files(overwrite = TRUE)

## Make figures
hallmarks_bar <- plot_hallmarks_barplot()
venn <-plot_venn_diagram()
fig <- plot_go_figure(overwrite = default_overwrite)


###########################################################################
###########################################################################
###                                                                     ###
###                       SECTION 1:                                    ###
###                       FUNCTIONS                                     ###
###                                                                     ###
###########################################################################
###########################################################################

plot_hallmarks_barplot <- function(){
  all_gen_assoc_targets <- target_annotations %>% select(targetId, all_of(d)) %>% filter_if(is.numeric, any_vars(. > 0.5)) %>% pull(targetId)
  
  go_compare <- target_go %>% 
    filter(targetId %in% all_gen_assoc_targets) %>% 
    mutate(group = "all") %>% 
    bind_rows(target_go %>% 
                filter(targetId %in% targetage) %>% 
                mutate(group = "targetage")) %>%
    group_by(group, targetId, goHallmark) %>% 
    count() %>% 
    pivot_wider(names_from = goHallmark, values_from = n, values_fill = 0) %>% 
    rename(other = `NA`) %>% 
    pivot_longer(-c("group", "targetId"))
  
  ## number of targets with at least one annotation in each hallmark
  go_compare_hallmarks <- go_compare %>% 
    ungroup() %>% 
    group_by(group) %>% 
    mutate(total = length(unique(targetId))) %>% 
    ungroup() %>% 
    group_by(group, total, name) %>% 
    summarise(n = sum(value>0)) %>% 
    mutate(proportion = n/total)
  
  ## bar plot 
  ## but need to make sure I'm only using >0.5 genes 
  hallmarks_bar <- ggplot(subset(go_compare_hallmarks, name!="other"), aes(x=proportion, y = reorder(name, proportion), fill = group, label = paste0(" ", scales::percent(proportion, accuracy = 0.1)))) + 
    geom_col(position = "dodge") + 
    theme_bw() + 
    geom_text(position = position_dodge(width = 0.9), hjust = 0, size = 3) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 10, colour = "black"), axis.title.y = element_blank(), legend.position = c(0.8, 0.1), legend.title = element_blank()) + 
    scale_x_continuous(limits = c(0, 0.15), labels = scales::percent_format(accuracy = 1), expand = c(0, NA)) + labs(x = "Percentage of targets", subtitle = "Proportion of targets mapped to\nHallmarks of Ageing GO terms")
  
  return(hallmarks_bar)
}

## Prepare files for enrichment analysis
prepare_enrichment_files <- function(overwrite = default_overwrite, save_dir = paper_save_dir){
  if (overwrite){
    genage_ids <- target_annotations %>% 
      select(targetId, GenAge) %>%
      unnest(cols = c(GenAge)) %>% 
      filter(!is.na(entrez.gene.id)) %>%
      pull(entrez.gene.id) %>%
      unique()
    
    cellage_ids <- target_annotations %>% 
      select(targetId, CellAge) %>%
      unnest(cols = c(CellAge)) %>% 
      filter(!is.na(CellAge.entrezid)) %>%
      pull(CellAge.entrezid) %>%
      unique()
    
    gene_sets_entrez <- list(TargetAge = targetage,
                             GenAgeHuman = genage_ids,
                             CellAge = cellage_ids)
    save(gene_sets_entrez, file = paste0(paper_save_dir,"gene_sets_entrez.Rda"))
  }
}

test_overlap <- function(setA, setB){
  ji = sum(setA %in% setB) / length(unique(c(setA, setB)))
  return(ji)
}

## Make a venn diagram 
plot_venn_diagram <- function(){
  library(ggvenn)
  genage_ids <- target_annotations %>% 
    select(targetId, GenAge) %>%
    unnest(cols = c(GenAge)) %>% 
    filter(!is.na(entrez.gene.id)) %>%
    pull(targetId) %>%
    unique()
  
  cellage_ids <- target_annotations %>% 
    select(targetId, CellAge) %>%
    unnest(cols = c(CellAge)) %>% 
    filter(!is.na(CellAge.entrezid)) %>%
    pull(targetId) %>%
    unique()
  
  gene_sets <- list(TargetAge = targetage,
                    GenAgeHuman = genage_ids,
                    CellAge = cellage_ids)
  venn <- ggvenn(gene_sets, fill_color = c("white", "white", "white"), show_percentage = FALSE, text_size = 4, set_name_size = 4, stroke_size = 0.5)
  return (venn)
}

plot_go_figure <- function(overwrite = default_overwrite, save_dir = paper_save_dir){
  library(patchwork)
  load(file = "scripts/enrichment_plot_1.Rda")
  load(file = "scripts/enrichment_plot_2.Rda")
  load(file = "scripts/enrichment_plot_3.Rda")
  layout <- "
14
14
25
25
35
35
"
  fig <- p_go_BP + p_go_MF + p_reactome + venn + hallmarks_bar  +
    plot_layout(design = layout,
                widths = c(1,1),
                #    heights = c(1,1,1,1,2),
                guides = "collect")  +
    plot_annotation(tag_levels = list(c("a", "", "", "b", "c"))) & 
    theme(plot.subtitle = element_text(size = 12), plot.tag.position = c(0, 1), legend.box="vertical", plot.tag = element_text(size = 14))
  
  if (overwrite){
    ggsave(fig, file = paste0(paper_save_dir, "enrichment_figure.png"), width = 8, height = 9, dpi = 300)
    ggsave(fig, file = paste0(paper_save_dir, "enrichment_figure.pdf"), width = 14, height = 13)
  }
  
  return(fig)
  
}

get_jaccard_cutoff_communities <- function(overwrite = default_overwrite, save_dir = paper_save_dir) {
  
  if (overwrite){
    
    jaccard_communities_results <- data.frame()
    
    ## Iterate over jaccard index cutoffs, get the number of clusters and sub-cluster communities
    ## Warning: computationally expensive
    for (curr_min_jaccard in seq(0, 1, 0.1)){
      print(curr_min_jaccard)
      g_jaccard <- count_communities(d, ard_leads, coloc_within, overlap_within, min_jaccard = curr_min_jaccard, detect_subgraph_communities = TRUE, plot = FALSE)
      jaccard_communities_results <- bind_rows(g_jaccard$cn %>% mutate(min_jaccard = curr_min_jaccard), jaccard_communities_results)
    }
    save(jaccard_communities_results, file = paste0(save_dir, "jaccard_communities.Rda"))
  }
  else {
    load(paste0(save_dir, "jaccard_communities.Rda"))
  }
  return(jaccard_communities_results)
}

plot_jaccard_figure <- function(overwrite = default_overwrite, save_dir = paper_save_dir){
  ji_n_clusters <- jaccard_communities_results %>% 
    mutate(type = "all") %>% 
    bind_rows(jaccard_communities_results %>% 
                filter(n_morbidities > 1) %>% 
                mutate(type = ">1 morbidity")) %>% 
    group_by(type, min_jaccard) %>% 
    summarise(n_clusters = length(unique(cluster))) %>% 
    ggplot(., aes(x=min_jaccard, y = n_clusters)) + 
    geom_line() + 
    geom_point() + 
    facet_wrap(~type, 
               labeller = as_labeller(c(`>1 morbidity` = "Clusters with >1 morbidity", `all` = "All clusters"))) + 
    theme_bw() + 
    labs(x = "Jaccard index threshold", 
         y = "Number of disconnected components (clusters)", 
         subtitle = "Number of clusters detected with varying overlap stringency") + 
    theme(panel.grid.minor.x = element_blank(), 
          panel.grid.major.x = element_blank()) + 
    scale_x_continuous(breaks = seq(0, 1, 0.2)) + 
    scale_y_continuous(limits = c(0, 3500), expand = c(0, NA)) 
  
  ji_n_communities <-  jaccard_communities_results %>% 
    mutate(type = "all") %>% 
    bind_rows(jaccard_communities_results %>% 
                filter(n_morbidities > 1) %>% 
                mutate(type = ">1 morbidity")) %>%  
    group_by(type, min_jaccard, cluster) %>% 
    summarise(communities = length(unique(community))) %>% 
    ungroup() %>% 
    group_by(type, min_jaccard) %>% 
    mutate(outlier.high = communities > quantile(communities, .75) + 1.50*IQR(communities)) %>% 
    mutate(outlier.communities = ifelse(outlier.high == TRUE, communities, NA)) %>% 
    ungroup() %>% 
    ggplot(., aes(y= as.character(min_jaccard), x = communities)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(aes(x=outlier.communities), size  = .2, width = 0.4, height = 0.3, alpha = 0.7) +
    theme_bw()  + 
    labs(x = "Number of communities per disconnected component (cluster)", y = "Jaccard index threshold", subtitle = "Number of communities detected in each cluster with varying overlap stringency") +
    facet_wrap(~type, labeller = as_labeller(c(`>1 morbidity` = "Clusters with >1 morbidity", `all` = "All clusters"))) + 
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + scale_x_continuous(breaks = seq(2,10,2))
  
  ## Assemble into figure
  ji_si_figure <- ji_n_clusters / ji_n_communities + plot_annotation(tag_levels = 'a')
  
  if (overwrite){
    ggsave(ji_si_figure, file = paste0(save_dir, "si_jaccard_index.png"), width = 8, height = 7.5, dpi = 300)
    ggsave(ji_si_figure, file = paste0(save_dir,"si_jaccard_index.pdf"), width = 8, height = 7.5)
  }
  return(ji_si_figure)
}

get_n_communities_per_cluster <- function(){
  chosen <- jaccard_communities_results %>% filter( min_jaccard == 0)
  N_clusters_with_communities <- chosen %>%
    mutate(type = "all") %>% 
    bind_rows(chosen %>% filter(n_morbidities > 1) %>% mutate(type = ">1 morbidity")) %>% 
    group_by(type, cluster) %>% 
    summarise(communities = length(unique(community))) %>% 
    ungroup() %>% 
    group_by(type,communities) %>% 
    count() %>% 
    ungroup() %>% 
    group_by(type) %>% 
    mutate(total = sum(n), proportion = round((n/total)*100,2))
  return(N_clusters_with_communities)
}


