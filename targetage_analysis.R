library(tidyverse)

## Global settings
default_save_dir = "data/analysis/"  # save location for intermediate Rda files
default_overwrite = FALSE            # if FALSE will load previously saved Rda files

###########################################################################
###########################################################################
###                                                                     ###
###                       SECTION 2:                                    ###
###                        ANALYSIS                                     ###
###                                                                     ###
###########################################################################
###########################################################################

### Age-related diseases/traits
diseases <- read.csv("data/full_disease_list.csv")

### Lead variants for genetic associations from all GWASs for age-related diseases/traits
ard_leads <- get_lead_variants(overwrite = default_overwrite)

### Diseases with GWAS evidence (some in `diseases` don't have any)
d <- ard_leads$morbidity %>% unique()
#save(d, file = "TargetAge/data/diseases_with_associations.Rda")

### Overlaps between genetic associations
overlaps <- get_overlaps(overwrite = default_overwrite)

### Target details and genetic associations with each phenotype
target_annotations <- get_associations_annotations(overwrite = default_overwrite, diseases = d)

### Graph analysis
g_all <- create_targetage_graph(overwrite = default_overwrite, save_dir = default_save_dir)

### Calculate overlap between each combination of diseases 
ft_all <- get_pairwise_overlaps(g_all)
qq_plot <- get_qq_plot(ft_all)

## need to update this figure
get_qq_plot <- function(pairwise_overlaps){
  x <- 10^(-seq(0,8,0.1))
  y <- runif(100000)
  df <- data.frame(x = -log10(x), y = -quantile(log10(pairwise_overlaps$p.value), probs = x), control = -quantile(log10(y), probs = x))
  gg <- ggplot(data = df) + geom_point(aes(x = x, y = y)) + geom_abline(slope = 1)
  gg
}

## Make a heatmap
shared_heatmap <- plot_shared_heatmap(ft_all, overwrite = FALSE)

### Get L2G from Open Targets Genetics
l2g_all_joined <- get_L2G(overwrite = default_overwrite)

### Get details of xQTL colocalisation 
l2g_qtl_coloc <- get_L2G_coloc(overwrite = default_overwrite)

### Table of clusters with number of nodes and morbidities
tbl <- g_all$cn %>% 
  group_by(cluster, n_morbidities) %>% 
  summarise(nodes = length(id), morbidities = paste0(unique(morbidity), collapse = ", ")) %>% 
  arrange(-n_morbidities) %>% filter(!is.na(cluster))

## Top gene per lead variant 
top_l2g <- l2g_all_joined %>%
  group_by(studyId, lead_variantId) %>% 
  top_n(1,yProbaModel ) %>% 
  ## 3 variants have more than one gene with identical L2G 
  ## in this case choose the closest by distanceToLocus
  ## all 3 have L2G < 0.5 and therefore aren't in TargetAge anyway
  top_n(-1, distanceToLocus)
save(top_l2g, file = paste0(save_dir, "top_l2g.Rda"))


## Add top gene to each node
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

targetage <- top_cluster_genes_summarised$gene.id
save(targetage, file = paste0(default_save_dir, "targetage_geneids.Rda"))

write.table(top_cluster_genes_summarised, file = paste0(save_dir, "top_L2G_genes.txt"), sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)

#top_l2g <- read.table(file = "top_L2G_genes.txt")

## Look at tractability

targetage_tract <- target_annotations %>% filter(targetId %in% targetage) %>% jsonlite::flatten()
tract_plot_data <- targetage_tract %>% 
  mutate(tractability.other_modalities.top_category = ifelse(tractability.other_modalities.categories.clinical_precedence > 0, "Clinical_Precedence", NA)) %>% 
  select(targetId, ends_with("top_category")) %>% 
  pivot_longer(-targetId, names_to = "modality", values_to = "top_category") %>% 
  count(modality, top_category) %>% 
  mutate(modality = gsub(".*[.]([^.]+)[.].*", "\\1", modality), top_category = gsub("_", " ", top_category)) %>% 
  mutate(top_category = gsub("le ab", "le \n", top_category), top_category = gsub(" sm", "", top_category), top_category = gsub(" ab", "", top_category)) 

tractability_plot <- ggplot(tract_plot_data, aes(x=modality, y = n, fill = top_category, label = ifelse(n>10, paste0(top_category, " (", n, ")"), NA))) + 
  geom_col() + 
  geom_text(colour = "white", position = position_stack(vjust = 0.5), size = 3) + 
  theme_bw() + 
  theme(legend.position = "bottom", panel.grid=element_blank(), axis.title.x = element_blank()) + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_discrete(limits = c("smallmolecule", "antibody", "other_modalities"), labels = c("Small Molecule", "Antibody", "Other Modalities")) + 
  guides(fill = FALSE) + 
  labs(y = "Number of targets", subtitle = "Open Targets Target Tractability Assessment")

tdb <- readxl::read_excel("data/targetDB/Export_644_entries_15Jun2021_105820.xlsx", skip = 1)
tdb_plot <- tdb %>% 
  group_by(Tractable, Pharos_class) %>% 
  count() %>%
  mutate(Tractable = factor(Tractable, levels = c("Tractable", "Challenging", "Intractable")),
         Pharos_class = factor(Pharos_class, levels = c("Tdark", "Tbio", "Tchem", "Tclin"))) %>% 
  ggplot(., aes(x = Tractable, fill = Pharos_class, y = n, label = ifelse(n>3, paste0(Pharos_class, " (", n, ")"), NA))) + 
  geom_col() + 
  theme_bw() + 
  scale_fill_discrete("Pharos class") + 
  scale_y_continuous(expand = expansion(mult = c(0, .1))) + 
  geom_text(colour = "white", position = position_stack(vjust = 0.5), size = 3) + 
  guides(fill = FALSE) +
  theme(panel.grid = element_blank()) +
  labs(y = "Number of targets", x = "Predicted Tractability Classification", subtitle = "TargetDB Target Tractability Assessment")

tract_fig <- tractability_plot / tdb_plot + plot_annotation(tag_levels = 'a')
ggsave(tract_fig, width = 7, height = 8, dpi = 300, file = paste0(default_save_dir, "fig3.png"))
ggsave(tract_fig, width = 7, height = 8, file = paste0(default_save_dir, "fig3.pdf"))

###########################################################################
###########################################################################
###                                                                     ###
###                       SECTION 1:                                    ###
###                       FUNCTIONS                                     ###
###                                                                     ###
###########################################################################
###########################################################################


## Processing data from Open Targets 
get_lead_variants <- function(overwrite = FALSE, save_dir = default_save_dir){
  if (overwrite){
    library(arrow)
    
    min_n_cases = 2000    ## min n cases for case control studies
    min_n_initial = 2000  ## min sample size for continuous traits
    
    ## All lead variant data for ARDs
    ard_leads_all <- arrow::read_parquet("data/targetage/ard_leads.parquet/part-00000-847cd340-cf61-47da-a25d-9d36d195ef39-c000.snappy.parquet")
    
    ## Filter out small GWAS studies
    ard_leads <- filter(ard_leads_all, (n_cases >= min_n_cases) | (is.na(n_cases) & n_initial >= min_n_cases))
    
    ## Save output 
    save(ard_leads, file = paste0(save_dir, "ard_leads_filtered.Rda"))
  }
  else {
    # Load preprepared leads 
    load(paste0(save_dir, "ard_leads_filtered.Rda"))
  }
  return(ard_leads)
}

get_overlaps <- function(overwrite = FALSE, save_dir = default_save_dir){
  if (overwrite){
    library(arrow)
    library(data.table) 
    
    min_n_cases = 2000    ## min n cases for case control studies
    min_n_initial = 2000  ## min sample size for continuous traits
    
    ### Colocalisation data 
    coloc <- arrow::read_parquet("data/targetage/coloc_ard_leads.parquet/part-00000-7408f34e-85ab-43bb-a12f-f5e78b7d314d-c000.snappy.parquet", as_tibble=TRUE)
    
    ### Overlap data
    overlap <- arrow::read_parquet("data/targetage/overlap_ard_leads.parquet/part-00000-968f5e84-1221-4b35-8a99-69fccda34afb-c000.snappy.parquet",as_data_frame = TRUE)
    
    ### Filter out small GWAS studies
    overlap <- overlap %>% 
      filter((n_cases >= min_n_cases) | (is.na(n_cases) & n_initial >= min_n_cases)) %>%
      filter((right_n_cases >= min_n_cases) | (is.na(right_n_cases) & right_n_initial >= min_n_cases)) %>%
      filter(!is.na(right_is_lead)) %>% ## one variant is NA for this: GCST006288 11_86942946_G_A
      as.data.table() 
    
    ## All the colocalisation where the right hand variant is the lead variant (so we have one node per associated locus)
    lead_coloc <- coloc %>% filter(right_is_lead == TRUE) 
    
    ## All the colocalisation hits where the right hand study is an xQTL study, for later
    qtl_coloc <- coloc %>% filter(right_type %in% c("eqtl", "pqtl"))
    
    ## All the hits within the ARDs (i.e. both studies are for ARDs)
    coloc_within <- lead_coloc %>% 
      filter(!is.na(right_morbidity)) %>% 
      filter((n_cases >= min_n_cases) | (is.na(n_cases) & n_initial >= min_n_cases)) %>%
      filter((right_n_cases >= min_n_cases) | (is.na(right_n_cases) & right_n_initial >= min_n_cases))
    overlap_within <- overlap[!is.na(right_morbidity)] %>% as.data.frame()
    
    ## All the hits with traits other than our ARDs of interest
    coloc_without <- lead_coloc %>% 
      filter(is.na(right_morbidity)) %>% 
      filter((n_cases >= min_n_cases) | (is.na(n_cases) & n_initial >= min_n_cases)) %>%
      filter((right_n_cases >= min_n_cases) | (is.na(right_n_cases) & right_n_initial >= min_n_cases))
    overlap_without <- overlap[is.na(right_morbidity)] %>% as.data.frame()
    
    ## overlaps
    overlaps <- list("overlap_within" = overlap_within, "coloc_within" = coloc_within)
    overlaps_without <- list("overlap_without" = overlap_without, "coloc_without" = coloc_without)
    
    ## save output 
    save(overlaps, file = paste0(save_dir, "overlaps_within.Rda"))  
    save(overlaps_without, file = paste0(save_dir, "overlaps_without.Rda"))
    save(qtl_coloc, file = paste0(save_dir, "qtl_coloc"))
  }
  else {
    load(paste0(save_dir, "overlaps_within.Rda"))
  }
  return (overlaps)
}

get_associations_annotations <- function(overwrite = FALSE, save_dir = default_save_dir, diseases){
  if (overwrite){
    library(arrow)
    ## Target details from Open Targets
    targets <- arrow::read_parquet("data/targetage/ard_annotations.parquet/part-00000-2d7e85c7-0644-4728-9921-dbdf119c427d-c000.snappy.parquet", as_tibble=TRUE)
    
    ## Associations
    associations <- arrow::read_parquet("data/targetage/ard_associations.parquet/part-00000-80ca0def-de33-4d7c-a012-dcbe8813a3d4-c000.snappy.parquet", as_tibble=TRUE) %>% 
      select(morbidity, targetId, overallDatatypeHarmonicVector) %>% 
      unnest(overallDatatypeHarmonicVector)
    
    ## genetic associations
    genetic_associations <- associations %>% 
      filter(datatypeId=="genetic_association", morbidity %in% diseases) %>% 
      select(morbidity, targetId, datatypeHarmonicScore) %>% 
      mutate(datatypeHarmonicScore = round(datatypeHarmonicScore, 2)) %>% 
      pivot_wider(names_from = "morbidity", values_from = "datatypeHarmonicScore") 
    
    literature_counts <- associations %>% 
      filter(datatypeId=="literature") %>% 
      select(morbidity, targetId, literatureEvidenceCount = datatypeEvidenceCount)
    
    ## Which morbidity has the largest number of literature associations
    max_literature_counts <- literature_counts %>% 
      group_by(targetId) %>% 
      slice(which.max(literatureEvidenceCount)) %>% 
      rename(maxLiteratureEvidencePhenotype = morbidity, maxLiteratureEvidenceCount = literatureEvidenceCount)
    
    ## go terms mapped to hallmarks of ageing 
    go_terms <- read.csv("data/full_goterm_list.csv")
    
    target_go <- targets %>% 
      select(targetId, go) %>% 
      unnest_wider(go) %>% 
      unnest(c(id, value), keep_empty = TRUE) %>%
      rename(goId = id) %>%
      left_join(go_terms) 
    
    save(target_go, file = paste0(save_dir, "targets_all_go.Rda"))
    
    ## Read in data from GenAge, CellAge
    load(paste0(save_dir,"genage.Rda"))
    load(paste0(save_dir,"cellage.Rda"))
    
    # number of protein-coding genes in OT (need to find out)
    in_TA = 646
    in_GA = 307
    in_TA_GA = 19
    n_OT = 25000
    
    genage_dat <- data.frame(
      "in_TA" = c(in_TA_GA, in_TA - in_TA_GA),
      "not_in_TA" = c(in_GA - in_TA_GA, n_OT - in_TA - in_GA + in_TA_GA),
      row.names = c("in_GA", "not_in_GA"),
      stringsAsFactors = FALSE
    )
    fisher.test(genage_dat)$p.value 
    
    in_CA = 279
    in_TA_CA = 21
    
    cellage_dat <- data.frame(
      "in_TA" = c(in_TA_CA, in_TA - in_TA_CA),
      "not_in_TA" = c(in_CA - in_TA_CA, n_OT - in_TA - in_CA + in_TA_CA),
      row.names = c("in_CA", "not_in_CA"),
      stringsAsFactors = FALSE
    )
    fisher.test(genage_dat)$p.value 
         
    ## Target annotations of interest to us
    target_annotations <- targets %>% 
      select(targetId, targetSymbol, targetName, bioType, hgncId, chemicalProbes, symbolSynonyms, tractability, safety, tep) %>% 
      jsonlite::flatten() %>% 
      left_join(max_literature_counts) %>%
      left_join(genetic_associations) %>%
      left_join(cellage %>% group_nest(targetSymbol, .key = "CellAge")) %>% 
      left_join(genage %>% group_nest(targetSymbol, .key = "GenAge")) %>%
      left_join(target_go %>% 
                  filter(!is.na(goHallmarkId)) %>% 
                  group_nest(targetId, .key = "ageingHallmarks"))
    
    save(target_annotations, file = paste0(save_dir, "target_annotations.Rda"))
    
  } else {
    load(file = paste0(save_dir, "target_annotations.Rda"))
  }
  
  return(target_annotations)
}

## Just to get nice colours for the graphs
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

## Prepares edges and nodes for the graph for a given ARD(s)
gwas_graph <- function(curr_morbidities, curr_variants, curr_coloc, curr_overlap, min_jaccard = 0, plot = TRUE) {
  library(visNetwork)
  library(igraph)
  
  g <- list()
  
  # Different colours for ARDs
  g$colours <-
    data.frame(color = gg_color_hue(length(curr_morbidities)), morbidity = curr_morbidities)
  # Get variants and overlaps just for the ARD(s) of interest
  m_variants <- curr_variants %>% filter(morbidity %in% curr_morbidities)
  m_overlap <-
    curr_overlap %>% filter(morbidity %in% curr_morbidities &
                              right_morbidity %in% curr_morbidities)
  m_coloc <-
    curr_coloc %>% filter(morbidity %in% curr_morbidities &
                            right_morbidity %in% curr_morbidities)
  
  # Edges: where two studyId+variantId have at least one overlapping tag variant
  g$edges <-
    m_overlap %>%
    # minimum overlap criteria 
    filter((LR_overlap / (right_distinct + lead_distinct + LR_overlap)) >= min_jaccard) %>% 
    mutate(
      from = paste(studyId, lead_variantId, sep = "."),
      to = paste(right_studyId, right_variantId, sep = ".")
    ) %>%
    # remove self-edges
    filter(from != to) %>%
    # dashed edge if both have sum stats (therefore will have coloc data)
    # we will use these edges for visualisation but not for calculating clusters 
    mutate(dashes = ifelse(has_sumstats == TRUE & right_has_sumstats, TRUE, FALSE)) %>%
    select(from, to, morbidity, dashes) %>%
    left_join(g$colours, by = "morbidity") 
  
  # Edges: where two studyId+variantId are colocalised (edges in black)
  g$coloc_edges <-
    m_coloc %>%
    mutate(
      from = paste(studyId, lead_variantId, sep = "."),
      to = paste(right_studyId, right_variantId, sep = ".")
    ) %>%
    filter(from != to) %>%
    select(from, to, morbidity) %>%
    mutate(color = "black", dashes = FALSE)
  
  # Each node is a studyId+variantId
  # Make a label with the variant, reported trait, and n_cases
  g$nodes <- m_variants %>%
    mutate(id = paste(studyId, lead_variantId, sep = ".")) %>%
    select(id,
           studyId,
           lead_variantId,
           morbidity,
           has_sumstats,
           n_cases,
           trait_reported,
           direction) %>%
    unique() %>%
    left_join(g$colours, by = "morbidity") %>%
    mutate(label = paste0(lead_variantId, "\n", trait_reported, " (", n_cases, ") ", direction)) %>%
    group_by(id) %>%
    mutate(n=length(morbidity)) %>% 
    mutate(morbidity = paste0(morbidity, collapse="+")) %>%
    mutate(color = ifelse(n==1, color, "black")) %>%
    slice(1) %>%
    select(-n)
  
  if (plot){
    # Make the plot
    g$plot <-
      visNetwork(
        g$nodes,
        bind_rows(g$edges, g$coloc_edges),
        main = curr_morbidities,
        height = "1000px",
        width = "100%"
      ) %>% 
      visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
      visIgraphLayout()
  }
  return(g)
}

## Makes the graph and counts the number of disconnected clusters and communities within a disconnected cluster
count_communities <- function(curr_morbidities, detect_subgraph_communities = FALSE, ...) {
  ## get relevant edges and nodes
  g <- gwas_graph(curr_morbidities, ...)
  if (nrow(g$edges) > 0) {
    ## make an igraph object
    g$ig <- graph_from_data_frame(
      bind_rows(
        subset(g$edges, dashes == FALSE),  ## overlap edges one or both studies don't have sumstats
        g$coloc_edges),                    ## colocalisation edges otherwise (both have sumstats)
      directed = FALSE)
    ## get clusters
    clusters <- clusters(g$ig)
    overall_communities <- cluster_louvain(g$ig, weights = NA)
    V(g$ig)$cluster <- clusters$membership
    V(g$ig)$overall_community <- overall_communities$membership
    
    # detect communities in connected component subgraphs
    if (detect_subgraph_communities == TRUE){
      V(g$ig)$community = 0
      for (component in unique(clusters$membership)) {
        g_component_subgraph <-
          induced_subgraph(g$ig, which(clusters$membership == component))
        m <- cluster_louvain(g_component_subgraph)
        V(g$ig)[V(g_component_subgraph)$name]$community = m$membership + max(V(g$ig)$community)
      }
    }
    
    ## add cluster ID, cluster size, and n_morbidities to nodes dataframe
    n <- get.data.frame(g$ig, what = "vertices")
    g$cn <- g$nodes %>%
      left_join(n, by = c("id" = "name")) %>%
      group_by(cluster) %>%
      mutate(c_size = ifelse(is.na(cluster), NA, length(id)),
             n_morbidities = length(unique(morbidity)))
    g$clusters <-
      g$cn %>% 
      group_by(cluster, n_morbidities, morbidity) %>% 
      count() %>% 
      pivot_wider(names_from = "morbidity", values_from = n)
    if (length(curr_morbidities) == 2) {
      ## remove nodes with both morbidity
      ## get number of nodes for each morbidity in each cluster
      tmp <-
        g$cn %>% 
        filter(str_detect(morbidity, "\\+", negate = TRUE)) %>% 
        group_by(cluster, n_morbidities, morbidity) %>% 
        count()
      ## how many nodes aren't in a cluster i.e. don't overlap with any other nodes (i.e. unreplicated) for each morbidity
      singles <-
        tmp %>% 
        ungroup() %>% 
        filter(is.na(cluster)) %>% 
        select(morbidity, single = n)
      ## for clusters, how many consist of just one morbidity (distinct) and how many contain both (shared)
      shared <- tmp %>%
        # remove single nodes (accounted for in `singles` above)
        filter(!is.na(cluster)) %>% 
        select(-n) %>%
        # clusters with more than 2 morbidities include nodes from GWAS studies involving both morbidity A and B together, so treat these as 2 (i.e. `shared`)
        mutate(n_morbidities = ifelse(n_morbidities > 1, 2, n_morbidities)) %>% 
        group_by(morbidity, n_morbidities) %>%
        count() %>%
        mutate(n_morbidities = recode(n_morbidities, `1` = "distinct", `2` = "shared")) %>% 
        pivot_wider(names_from = n_morbidities, values_from = n) %>%
        full_join(singles)
      if (nrow(shared) == 0) {
        shared <-
          data.frame(morbidity = curr_morbidities) %>% 
          mutate(distinct = 0, shared = 0) %>% left_join(singles)
      }
      g$overlap <- shared %>%
        ungroup() %>%
        rename(morbidity_A = morbidity) %>%
        mutate(morbidity_B = rev(.$morbidity_A)) %>%
        inner_join(
          shared %>% 
            rename(morbidity_B = morbidity),
          by = "morbidity_B",
          suffix = c("_A", "_B")
        )
    }
  } else {
    if (length(curr_morbidities) == 2) {
      singles <-
        g$nodes %>% 
        group_by(morbidity) %>% 
        count(name = "single")
      shared <-
        data.frame(morbidity = curr_morbidities) %>% 
        mutate(distinct = 0, shared = 0) %>% 
        left_join(singles)
      g$overlap <- shared %>%
        rename(morbidity_A = morbidity) %>%
        mutate(morbidity_B = rev(.$morbidity_A)) %>%
        inner_join(
          shared %>% 
            rename(morbidity_B = morbidity),
          by = "morbidity_B",
          suffix = c("_A", "_B")
        )
    }
  }
  return(g)
}

## Create a graph with all morbidities
create_targetage_graph <- function(overwrite = FALSE, save_dir, ...){
  library(visNetwork)
  library(igraph)
  
  if (overwrite) {
    g_all <- count_communities(d, ard_leads, overlaps$coloc_within, overlaps$overlap_within, min_jaccard = 0, detect_subgraph_communities = TRUE, plot = FALSE)
    save(g_all, file = paste0(save_dir, "graph_all_morbidities.Rda"))
  } else {
    load(paste0(save_dir, "graph_all_morbidities.Rda"))
  }
  return(g_all)
}

get_pairwise_overlaps <- function(g){
  ## get morbidities in each cluster (excluding nodes with multiple morbidities)
  tmp_mat <- g$clusters %>% 
    ungroup() %>% 
    filter(!is.na(cluster)) %>% ## unclustered nodes
    select(-cluster, -n_morbidities, -contains("+"))  %>% 
    as.matrix()
  morbidities_per_cluster <- lapply(seq(nrow(tmp_mat)), function(i) names(which(tmp_mat[i,] > 0)))
  # get all possible pairwise combinations 
  pairwise_combinations <- gtools::combinations(n = ncol(tmp_mat), r = 2, colnames(tmp_mat))
  # count how many times each pair appears together in a cluster 
  count <- apply(pairwise_combinations, 1, function(x){
    shared_clusters <- lapply(x, function(y) sapply(morbidities_per_cluster, `%in%`, x = y))
    sum(Reduce(`&`, shared_clusters))
  })
  #n_pairwise_shared_clusters <- cbind(data.frame(all.pairs, count))[count > 0,] 
  n_pairwise_shared_clusters <- cbind(data.frame(all.pairs, count)) 
  
  # join with total number of clusters for each morbidity
  overall_n_clusters <- g$clusters %>% 
    ungroup() %>% 
    filter(!is.na(cluster)) %>% 
    select(-n_morbidities, -contains("+")) %>% 
    pivot_longer(-cluster) %>% 
    group_by(name) %>% 
    summarise(X_overall = sum(value>0, na.rm = TRUE)) %>% 
    mutate(Y_overall = X_overall)
  pairwise_shared <- n_pairwise_shared_clusters %>% 
    data.frame() %>% 
    left_join(select(overall_n_clusters, name, X_overall) , by = c("X1" = "name")) %>% 
    left_join(select(overall_n_clusters, name, Y_overall), by = c("X2" = "name")) %>% 
    rename(X = X1, Y = X2) %>%
    mutate(total = nrow(tmp_mat))
  
  ## calculate significance of overlap
  n_tests <- length(pairwise_combinations)
  
  ## Do I correct using the total number of pairs (possible tests) or the number of tests I actually perform? I guess getting 
  ft_all <- pairwise_shared %>% 
    rowwise() %>% 
    mutate(p.value = get_fishers_exact(overlap = count, morb_x_overall = X_overall, morb_y_overall = Y_overall, total = total)) %>%
    ungroup() %>% 
    mutate(p.adjust = p.adjust(p.value, method = "BH", n = n_tests))   # n tests
  # mutate(p.adjust = p.adjust(p.value, method = "BH", n = nrow(.)))   # n possible pairs
  
  return(ft_all)
}

## One-sided fishers exact test for overlapping clusters between 
## pairwise combination of morbidities
get_fishers_exact <- function(overlap, morb_x_overall, morb_y_overall, total, alternative = "greater"){
  dat <- data.frame(
    "morb_x" = c(overlap, morb_x_overall - overlap),
    "not_morb_x" = c(morb_y_overall - overlap, total - morb_x_overall - morb_y_overall + 3),
    row.names = c("morb_y", "not_morb_y"),
    stringsAsFactors = FALSE
  )
  fisher.test(dat, alternative = alternative)$p.value 
}


## Calculate overlap between each pair of diseases in a graph
calculate_pairwise_overlap <- function(g){
  unreplicated_associations <- g_all$cn %>% filter(is.na(cluster)) %>% count(morbidity)
}

## Initialise OT Genetics Portal API
initialise_api <- function(){
  cli <- GraphqlClient$new(
    url = "https://genetics-api.opentargets.io/graphql"
  )
  return(cli)
}

## GraphQL Query to retrive L2G
initialise_queries <- function(){
  qry <- Query$new()
  qry$query('l2g_query', 'query l2gQuery($studyId: String!, $variantId: String!){
  studyLocus2GeneTable(studyId: $studyId, variantId: $variantId){
    rows {
      gene {
        id
        symbol
      }
      hasColoc
      yProbaModel
      yProbaDistance
      yProbaInteraction
      yProbaMolecularQTL
      yProbaPathogenicity
      distanceToLocus
    }
  }
  variantInfo(variantId: $variantId){
      mostSevereConsequence
  }
}')
  return(qry)
}

## Do the API call
fetch_l2g <- function(df, variables){
  result <- fromJSON(cli$exec(qry$queries$l2g_query, variables, flatten = TRUE))$data
  
  l2g_result <- result$studyLocus2GeneTable %>% bind_cols(result$variantInfo) %>% bind_cols(df) 
  return(l2g_result)
}

## Retrieve L2G scores from Open Targets Genetics
## Warning: This will take a while!
get_L2G <- function(overwrite = default_overwrite, save_dir = default_save_dir){

  if(overwrite) {
    library(ghql)
    library(jsonlite)
    
    # studies + variants
    studies_variants <- g_all$cn %>% filter(n_morbidities > 1) %>% ungroup() %>% select(studyId, lead_variantId) %>% unique()

    cli <- initialise_api()
    qry <- initialise_queries()
    
    ## split the data frame into smaller chunks (1000 rows)
    ## I don't really know if we need to do this but just in case
    ## I don't want to break the server with too many successive calls
    
    n <- 100
    nr <- nrow(studies_variants)
    studies_variants_split <-
      studies_variants %>% 
      split(., rep(1:ceiling(nr / n), each = n, length.out = nr))
    
    # Somewhere to hold the results
    l2g_all <- vector(mode = "list", length = length(studies_variants_split))
    
    ## Do the first chunk on its own to check 
    l2g_all[[1]] <-  studies_variants_split[[1]]  %>% 
      group_by(studyId, lead_variantId) %>% 
      group_split() %>% 
      ## API call for each studyID + variantID
      purrr::map(~fetch_l2g(df = ., variables = list(studyId = .$studyId, variantId = .$lead_variantId))) %>%
      bind_rows() 
    
    ## Do all the other chunks
    for (i in seq(2, length(l2g_all))){
      print(i)
      l2g_all[[i]] <- studies_variants_split[[i]] %>% 
        group_by(studyId, lead_variantId) %>% 
        group_split() %>% 
        ## API call for each studyID + variantID
        purrr::map(~fetch_l2g(df = ., variables = list(studyId = .$studyId, variantId = .$lead_variantId))) %>%
        bind_rows() 
      Sys.sleep(3)
    }
    
    ## Collapse into one dataframe, remove low confidence (L2G<0.05)
    l2g_all_joined <- l2g_all %>% bind_rows() %>% jsonlite::flatten() %>% filter(yProbaModel > 0.05)
    
    ## Use more up-to-date gene symbol from the Platform 
    l2g_all_joined <- l2g_all_joined %>% left_join(targets_annotations %>% select(gene.id = targetId, targetSymbol))
    
    ## Save output
    save(l2g_all_joined, file=paste0(save_dir, "ltg_all.Rda"))
    
  } else {
    load(paste0(save_dir, "ltg_all.Rda"))
  }
  
  return(ltg_all_joined)
  
}

get_L2G_coloc <- function(overwrite = default_overwrite, save_dir = default_save_dir){
  if (overwrite){
    load(paste0(save_dir, "qtl_coloc"))
    qtl_to_join <- qtl_coloc %>%
      select(studyId, lead_variantId, right_type, right_studyId, coloc_h4, starts_with("lead_var_right_"), right_phenotype, right_bio_feature, right_gene_id) %>%
      mutate(direction = ifelse(lead_var_right_study_beta > 0, "increased", "decreased")) %>% 
      mutate(trait = ifelse(right_type == "eqtl", "expression", ifelse(right_type == "pqtl", "abundance", NA))) %>% 
      unite(molecular_trait, c(direction, trait), sep = " ", remove = TRUE)
    
    l2g_qtl_coloc <- l2g_all_joined %>% 
      filter(hasColoc == TRUE) %>% 
      left_join(qtl_to_join, c("studyId", "lead_variantId", "gene.id" = "right_gene_id"))
    
    save(l2g_qtl_coloc, file = paste0(save_dir,"l2g_qtl_coloc.Rda"))
  } else {
    load(paste0(save_dir,"l2g_qtl_coloc.Rda"))
  }
  return(l2g_qtl_coloc)
}

## make the heatmap figure showing overlaps between phenotypes
plot_shared_heatmap <- function(pairwise_overlaps, overwrite = default_overwrite, save_dir = default_save_dir){
  library(ggdendro)
  
  textcol <- "grey40"
  
  df <- pairwise_overlaps %>% 
    mutate(proportion_X = count/X_overall,
           face = ifelse(p.value >= 0.05, "plain", ifelse(p.adjust < 0.05, "bold.italic", "bold"))) %>%
    mutate_at((c("X", "Y")), str_replace_all, pattern ="chronic obstructive pulmonary disease", replacement = "COPD") 
  
  
  # make the bottom half of the matrix
  swapped_df <- 
    df %>% 
    rename(tmpX = Y,
           tmpY = X,
           tmpX_overall = Y_overall, 
           tmpY_overall = X_overall) %>%
    rename_all(~stringr::str_replace(.,"^tmp","")) %>%
    mutate(proportion_X = count/X_overall)
  
  full_df <- df %>% 
    bind_rows(swapped_df)
  
  # ## clustering
  # pre_mat <- full_df %>% 
  #   select(X, Y, proportion_X) %>%
  #   pivot_wider(names_from = Y, values_from = proportion_X) %>%
  #   replace(is.na(.), 0)
  # 
  # mat <- as.matrix(pre_mat[,-1])
  # rownames(mat) <- pre_mat$X
  # shared.dendro <- as.dendrogram(hclust(d = dist(x = mat)))
  # 
  # # Create dendro
  # dendro.plot <- ggdendrogram(data = shared.dendro, rotate = TRUE)
  # 
  # ## reorder to the clustered order
  # dendro.order <- order.dendrogram(shared.dendro)
  # 
  # full_df$Y <- factor(x = full_df$Y,
  #                levels = unique(pre_mat$X)[dendro.order],
  #                ordered = TRUE)
  # full_df$X <- factor(x = full_df$X,
  #                levels = unique(pre_mat$X)[dendro.order],
  #                ordered = TRUE)
  
  
  ## heatmap 
  shared_heatmap <- full_df %>%
    arrange(X_overall) %>% 
    mutate(X=factor(X, levels=unique(X)),
           Y=factor(Y, levels = unique(X))) %>%   
    mutate(proportion_X = ifelse(proportion_X == 0, NA, proportion_X)) %>%
    # mutate_at((c("morbidity_A", "morbidity_B")), str_replace_all, pattern ="chronic obstructive pulmonary disease", replacement = "COPD") %>%
    ggplot(., aes(x=X, y = Y, fill = proportion_X)) +
    geom_tile(colour = textcol) +
    geom_tile(data = subset(full_df, p.value < 0.05), colour = textcol, size = 0.65) +
    #   geom_text(aes(label = count, colour = proportion_X > 0.45, fontface = face), size = 2.1) + 
    geom_text(aes(label = ifelse(p.adjust < 0.05, paste0(count, "*"), count), colour = proportion_X > 0.45), size = 2.1) + 
    labs(subtitle = "Independent genetic associations shared between\nage-related diseases and traits") +
    coord_fixed() +
    theme_bw(base_size = 10) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_color_manual(values = c("black", "white"), breaks = c(FALSE, TRUE)) +
    scale_fill_gradient(
      "proportion",
      low = "#CBDEF0",
      high = "#08306B",
      na.value = "white",
      guide = guide_colorbar(frame.colour = "black", frame.linewidth = 0.8)
    ) +
    guides(color = FALSE) + 
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(
        colour = textcol,
        angle = 45,
        hjust = 1
      ),
      axis.text.y = element_text(vjust = 0.2,
                                 colour = textcol),
      axis.ticks = element_line(size = 0.4),
      axis.title = element_blank(),
      plot.background = element_blank(),
      plot.margin = margin(0.7, 0.4, 0.1, 0.2, "cm"),
      plot.title = element_text(
        colour = textcol,
        hjust = 0,
        size = 14,
        face = "bold"
      ),
      panel.border = element_rect(colour = "black",
                                  size = 0.8)
    )
  
  if (overwrite) {
    ggsave(shared_heatmap, file = paste0(save_dir, "Fig1.png"), width = 8, height = 7.5, dpi = 300)
    ggsave(shared_heatmap, file = paste0(save_dir, "Fig1.pdf"), width = 8, height = 7.5)
  }
  
  return(shared_heatmap)
}







wide_associations <- associations %>% 
  select(-weight, -datatypeEvidenceCount) %>% 
  pivot_wider(names_from = datatypeId, values_from = datatypeHarmonicScore) %>% 
  left_join(associations %>% 
              filter(datatypeId=="literature") %>% 
              select(targetId, morbidity, literatureEvidenceCount = datatypeEvidenceCount))
  




