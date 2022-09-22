
###########################################################################
###########################################################################
###                                                                     ###
###                       SECTION 1:                                    ###
###               Process Input Data FUNCTIONS                          ###
###                                                                     ###
###########################################################################
###########################################################################

## Processing data from Open Targets 
get_lead_variants <- function(overwrite = FALSE, input_file, save_dir = default_save_dir){
  if (overwrite){
    library(arrow)
    
    min_n_cases = 2000    ## min n cases for case control studies
    min_n_initial = 2000  ## min sample size for continuous traits
    
    ## All lead variant data for ARDs
    ard_leads_all <- arrow::read_parquet(input_file)
    
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

## Calculate overlap between each pair of diseases in a graph
get_overlaps <- function(overwrite = FALSE, coloc_input_file, overlap_input_file, save_dir = default_save_dir){
  if (overwrite){
    library(arrow)
    library(data.table) 
    
    min_n_cases = 2000    ## min n cases for case control studies
    min_n_initial = 2000  ## min sample size for continuous traits
    
    ### Colocalisation data 
    coloc <- arrow::read_parquet(coloc_input_file, as_tibble=TRUE)
    
    ### Overlap data
    overlap <- arrow::read_parquet(overlap_input_file, as_data_frame = TRUE)
    
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

get_go_terms <- function(targets, go_input_file = go_input_file, save_dir = default_save_dir, overwrite = default_overwrite){
  ## go terms mapped to hallmarks of ageing 
  go_terms <- read.csv(go_input_file)
  
  target_go <- targets %>% 
    select(targetId, go) %>% 
    unnest_wider(go) %>% 
    unnest(c(id), keep_empty = TRUE) %>%
    select(targetId, id) %>% 
    rename(goId = id) %>%
    left_join(go_terms) 
  if (overwrite){
    save(target_go, file = paste0(save_dir, "targets_all_go.Rda"))
  } else {
    load(paste0(save_dir, "targets_all_go.Rda"))
  }
  return(target_go)
}

get_associations_annotations <- function(overwrite = FALSE, annotations_input_file, associations_input_file, go_input_file, save_dir = default_save_dir, diseases){
  if (overwrite){
    library(arrow)
    
    ## Read in data from GenAge, CellAge
    load(paste0(save_dir,"genage.Rda"))
    load(paste0(save_dir,"cellage.Rda"))
    
    ## Target details from Open Targets
    targets <- arrow::read_parquet(annotations_input_file, as_tibble=TRUE)
    
    ## Get ids of the top level diseases to get association scores
    top_level_disease_ids <- diseases %>% 
      filter(diseaseId == specificDiseaseId) %>% 
      pull(diseaseId)
    
    ## Associations
    associations <- arrow::read_parquet(associations_input_file, as_tibble=TRUE) %>% 
      filter(diseaseId %in% top_level_disease_ids) %>% 
      mutate_if(is.numeric, round, 2) 
    
    ## Which morbidity has the largest number of literature associations
    max_literature_counts <- associations %>% 
      select(targetId, morbidity, literatureCount) %>% 
      group_by(targetId) %>% 
      slice(which.max(literatureCount)) %>% 
      rename(maxLiteratureEvidencePhenotype = morbidity, maxLiteratureEvidenceCount = literatureCount)
    
    ## genetic associations
    genetic_associations <- associations %>%
      select(targetId, morbidity, genetic_association) %>%
      pivot_wider(names_from = morbidity, values_from = genetic_association) %>%
      left_join(max_literature_counts, by = "targetId")
    
    ## HGNC and PDB ids
    xref_ids <- targets %>% 
      select(targetId, dbXrefs) %>% 
      unnest(dbXrefs) %>% 
      filter(source %in% c("HGNC", "PDB")) %>% 
      group_by(targetId, source) %>% 
      summarise(
        source = toString(unique(source)),
        id = toString(unique(id))
      ) %>% pivot_wider(names_from = source, values_from = id)
    
    ## GO terms
    target_go <- get_go_terms(targets, go_input_file, save_dir, overwrite)
    
    ## Target annotations of interest to us
    target_annotations <- targets %>% 
      select(targetId, targetSymbol, functionDescriptions, targetName, biotype, targetClass, dbXrefs, chemicalProbes, symbolSynonyms, tractability, safetyLiabilities, tep) %>% 
      jsonlite::flatten() %>% 
      left_join(max_literature_counts) %>%
      left_join(genetic_associations) %>%
      left_join(cellage %>% group_nest(targetSymbol, .key = "CellAge")) %>% 
      left_join(genage %>% group_nest(targetSymbol, .key = "GenAge")) %>%
      left_join(target_go %>% 
                  filter(!is.na(goHallmarkId)) %>% 
                  group_nest(targetId, .key = "ageingHallmarks")) %>%
      left_join(xref_ids, by = "targetId")
    
    save(target_annotations, file = paste0(save_dir, "target_annotations.Rda"))
    
  } else {
    load(file = paste0(save_dir, "target_annotations.Rda"))
  }
  return(target_annotations)
}

## Counts for each hallmark
get_hallmark_counts <- function(targetage_hallmarks, save_dir = default_save_dir, overwrite = default_overwrite){
  # High level GO Terms
  go_terms_high <- read.csv("data/goterm_list.csv", col.names = c("goHallmarkId", "goHallmark", "goHallmarkTerm"))
  # including descendant terms
  go_terms <- read.csv("data/full_goterm_list.csv")
  # the number of terms under each high level term
  total_terms <- go_terms %>% group_by(goHallmarkId, goHallmark, goHallmarkTerm) %>% summarise(n_terms = length(goTerm))
  # the number of targets in each high level GO term 
  term_targets <- targetage_hallmarks %>% select(targetId, goHallmarkId) %>% unique() %>% group_by(goHallmarkId) %>% summarise(n_targets = length(targetId))
  # the number of targets in each Hallmark
  hallmark_targets <- targetage_hallmarks %>% select(targetId, goHallmark) %>% unique() %>% group_by(goHallmark) %>% summarise(total_targets = length(targetId))
  hallmarks_terms_targets <- go_terms_high %>% 
    left_join(total_terms) %>% 
    left_join(term_targets) %>% 
    left_join(hallmark_targets) %>% 
    select(goHallmark,everything())
  
  if (overwrite){
    write.csv(hallmarks_terms_targets, file = paste0(save_dir, "targetage_hallmark_counts.csv"), row.names = FALSE, quote = FALSE)
  }
  return (hallmarks_terms_targets)
}

## Warning: the following section uses ClusterProfiler which masks many important tidyverse functions
## so don't run it unless you have to
convert_ensembl_entrez <- function(ids){
  clusterProfiler::bitr(ids, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID
}

get_entrez_gene_sets <- function(overwrite = default_overwrite, save_dir = default_save_dir){
  if (overwrite){
    load(paste0(default_save_dir,"genage.Rda"))
    load(paste0(default_save_dir,"cellage.Rda"))
    
    targetage_entrez <- convert_ensembl_entrez(unique(targetage_annotations$targetId))
    gene_sets_entrez <- list(TargetAge = targetage_entrez,
                             GenAgeHuman = genage$entrez.gene.id %>% unique(),
                             CellAge = cellage$CellAge.entrezid %>% unique()
    )
    save(gene_sets_entrez, file = paste0(save_dir, "gene_sets_entrez.Rda"))
  } else {
    load(paste0(save_dir, "gene_sets_entrez.Rda"))
  }
  return (gene_sets_entrez)
}


get_gene_sets <- function(overwrite = default_overwrite, save_dir = default_save_dir){
  if (overwrite){
    load(paste0(default_save_dir,"genage.Rda"))
    load(paste0(default_save_dir,"cellage.Rda"))
    hallmarks_genes <- read.csv("data/full_goterm_genes.csv")
    
    gene_sets <- list(TargetAge = targetage_annotations$targetSymbol %>% unique(),
                      CellAge = cellage$targetSymbol %>% unique(),
                      GenAge = genage$targetSymbol %>% unique(),
                      Hallmarks = hallmarks_genes$SYMBOL %>% unique(),
                      ARDs = target_annotations$targetSymbol %>% unique()
    )
    
    save(gene_sets, file = paste0(save_dir, "gene_sets.Rda"))
  } else {
    load(paste0(save_dir, "gene_sets.Rda"))
  }
  return(gene_sets)
}

## Get and unnest an annotation 
get_subset_annotations <- function(annotations, genes = targetage, field){
  annotations %>% 
    filter(targetId %in% genes) %>% 
    dplyr::select(targetId, targetSymbol, all_of(field)) %>%
    hoist(field) %>% 
    unnest(field)
}


###########################################################################
###########################################################################
###                                                                     ###
###                       SECTION 2:                                    ###
###                    Graph analysis                                   ###
###                                                                     ###
###########################################################################
###########################################################################

## Just to get nice colours for the graphs
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

## Prepares edges and nodes for the graph for a given ARD(s)
gwas_graph <- function(curr_morbidities, curr_variants, curr_coloc, curr_overlap, min_jaccard = 0, plot = TRUE, id = NULL) {
  library(visNetwork)
  library(igraph)
  
  g <- list()
  
  g$id <- id
  
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
    g_all <- count_communities(d, ard_leads, overlaps$coloc_within, overlaps$overlap_within, min_jaccard = 0, detect_subgraph_communities = TRUE, plot = FALSE, id = "targetage")
    save(g_all, file = paste0(save_dir, "graph_all_morbidities.Rda"))
  } else {
    load(paste0(save_dir, "graph_all_morbidities.Rda"))
  }
  return(g_all)
}

# Get information on how many clusters have multiple communities 
get_n_multicommunity_clusters <- function(g_all){
  clustered <- g_all$cn %>% filter(!is.na(cluster))
  all <- bind_rows(clustered %>% mutate(set = "all"),
                   clustered %>% filter(n_morbidities > 1) %>% mutate(set = "multimorbidity"))
  all %>% 
    select(set, cluster, community) %>%
    unique() %>% 
    count(set, cluster) %>%
    group_by(set) %>% 
    summarise(total = length(cluster),
              multicommunity = sum(n>1),
              two_communities = sum(n==2),
              single_community = sum(n==1)) %>%
    mutate_at(.vars = c("multicommunity", "two_communities", "single_community"), 
              .funs = c(p = ~./total))
}

## Get total replicated and single nodes
get_n_unique_genetic_hits <- function(g, curr_morbidity){
  if (nrow(g$edges)>0){
    ## get number of nodes in each cluster
    replicated <- sum(!is.na(unique(g$cn$cluster)))
    # how many nodes aren't in a cluster i.e. don't overlap with any other nodes (i.e. unreplicated) 
    singles <- g$cn %>% filter(is.na(cluster)) %>% nrow()
  } else {
    replicated <- 0
    singles <- nrow(g$nodes)
  }
  n_genetic_hits <- data.frame(morbidity = curr_morbidity, replicated = replicated, single = singles, total = replicated + singles)
  return(n_genetic_hits)
}

## For each disease, get the graph and count clusters
get_diseases_individually <- function(ard_leads, overlaps, overwrite = FALSE, save_dir = default_save_dir){
  individual_graphs <- list()
  individual_communities <- data.frame()
  if (overwrite){
    d <- ard_leads$morbidity %>% unique()
    for (i in seq_along(d)){
      individual_graphs[[i]] <- count_communities(d[i], ard_leads, overlaps$coloc_within, overlaps$overlap_within, min_jaccard = 0, plot = FALSE, detect_subgraph_communities = FALSE)
      n_unique_genetic_hits <- get_n_unique_genetic_hits(individual_graphs[[i]], d[i])
      individual_communities <- individual_communities %>% bind_rows(n_unique_genetic_hits)
    }
    save(individual_graphs, file = paste0(save_dir, "individual_disease_graphs.Rda"))
    save(individual_communities, file = paste0(save_dir, "individual_disease_n_communities.Rda"))
  } else {
    load(paste0(save_dir, "individual_disease_graphs.Rda"))
    load(paste0(save_dir, "individual_disease_n_communities.Rda"))
  }
  return(list("individual_graphs" = individual_graphs, 
              "individual_communities" = individual_communities))
}

## Calculate overlap between each pair of diseases in a graph
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
  n_pairwise_shared_clusters <- cbind(data.frame(pairwise_combinations, count)) 
  
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
  n_tests <- nrow(pairwise_combinations)
  
  ## Correct for multiple testing  
  ft_all <- pairwise_shared %>% 
    rowwise() %>% 
    mutate(p.value = get_fishers_exact(overlap = count, morb_x_overall = X_overall, morb_y_overall = Y_overall, total = total)) %>%
    ungroup() %>% 
    mutate(p.adjust = p.adjust(p.value, method = "BH", n = n_tests))   # n tests
  # mutate(p.adjust = p.adjust(p.value, method = "BH", n = nrow(.)))   # n possible pairs
  
  return(ft_all)
}

## Enrich hallmarks of ageing
enrich_hallmarks <- function(gene_sets){
  hallmark_term_genes <- read.csv("data/full_goterm_genes.csv") %>% 
    filter(!is.na(EVIDENCE)) %>% 
    dplyr::select(goId = GO, gene = SYMBOL) %>% 
    unique() %>%
    group_by(goId) %>% mutate(n_reference = length(goId))
  
  hallmark_genes <- read.csv("data/full_goterm_list.csv") %>% 
    left_join(hallmark_term_genes) %>%
    dplyr::select(goHallmark, gene) %>%                                   # overall hallmark
    # select( goHallmark, goHallmarkId, goHallmarkTerm, gene) %>%    # high terms
    #dplyr::select( goHallmark, goHallmarkId, goHallmarkTerm, goTerm, goId, gene) %>%    # all terms
    unique() %>%
    group_by(goHallmark) %>% 
    #  group_by(goHallmarkId, goHallmark) %>% 
    #group_by(goHallmark, goHallmarkId, goHallmarkTerm, goId, goTerm) %>% 
    mutate(total_reference = length(gene)) %>%
    filter(total_reference >= 12)
  
  sets = bind_rows(data.frame(set = "TargetAge", gene = gene_sets$TargetAge),
                   data.frame(set = "CellAge", gene = gene_sets$CellAge),
                   data.frame(set = "GenAge", gene = gene_sets$GenAge)) %>%
    group_by(set) %>%  
    mutate(total_query = length(set))
  
  enriched <- inner_join(sets, hallmark_genes) %>% 
     group_by(set, total_query, goHallmark, total_reference) %>% 
    #  group_by(set, total_query, goHallmark, goHallmarkId,  goHallmarkTerm, total_reference) %>% 
    #group_by(set, total_query, goHallmark, goHallmarkId,  goHallmarkTerm, goId, goTerm, total_reference) %>% 
    summarise(overlap = length(gene), genes = paste(gene, collapse = ", ")) %>% 
    mutate(background = length(unique(hallmark_genes$gene))) %>% 
    rowwise() %>% 
    mutate(p.value = get_fishers_exact(overlap = overlap, morb_x_overall = total_query, morb_y_overall = total_reference, total = 19938)) %>%
    mutate(p.adjust = p.adjust(p.value, method = "BH", n = length(unique(hallmark_genes$goHallmark))))   # n Hallmark sets
    #  mutate(p.adjust = p.adjust(p.value, method = "BH", n = length(unique(hallmark_genes$goHallmark))))   # n Hallmark sets
    #mutate(p.adjust = p.adjust(p.value, method = "BH", n = length(unique(hallmark_genes$goTerm))))   # n Hallmark sets
  
  return(enriched)
}

## Enrich hallmarks of ageing
enrich_small_hallmarks <- function(gene_sets){
  hallmark_term_genes <- read.csv("data/full_goterm_genes.csv") %>% 
    filter(!is.na(EVIDENCE)) %>% 
    dplyr::select(goId = GO, gene = SYMBOL) %>% 
    unique() %>%
    group_by(goId) %>% mutate(n_reference = length(goId))
  
  hallmark_genes <- read.csv("data/full_goterm_list.csv") %>% 
    left_join(hallmark_term_genes) %>%
    #dplyr::select(goHallmark, gene) %>%                                   # overall hallmark
     dplyr::select( goHallmark, goHallmarkId, goHallmarkTerm, gene) %>%    # high terms
    #dplyr::select( goHallmark, goHallmarkId, goHallmarkTerm, goTerm, goId, gene) %>%    # all terms
    unique() %>%
    #group_by(goHallmark) %>% 
      group_by(goHallmarkId, goHallmark) %>% 
    #group_by(goHallmark, goHallmarkId, goHallmarkTerm, goId, goTerm) %>% 
    mutate(total_reference = length(gene)) %>%
    filter(total_reference >= 12)
  
  sets = bind_rows(data.frame(set = "TargetAge", gene = gene_sets$TargetAge),
                   data.frame(set = "CellAge", gene = gene_sets$CellAge),
                   data.frame(set = "GenAge", gene = gene_sets$GenAge)) %>%
    group_by(set) %>%  
    mutate(total_query = length(set))
  
  enriched <- inner_join(sets, hallmark_genes) %>% 
    #group_by(set, total_query, goHallmark, total_reference) %>% 
      group_by(set, total_query, goHallmark, goHallmarkId,  goHallmarkTerm, total_reference) %>% 
    #group_by(set, total_query, goHallmark, goHallmarkId,  goHallmarkTerm, goId, goTerm, total_reference) %>% 
    summarise(overlap = length(gene), genes = paste(gene, collapse = ", ")) %>% 
    #mutate(background = length(unique(hallmark_genes$gene))) %>% 
    mutate(background = 19938) %>% 
    rowwise() %>% 
    mutate(p.value = get_fishers_exact(overlap = overlap, morb_x_overall = total_query, morb_y_overall = total_reference, total = background)) %>%
  #  mutate(p.adjust = p.adjust(p.value, method = "BH", n = length(unique(hallmark_genes$goHallmark))))   # n Hallmark sets
    mutate(p.adjust = p.adjust(p.value, method = "BH", n = length(unique(hallmark_genes$goHallmark))))   # n Hallmark sets
  #mutate(p.adjust = p.adjust(p.value, method = "BH", n = length(unique(hallmark_genes$goTerm))))   # n Hallmark sets
  
  return(enriched)
}

###########################################################################
###########################################################################
###                                                                     ###
###                       SECTION 3:                                    ###
###                     Query OT APIs                                   ###
###                                                                     ###
###########################################################################
###########################################################################
## Initialise OT Genetics Portal API
initialise_api <- function(){
  library(ghql)
  library(jsonlite)
  cli <- GraphqlClient$new(
    url = "https://api.genetics.opentargets.org/graphql"
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
      bind_rows() %>%
      # remove low confidence (L2G<0.05)
      filter(yProbaModel > 0.05)
    
    ## Do all the other chunks
    for (i in seq(1, length(l2g_all))){
      print(i)
      l2g_all[[i]] <- studies_variants_split[[i]] %>% 
        group_by(studyId, lead_variantId) %>% 
        group_split() %>% 
        ## API call for each studyID + variantID
        purrr::map(~fetch_l2g(df = ., variables = list(studyId = .$studyId, variantId = .$lead_variantId))) %>%
        bind_rows() %>%
        filter(yProbaModel > 0.05)
      Sys.sleep(3)
    }
    
    ## Collapse into one dataframe
    l2g_all_joined <- l2g_all %>% bind_rows() %>% jsonlite::flatten() 
    
    ## Use more up-to-date gene symbol from the Platform 
    l2g_all_joined <- l2g_all_joined %>% left_join(target_annotations %>% select(gene.id = targetId, targetSymbol))
    
    ## Save output
    save(l2g_all_joined, file=paste0(save_dir, "ltg_all.Rda"))
    
  } else {
    load(paste0(save_dir, "ltg_all.Rda"))
  }
  
  return(l2g_all_joined)
  
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

###########################################################################
###########################################################################
###                                                                     ###
###                       SECTION 4:                                    ###
###                      Enrichment                                     ###
###                                                                     ###
###########################################################################
###########################################################################



perform_enrichment <- function(gene_sets_entrez, overwrite = default_overwrite, save_dir = default_save_dir){
  if (overwrite){
    library(clusterProfiler)
    library(ReactomePA)
    library(dplyr) 
    library(patchwork)
    library(ggplot2)
    
    cluster_compa_reactome <- compareCluster(geneCluster = gene_sets_entrez, fun = "enrichPathway")
    cluster_compa_go_BP <- compareCluster(geneCluster = gene_sets_entrez, fun = "enrichGO", OrgDb = "org.Hs.eg.db", ont = "BP")
    cluster_compa_go_MF <- compareCluster(geneCluster = gene_sets_entrez, fun = "enrichGO", OrgDb = "org.Hs.eg.db", ont = "MF")
    
    reactome_df <- as.data.frame(cluster_compa_reactome) %>% mutate(set = "Reactome Pathways")
    go_bp_df <- as.data.frame(cluster_compa_go_BP) %>% mutate(set = "GO Biological Processes")
    go_mf_df <- as.data.frame(cluster_compa_go_MF) %>% mutate(set = "GO Molecular Functions")
    
    all_enriched <- bind_rows(go_bp_df, go_mf_df, reactome_df)
    
    save(all_enriched, file = paste0(save_dir, "enrichment_output.Rda"))
  } else {
    load(paste0(save_dir, "enrichment_output.Rda"))
  }
  return(all_enriched)
}


###########################################################################
###########################################################################
###                                                                     ###
###                       SECTION 4:                                    ###
###                 Summaries/figures for paper                         ###
###                                                                     ###
###########################################################################
###########################################################################

## Look at ancestry of GWAS studies
get_ancestries <- function(ard_leads) {
  ancestries <-
    ard_leads %>% 
    select(studyId, ancestry_initial) %>% 
    unique() %>% 
    unnest(ancestry_initial) %>% 
    separate(ancestry_initial,
             into = c("ancestry", "n"),
             sep = "=") %>% 
    group_by(studyId) %>% 
    mutate(n = as.numeric(n),
           total = sum(n),
           proportion = n / total)
  majority_population <-
    ancestries %>% 
    group_by(studyId) %>% 
    top_n(1, proportion) %>%
    mutate(single_population = ifelse(proportion == 1, TRUE, FALSE))
  return(majority_population)
}

## Get the details for a table on number, size, etc of GWAS
get_genetics_table <- function(ard_leads, diseases, g_all, individual_diseases){
  ard_studies <- ard_leads %>% 
    select(studyId, 
           n_cases, 
           n_initial) %>%
    unique() %>%
    mutate(diseaseName = "Total", morbidity = "Total")
  
  study_size <- ard_leads %>% 
    select(morbidity, 
           diseaseName, 
           studyId, 
           n_cases, 
           n_initial) %>% 
    group_by(diseaseName) %>% 
    bind_rows(ard_studies) %>% 
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
    mutate(diseaseName = "Total", morbidity="Total")
  
  gwas_tbl <- 
    ard_leads %>% 
    group_by(diseaseId, 
             diseaseName,
             morbidity, 
             specificDiseaseId, 
             specificDiseaseName, 
             studyId)  %>% 
    summarise(has_sumstats = max(has_sumstats), 
              n_variants = length(studyId)) %>% 
    ungroup() %>% 
    group_by(diseaseId, diseaseName, morbidity, specificDiseaseId, specificDiseaseName) %>% 
    summarise(n_gwas_studies = length(specificDiseaseName), 
              has_sumstats = sum(has_sumstats), 
              n_variants = sum(n_variants))
  
  summarised_gwas_tbl <- 
    gwas_tbl %>% 
    ungroup() %>% 
    group_by(diseaseName, morbidity,) %>% 
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
    full_join(bind_rows(individual_diseases$individual_communities,
                        get_n_unique_genetic_hits(g_all, "Total"))) %>%
    select(-morbidity) %>%
    select(diseaseName, n_gwas_studies, has_sumstats, min, max, median, n_variants, everything())
  
  full_gwas_tbl <- gwas_tbl %>% 
    full_join(diseases) %>% 
    ungroup() %>% 
    select(-therapeuticAreas, -morbidity) %>%
    filter(!is.na(n_gwas_studies) | diseaseId == specificDiseaseId) %>%
    arrange(diseaseName) 
  return(list("summarised" = summarised_gwas_tbl, "full" = full_gwas_tbl))
}

## Upload to google drive
update_gs_genetics_table <- function(genetics_tables){
  library(googlesheets4)
  gs4_create("SITable1", sheets = genetics_tables$summarised)
  gs4_create("SITable2", sheets = genetics_tables$full)
}

## QQ plot
get_qq_plot <- function(pairwise_overlaps, overwrite = default_overwrite, save_dir = default_save_dir){
  x <- 10^(-seq(0,8,0.1))
  y <- runif(100000)
  df <- data.frame(x = -log10(x), y = -quantile(log10(pairwise_overlaps$p.value), probs = x), control = -quantile(log10(y), probs = x))
  gg <- ggplot(data = df) + 
    geom_point(aes(x = x, y = y)) + 
    geom_abline(slope = 1) + 
    labs(x = "Expected quantiles", y = "Observed quantiles", title = "QQ plot of ARD overlap significance") + 
    theme_bw() +
    theme(panel.grid = element_blank())
  if (overwrite){
    ggsave(gg, file = paste0(save_dir, "QQ.png"), width = 4.5, height = 4.5, dpi = 300)
    ggsave(gg, file = paste0(save_dir,"QQ.pdf"), width = 4.5, height = 4.5)
  }
  gg
  
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


calculate_pairwise_overlap <- function(g){
  unreplicated_associations <- g_all$cn %>% filter(is.na(cluster)) %>% count(morbidity)
}


## make the heatmap figure showing overlaps between phenotypes
plot_shared_heatmap <- function(pairwise_overlaps, detailed = TRUE, overwrite = default_overwrite, save_dir = default_save_dir){
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
    ## bold around statistically significant overlap 
    #  geom_tile(data = subset(full_df, p.value < 0.05), colour = textcol, size = 0.65) + 
    #   geom_text(aes(label = count, colour = proportion_X > 0.45, fontface = face), size = 2.1) + 
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
  
  ## label significant overlaps
  if (detailed){
    shared_heatmap = shared_heatmap + 
      geom_text(aes(label = ifelse(p.adjust < 0.05, paste0(count, "*"), count), colour = proportion_X > 0.45), size = 2.1)
    ext = "SI"
    fig_height = 7.5
    fig_width = 8
  } else {
    shared_heatmap = shared_heatmap + 
      geom_text(aes(label = ifelse(p.adjust < 0.05, "*", NA), colour = proportion_X > 0.45), size = 2.1)
    ext = ""
    fig_height = 6.5
    fig_width = 7
  }
  
  
  if (overwrite) {
    ggsave(shared_heatmap, file = paste0(save_dir, ext, "Fig1.png"), width = fig_width, height = fig_height, dpi = 300)
    ggsave(shared_heatmap, file = paste0(save_dir, ext, "Fig1.pdf"), width = fig_width, height = fig_height)
    ggsave(shared_heatmap, file = paste0(save_dir, ext, "Fig1.svg"), width = fig_width, height = fig_height)
    
  }
  
  return(shared_heatmap)
}

## Test significance of gene set overlaps 
test_overlap <- function(gene_sets, set1 = "TargetAge", set2){
  in_set1 = length(gene_sets[set1][[1]])
  in_set2 = length(gene_sets[set2][[1]])
  n_OT = 19938   ## Version 21.06 (protein-coding)
  in_both = sum(gene_sets[set1][[1]] %in% gene_sets[set2][[1]])
  
  dat <- data.frame(
    "in_set1" = c(in_both, in_set1 - in_both),
    "not_in_set1" = c(in_set2 - in_both, n_OT - in_set1 - in_set2 + in_both),
    row.names = c("in_set2", "not_in_set2"),
    stringsAsFactors = FALSE
  )
  pval = fisher.test(dat)$p.value
  return(list(contingency = dat, pval = pval))
}

## plot hallmarks barcharts
plot_hallmarks_barplot <- function(enriched_hallmarks, save_dir = default_save_dir, overwrite = default_overwrite){
  hallmarks_barplot <- enriched_hallmarks %>% 
    mutate(set = factor(set, levels = c("TargetAge", "GenAge", "CellAge"))) %>% 
    ggplot(., aes(fill = set, x=overlap, y = reorder(goHallmark, -overlap))) + 
    geom_col(position = "dodge") + 
    theme_classic() + 
    labs(x = "Number of genes", subtitle = "Hallmarks of Ageing Genes") + 
    theme(panel.grid = element_blank(), 
          axis.title.y = element_blank(), 
          #     legend.position = "bottom", 
          legend.position = c(0.8, 0.85),
          legend.title = element_blank(),
          legend.text = element_text(margin = margin(r = 10, unit = "pt"))) + 
    scale_x_continuous(expand = c(0,0))
  if (overwrite){
    ggsave(hallmarks_barplot, width = 5, height = 4.1, dpi = 600, file = paste0(save_dir, "hallmarks_bar.png"))
    ggsave(hallmarks_barplot, width = 5, height = 4.1, file = paste0(save_dir, "hallmarks_bar.pdf"))
  }
  return(hallmarks_barplot)
}

plot_mini_hallmarks_barplot<- function(enriched_hallmarks, save_dir = default_save_dir, overwrite = default_overwrite){
  hallmarks_barplot_ta <- 
    enriched_hallmarks %>% 
    mutate(perc = overlap/total_query) %>% 
    filter(set=="TargetAge") %>% 
    ggplot(., aes(fill = set, x=overlap, y = reorder(goHallmark, -overlap))) + 
    geom_col(position = "dodge") + 
    theme_classic() + 
    theme(panel.grid = element_blank(), 
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "none", 
          legend.title = element_blank()) + 
    scale_x_continuous(expand = c(0,0))
  if (overwrite){
    ggsave(hallmarks_barplot_ta, width = 3, height = 2, dpi = 600, file = paste0(save_dir, "hallmarks_mini.png"))
    ggsave(hallmarks_barplot_ta, width = 3, height = 2, file = paste0(save_dir, "hallmarks_mini.pdf"))
  }
  return(hallmarks_barplot_ta)
}
