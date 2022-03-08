source("graph_functions")

if (rerun) {
  for (i in seq_along(d)) {
    individual_graphs[[i]] <-
      count_communities(d[i], ard_leads, coloc_within, overlap_within, plot = FALSE)
    if (nrow(individual_graphs[[i]]$edges) > 0) {
      ## get number of nodes in each cluster
      replicated <-
        sum(!is.na(unique(individual_graphs[[i]]$cn$cluster)))
      # how many nodes aren't in a cluster i.e. don't overlap with any other nodes (i.e. unreplicated)
      singles <-
        individual_graphs[[i]]$cn %>% filter(is.na(cluster)) %>% nrow()
      ## for clusters, how many consist of just one morbidity (distinct) and how many contain both (shared)
    } else {
      replicated <- 0
      singles <- nrow(individual_graphs[[i]]$nodes)
    }
    individual_communities <-
      individual_communities %>% bind_rows(
        data.frame(
          morbidity = d[i],
          replicated = replicated,
          single = singles,
          total = replicated + singles
        )
      )
  }
  save(individual_graphs, file = "individual_disease_graphs.Rda")
  save(individual_communities, file = "individual_disease_n_communities.Rda")
} else {
  load("individual_disease_graphs.Rda")
  load("individual_disease_n_communities.Rda")
}