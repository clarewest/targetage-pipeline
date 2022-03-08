rerun = FALSE

if(rerun){
  combs <- combn(d, 2, simplify = FALSE)
  pairwise_graphs <- list()
  for (i in seq_along(combs)){
    pairwise_graphs[[i]] <- count_communities(combs[[i]], ard_leads, coloc_within, overlap_within, plot = FALSE)
  }
  overlaps <- pairwise_graphs %>% purrr::map("overlap")
  overlaps_df <- pairwise_graphs %>% purrr::map(8) %>% dplyr::bind_rows()
  
  save(overlaps_df, file = "genetic_overlaps_df.Rda")
  save(pairwise_graphs, file = "genetic_overlaps_graphs.Rda")
} else {
  load("genetic_overlaps_df.Rda")
  load("genetic_overlaps_graphs.Rda")
}