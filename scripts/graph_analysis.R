source("graph_functions.R")

if (rerun){
  g_all <- count_communities(d, ard_leads, coloc_within, overlap_within, min_jaccard = 0, detect_subgraph_communities = TRUE, plot = FALSE)
  save(g_all, file = paste0(save_dir, "graph_all_morbidities.Rda"))
} else {
  load()
}
