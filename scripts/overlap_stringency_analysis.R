library(tidyverse)

source("graph_functions.R")

save_files = FALSE
save_dir = "data/analysis/"

load(paste0(save_dir, "ard_leads_filtered.Rda"))
load(paste0(save_dir, "colocalisation_within.Rda"))
load(paste0(save_dir, "overlap_within.Rda"))
d <- ard_leads$morbidity %>% unique()

jaccard_communities_results <- data.frame()

for (curr_min_jaccard in seq(0, 1, 0.1)){
  print(curr_min_jaccard)
  g_jaccard <- count_communities(d, ard_leads, coloc_within, overlap_within, min_jaccard = curr_min_jaccard, detect_subgraph_communities = TRUE, plot = FALSE)
  jaccard_communities_results <- bind_rows(g_jaccard$cn %>% mutate(min_jaccard = curr_min_jaccard), jaccard_communities_results)
}

jaccard_communities_results_coloc <- data.frame()
  
  for (curr_min_jaccard in seq(0, 1, 0.1)){
    print(curr_min_jaccard)
    g_jaccard <- count_communities(d, ard_leads, coloc_within, overlap_within, min_jaccard = curr_min_jaccard, detect_subgraph_communities = TRUE, plot = FALSE)
    jaccard_communities_results_coloc <- bind_rows(g_jaccard$cn %>% mutate(min_jaccard = curr_min_jaccard), jaccard_communities_results_coloc)
  }

if (save_files){
  save(jaccard_communities_results, file = "jaccard_communities.Rda")
  save(jaccard_communities_results_coloc, file = "jaccard_communities_coloc.Rda")
}

jaccard_communities_results <- jaccard_communities_results %>% mutate(edges = "finemapping") 
jaccard_communities_results_coloc <- jaccard_communities_results_coloc %>% mutate(edges = "coloc only")
jaccard_communities_results <- bind_rows(jaccard_communities_results, jaccard_communities_results_coloc)

ji_n_clusters_comp <-  jaccard_communities_results_all  %>% mutate(type = "all") %>% bind_rows(jaccard_communities_results %>% filter(n_morbidities > 1) %>% mutate(type = ">1 morbidity"))  %>% group_by(edges, type, min_jaccard) %>% summarise(n_clusters = length(unique(cluster))) %>% ggplot(., aes(x=min_jaccard, y = n_clusters, colour = edges)) + geom_line() + geom_point() + facet_wrap(~type, labeller = as_labeller(c(`>1 morbidity` = "Clusters with >1 morbidity", `all` = "All clusters"))) + theme_bw() + labs(x = "Jaccard index threshold", y = "Number of disconnected components (clusters)", subtitle = "Number of communities detected with varying overlap stringency")

ji_n_clusters <-  jaccard_communities_results_coloc %>% mutate(type = "all") %>% bind_rows(jaccard_communities_results_coloc %>% filter(n_morbidities > 1) %>% mutate(type = ">1 morbidity"))  %>% group_by(type, min_jaccard) %>% summarise(n_clusters = length(unique(cluster))) %>% ggplot(., aes(x=min_jaccard, y = n_clusters)) + geom_line() + geom_point() + facet_wrap(~type, labeller = as_labeller(c(`>1 morbidity` = "Clusters with >1 morbidity", `all` = "All clusters"))) + theme_bw() + labs(x = "Jaccard index threshold", y = "Number of disconnected components (clusters)", subtitle = "Number of clusters detected with varying overlap stringency") + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank()) + scale_x_continuous(breaks = seq(0, 1, 0.2)) + scale_y_continuous(limits = c(0, 3500), expand = c(0, NA)) 

ji_n_communities <-  jaccard_communities_results_coloc %>% 
  mutate(type = "all") %>% 
  bind_rows(jaccard_communities_results_coloc %>% 
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


ji_si_figure <- ji_n_clusters / ji_n_communities + plot_annotation(tag_levels = 'a')

ggsave(ji_si_figure, file = "si_jaccard_index.png", width = 8, height = 7.5, dpi = 300)
ggsave(ji_si_figure, file = "si_jaccard_index.pdf", width = 8, height = 7.5)

chosen <- jaccard_communities_results %>% filter(edges == "coloc only", min_jaccard == 0)

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