library(clusterProfiler)
library(ReactomePA)
library(dplyr) 
library(patchwork)
library(ggplot2)

load("../gene_sets.Rda")
load("../gene_sets_entrez.Rda")

genes <- read.table(file = "../top_L2G_genes.txt", col.names = c("ENSEMBL", "symbol", "groups"))

### In order to perform the pathway comparison, a list must be generated and the gene ids should be entrez_id

# 1.- Convert the gene_symbols to entrez_ids by using clusterprofiler
gene_sets_list <- list(
  TargetAge = bitr(gene_sets_entrez$TargetAge, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID,
  GenAgeHuman = as.character(gene_sets_entrez$GenAgeHuman),
  CellAge = as.character(gene_sets_entrez$CellAge)
)

# 2.- Perform the analysis
cluster_compa_reactome <- compareCluster(geneCluster = gene_sets_list, fun = "enrichPathway")
cluster_compa_go_BP <- compareCluster(geneCluster = gene_sets_list, fun = "enrichGO", OrgDb = "org.Hs.eg.db", ont = "BP")
cluster_compa_go_MF <- compareCluster(geneCluster = gene_sets_list, fun = "enrichGO", OrgDb = "org.Hs.eg.db", ont = "MF")


enrichment_theme <-  function(gp) {
  gp + 
    #  scale_color_gradient(limits = c(0, 0.02), breaks = seq(0, 0.02, 0.01), low = "#ca0020", high = "#0571b0") + 
    #    scale_size_continuous(range = c(0, 0.25), breaks = seq(0, 0.25, 0.05), limits = c(0, 0.25)) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}


p_go_BP <- dotplot(cluster_compa_go_BP, font.size=10) %>% enrichment_theme() + scale_color_gradient(limits = c(0, 0.02), breaks = seq(0, 0.02, 0.01), low = "#ca0020", high = "#0571b0") +  labs(subtitle = "GO Biological Processes")
p_go_MF <- dotplot(cluster_compa_go_MF, font.size=10)  %>% enrichment_theme() + scale_color_gradient(limits = c(0, 0.02), breaks = seq(0, 0.02, 0.01), low = "#ca0020", high = "#0571b0") + guides(size = FALSE) + labs(subtitle = "GO Molecular Function")
p_reactome <- dotplot(cluster_compa_reactome, font.size=10) %>% enrichment_theme() + scale_color_gradient(limits = c(0, 0.02), breaks = seq(0, 0.02, 0.01), low = "#ca0020", high = "#0571b0") + guides(size = FALSE) + labs(subtitle = "Reactome Pathway")
p_all <- p_go_BP / p_go_MF / p_reactome / plot_layout(guides = "collect")
save(p_go_BP, file = "enrichment_plot_1.Rda")
save(p_go_MF, file = "enrichment_plot_2.Rda")
save(p_reactome, file = "enrichment_plot_3.Rda")



go_enrichment_simplified <- simplify(cluster_compa_go)
dotplot(go_enrichment_simplified)
dotplot(go_enrichment_simplified, 25)


# 3.- Visualizing the results (open the plot in a separate window to see properly)
p <- dotplot(cluster_compa_reactome) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dotplot(cluster_compa_go, showCategory = 25)

result <- cluster_compa_reactome@compareClusterResult
result %>% group_by(Cluster) %>% top_n(10, -p.adjust) %>% ggplot(., aes(y = Description, x = Count, fill = p.adjust)) + geom_col() + facet_wrap(~Cluster, ncol = 3) + theme_minimal()

## 5.- Accessing Reactome IDs and names
# results are in:
cluster_compa@compareClusterResult

# E.g.:
# Reactome term IDs and descriptions
cluster_compa_reactome@compareClusterResult$ID
cluster_compa_reactome@compareClusterResult$Description
# Gene IDs belonging to each term
cluster_compa_reactome@compareClusterResult$geneID

