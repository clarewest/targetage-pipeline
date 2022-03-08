textcol <- "grey40"

prop_shared <- individual_communities %>% 
  select(morbidity_A = morbidity, total_A = total) %>% 
  full_join(overlaps_df) %>% 
  mutate(proportion_A = shared_A/total_A) %>% 
  select(morbidity_A, morbidity_B, total_A, shared_A, proportion_A)

### making a heatmap
coloc_pre_heatmap <- prop_shared %>% 
  select(-shared_A, -total_A) %>% 
  pivot_wider(names_from = morbidity_A, values_from = proportion_A) %>%  
  replace(is.na(.), 0)

library("ggdendro")

# Run clustering
coloc_mat <- as.matrix(coloc_pre_heatmap[,-1])
rownames(coloc_mat) <- coloc_pre_heatmap$morbidity.y
coloc.dendro <- as.dendrogram(hclust(d = dist(x = coloc_mat)))

# Create dendro
dendro.plot <- ggdendrogram(data = coloc.dendro, rotate = TRUE)

# Preview the plot
print(dendro.plot)

## reorder to the clustered order
coloc.order <- order.dendrogram(coloc.dendro)

coloc_pre_heatmap_long <- coloc_pre_heatmap %>% 
  replace(. == 0, NA) %>%
  pivot_longer(-morbidity_B, names_to = "morbidity_A", values_to = "proportion") %>%
  left_join(prop_shared %>% 
              select(-proportion_A), 
            by = c("morbidity_A", "morbidity_B")) 

coloc_pre_heatmap_long$morbidity_B <- factor(x = coloc_pre_heatmap_long$morbidity_B,
                                             levels = coloc_pre_heatmap$morbidity_B[coloc.order],
                                             ordered = TRUE)
coloc_pre_heatmap_long$morbidity_A <- factor(x = coloc_pre_heatmap_long$morbidity_A,
                                             levels = coloc_pre_heatmap$morbidity_B[coloc.order],
                                             ordered = TRUE)

coloc_heatmap <- coloc_pre_heatmap_long %>%
  # mutate_at((c("morbidity_A", "morbidity_B")), str_replace_all, pattern ="chronic obstructive pulmonary disease", replacement = "COPD") %>%
  ggplot(., aes(x=morbidity_A, y = morbidity_B, fill = proportion)) +
  geom_tile(colour = textcol) +
  geom_text(aes(label = shared_A, colour = proportion > 0.45), size = 2.1) + 
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


ggsave(coloc_heatmap, file = "coloc_heatmap.png", width = 8, height = 7.5, dpi = 300)
