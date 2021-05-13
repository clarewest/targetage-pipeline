#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
library(shiny)
library(shinyjs)
library(shinyWidgets)
library(shinydashboard)
library(tidyverse)
library(DT)
library(reactable)
library(ggtext)
library(ggiraph)
library(visNetwork)
library(igraph)

#### Data ####
load("ltg_all.Rda")
load("graph_all_morbidities.Rda")
load("ard_leads_filtered.Rda")

ard_leads <- ard_leads %>% 
  mutate(specificDiseaseName = recode(specificDiseaseName, 
                                      `low density lipoprotein cholesterol measurement` = "LDL cholesterol measurement",
                                      `high density lipoprotein cholesterol measurement` = "HDL cholesterol measurement",
                                      `very low density lipoprotein cholesterol measurement` = "VLDL cholesterol measurement",
                                      `type II diabetes mellitus` = "type 2 diabetes")) %>%
  mutate(specificDiseaseName = ifelse(specificDiseaseName != morbidity, specificDiseaseName, NA))
top_l2g <- l2g_all_joined %>%
    group_by(studyId, lead_variantId) %>% 
    top_n(1,yProbaModel ) %>% ## need to remove duplicates for 5 variants 
    select(studyId, lead_variantId, gene.symbol, gene.id, L2G = yProbaModel, hasColoc, distanceToLocus)

top_genes <- top_l2g %>%
    mutate(gene.symbol = ifelse(hasColoc == TRUE, paste0('atop(bold("',gene.symbol,'")'), gene.symbol)) %>%
    top_n(1, L2G) %>%
    summarise(gene.symbol = paste0(gene.symbol, collapse = ", "),
              L2G = max(L2G))

tbl <- g_all$cn %>%
    select(id, cluster, n_morbidities, morbidity, has_sumstats, everything(), -color, -label, -c_size) %>% 
    filter(!is.na(cluster)) %>% 
    filter(n_morbidities > 1) %>% 
    arrange(-n_morbidities, cluster) %>%
    left_join(top_l2g) %>%
    left_join(ard_leads %>% select(studyId, lead_variantId, morbidity, specificDiseaseName))

groups <- bind_rows(g_all$coloc_edges, g_all$edges)  %>% rownames_to_column("group") %>%  pivot_longer(c("from", "to"), values_to = "id")


tbl_leads <- tbl %>%
    left_join(ard_leads %>% select(-c("morbidity", "diseaseId", "diseaseName", "specificDiseaseId",   "specificDiseaseName")) %>% unique())

tbl_genes <- tbl_leads  %>% 
    left_join(top_l2g) 

## Custom JS aggregate function for reactable that drops NA
uniqueDropNA <- JS("function(values) {
      const isNotNA = function(x) { return x != null && x !== 'NA' }
      let uniqueValues = [...new Set(values)]
      return uniqueValues.filter(isNotNA).join(', ')  
    }")

#### A function to view visNetwork graph of genetic association cluster
view_cluster <- function(nodes, g, curr_cluster){
    cluster_nodes <- nodes %>% filter(cluster %in% curr_cluster) %>%
      mutate(color = ifelse(color == "black", "#5E4B56", color)) %>% 
  #    mutate(color.background = ifelse(has_sumstats == TRUE, color, paste0(color, "90"))) %>%   ### make transparent if no sumstats
      mutate(color.border = ifelse(has_sumstats == TRUE, "black", color),
             borderWidth = 2) %>%
      rename(color.background = color) 
      
    cluster_edges <- bind_rows(g$edges %>% mutate(width = 1), g$coloc_edges %>% mutate(width = 3)) %>% 
        filter(from %in% cluster_nodes$id) %>%
        filter(to %in% cluster_nodes$id) 
      
    vn <- 
        visNetwork(
            cluster_nodes,
            cluster_edges,
            height = "1000px",
            width = "100%"
        ) %>% 
        #    visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
        visIgraphLayout()
    return(vn)
}

## Function to plot OR/beta for each variant in a cluster
plot_forest <- function(leads, curr_cluster){
   # labels <- leads %>% filter(cluster == curr_cluster) %>% select(cluster, id, lead_variantId, gene.symbol, lead_variantId, L2G) %>% mutate(labels = ifelse(is.na(gene.symbol), lead_variantId, paste0(gene.symbol, " (", round(L2G, 2), ")    ", lead_variantId))) 
    cols <- bind_rows(g_all$colours, c(color = "#5E4B56", morbidity = "combination"))
    labels <- leads %>% 
        filter(cluster == curr_cluster) %>% 
        select(cluster, id, lead_variantId, gene.symbol, lead_variantId, L2G, hasColoc) %>% 
        mutate(bestGene = ifelse(is.na(gene.symbol), lead_variantId, paste0(gene.symbol, " (", round(L2G, 2), ")"))) %>%
     #   mutate(labels = ifelse(is.na(gene.symbol), lead_variantId, paste0(gene.symbol, " (", round(L2G, 2), ")    ", lead_variantId)))  %>%
        mutate(name = ifelse(hasColoc == TRUE, glue::glue("<b>{bestGene}</b>   ({lead_variantId})"), glue::glue("{bestGene} ({lead_variantId})")))
    harmonise_oddsr <- function(x) (ifelse(x<1, 1/x, x))
    harmonise_beta <- function(x) (abs(x))
    
    plot_data <- 
        tbl_genes %>% 
        filter(cluster == curr_cluster) %>% 
        select(id, cluster, lead_variantId, beta, odds_ratio, contains("_ci_"), morbidity, specificDiseaseName, trait_reported, gene.symbol, has_sumstats, direction) %>% 
        rename(oddsr = odds_ratio) %>% 
        unique() %>%  ### duplicates arise where a single gwas is mapped to two EFO codes e.g. "OA of hip or knee"
        mutate_at(vars(contains("oddsr")), harmonise_oddsr) %>%
        mutate_at(vars(contains("beta")), harmonise_beta) %>%
        pivot_longer(contains("_ci_"), names_to = c("measure", "ci"), names_sep = "_ci_")  %>% 
        mutate(score = ifelse(measure=="beta", beta, oddsr)) %>% 
        pivot_wider(names_from = "ci", values_from = "value") %>%
        filter(!is.na(score)) %>%
     #   droplevels() %>%
        left_join(labels) %>%
        mutate(morbidity = ifelse(str_detect(morbidity, "\\+"), "combination", morbidity)) %>%
        mutate(morbidity = factor(morbidity, levels = cols$morbidity))
    origin = data.frame(measure = c("beta", "oddsr"), intercept = c(0, 1)) %>%
        filter(measure %in% plot_data$measure)
    
    main_plot <- plot_data %>% 
        ggplot(., aes(y = id, x=score, xmin = lower, xmax = upper, shape = direction)) + 
        geom_vline(data = origin, aes(xintercept = intercept)) + 
      geom_errorbarh_interactive(aes(colour = morbidity, tooltip = trait_reported, data_id = lead_variantId), height = 0.1) + 
        geom_point_interactive(aes(colour = morbidity, fill = morbidity, tooltip = trait_reported, data_id = lead_variantId)) + 
      geom_point_interactive(data = subset(plot_data, has_sumstats == TRUE), aes(tooltip = trait_reported, data_id = lead_variantId), colour = "black") + 
        theme_bw() + 
        scale_y_discrete(breaks = labels$id, labels = labels$name) + 
        scale_colour_manual(values = cols$color, breaks = cols$morbidity) + 
        scale_fill_manual(values = cols$color, breaks = cols$morbidity) + 
        scale_shape_manual(values = c(24, 25), breaks = c("+" ,"-"), labels = c("Positive", "Negative")) + 
        ggforce::facet_col(~measure, space = "free", scales = "free", strip.position = "right") + 
        theme(axis.title.y = element_blank(), 
              axis.title.x = element_blank(),
              panel.grid.minor = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.y = ggtext::element_markdown(),
              legend.position = "bottom",
              legend.title = element_blank(),
              legend.box="vertical") 
    }


## Regional association plot
plot_regional <- function(tbl, curr_cluster){
  cols <- bind_rows(g_all$colours, c(color = "#5E4B56", morbidity = "combination"))
  plot_data <- tbl %>%
#    plot_data <- tmp %>%
    filter(cluster == curr_cluster) %>% 
    separate(lead_variantId, into = c("lead_chrom", "lead_pos", "lead_ref", "lead_alt"), remove = FALSE) %>%
    mutate(direction = NA) %>% 
    mutate(direction = ifelse(is.na(beta), direction, ifelse(beta<0, "-", "+"))) %>% 
    mutate(direction = ifelse(is.na(odds_ratio), direction, ifelse(odds_ratio<1, "-", "+"))) %>%
    mutate(lead_pos = as.numeric(lead_pos))
  
  limits = c(round((min(as.numeric(plot_data$lead_pos))-250000),-3),
             round((max(as.numeric(plot_data$lead_pos))+250000),-3))
  
  plot_data %>%  
    ggplot(., aes(x = lead_pos, y = -log(pval), shape = direction)) + 
    geom_point_interactive(aes(fill = morbidity, colour = morbidity, tooltip = trait_reported, data_id = lead_variantId), alpha = 0.7) + 
    geom_point_interactive(data = subset(plot_data, has_sumstats == TRUE), aes(tooltip = trait_reported, data_id = lead_variantId), colour = "black") + 
    scale_colour_manual(values = cols$color, breaks = cols$morbidity) + 
    scale_fill_manual(values = cols$color, breaks = cols$morbidity) + 
    scale_shape_manual(values = c(24, 25), breaks = c("+" ,"-"), labels = c("Positive", "Negative")) + 
    theme_bw() + 
    scale_x_continuous(limits = limits, breaks = scales::breaks_width(250000), labels = scales::label_number(scale = 1/1000000, suffix = "Mb", accuracy = 0.01)) + 
    scale_y_continuous(limits = c(0, NA)) + 
    labs(x = paste("Position on chromosome", unique(plot_data$lead_chrom)),
         y = expression(-log[10]("p-value"))) +
    theme(panel.grid = element_blank(),
          legend.position = "bottom",
          legend.title = element_blank(),
          legend.box="vertical")
  }


## UI
ui <- dashboardPage(
    
    dashboardHeader(title = "TargetAgeGenetics"),
    dashboardSidebar(disable = T),
    
    # Dashboard ----
    dashboardBody(useShinyjs(),
                  # Row containing target table and details ----
                  fluidRow(
                      # Main panel column ----
                      column(width = 6,
                             # Box for table output
                             box(
                                 width = 12,
                                 reactable::reactableOutput("cluster_table"),
                             )
                      ), # end of column
                      column(width = 6,
                             # Box for table output
                             box(
                                 width = 12,
                                 visNetworkOutput("cluster_graph")
                             ),
                             box(
                                 width = 12,
                                 girafeOutput("forest_plot", height = "auto")
                             ),
                             box(
                               width = 12,
                               girafeOutput("regional_plot", height = "auto")
                             )
                      ) # end of column
            ) # end of fluidRow
    ) # end of dashboard body
) # end of UI


# Define server logic 
server <- function(input, output) {
    
    sticky_style <- list(position = "sticky", left = 0, background = "#fff", zIndex = 1,
                         borderRight = "1px solid #eee")
    
    output$cluster_table <- reactable::renderReactable({
        tbl %>% 
            select(-id) %>% 
            mutate(hasColoc = ifelse(hasColoc == TRUE, 'Yes', NA)) %>% 
            mutate(has_sumstats = ifelse(has_sumstats == TRUE, 'Yes', NA)) %>% 
            select(cluster, n_morbidities, morbidity, specificDiseaseName, gene.symbol, L2G, hasColoc, distanceToLocus, everything()) %>% 
            reactable(groupBy = "cluster", 
                      # paginateSubRows = FALSE,
                      searchable = TRUE,
                      columns = list(
                          ## Make cluster id sticky
                          cluster = colDef(      
                              style = sticky_style,
                              maxWidth = 90,
                              headerStyle = sticky_style),
                          ## Aggregate comma-separated list of morbidities
                          morbidity = colDef(
                              name = "Morbidities",
                              aggregate = "unique",
                              width = 225),
                          ## Specific disease name
                          specificDiseaseName = colDef(name = "EFO disease label"),
                          ## Summarise the number of morbidities
                          n_morbidities = colDef(name = "", 
                                                 aggregate = "unique",
                                                 format = list(
                                                     aggregated = colFormat(suffix = " morbidities")
                                                 )),
                          distanceToLocus = colDef(name = "Distance to Locus",
                                                   aggregate = "min",
                                                   format = list(cell = colFormat(separators = TRUE),
                                                                 aggregated = colFormat(separators = TRUE))),
                          ## Summarise whether any of the GWAS have summary statistics
                          has_sumstats = colDef(name = "Has SumStats",
                                                aggregate = uniqueDropNA),
                          ## Summarise whether there is colocalisation
                          hasColoc = colDef(name = "Evidence of colocalisation", 
                                            aggregate = uniqueDropNA),
                          ## Round L2G
                          L2G = colDef(aggregate = "max",
                                       format = colFormat(digits = 2)), 
                          ## Link to study page
                          studyId = colDef(html = TRUE,
                                                  cell = JS("
    function(cellInfo) {
      // Render as a link
      var url = 'https://genetics.opentargets.org/study/'+ cellInfo.value
      return '<a href=\"' + url + '\" target=\"_blank\">' + cellInfo.value + '</a>'
    }
  ")),
                          ## Link to gene prioritisation page for the study+variant
                          lead_variantId = colDef(html = TRUE,
                                                  cell = JS("
    function(cellInfo) {
      // Render as a link
      var url = 'https://genetics.opentargets.org/variant/' + cellInfo.value
      return '<a href=\"' + url + '\" target=\"_blank\">' + cellInfo.value + '</a>'
    }
  ")),
                          ## Summarise the mapped genes and number of variants
                          gene.symbol = colDef(name = "Top Gene", 
                                               aggregate = "frequency",
                                               html = TRUE,

                                               cell = JS("
    function(cellInfo) {
      // Render as a link
      var url = 'https://genetics.opentargets.org/study-locus/'+ cellInfo.row['studyId'] + '/' + cellInfo.row['lead_variantId']
      return '<a href=\"' + url + '\" target=\"_blank\">' + cellInfo.value + '</a>'
    }
  ")
                                               )
                      ),
                      selection = "multiple", onClick = "select",
                      #      defaultSelected = 1,
                      theme = reactableTheme(
                          rowSelectedStyle = list(backgroundColor = "#eee", boxShadow = "inset 2px 0 0 0 #ffa62d")
                      ))
    })
    
    selected <- reactiveValues()
    observe({
        selected$plot_height = 1000
        selected$rows <- getReactableState("cluster_table", "selected")
        if (!is.null(selected$rows)){
            selected$clusters <- unique(tbl$cluster[selected$rows])
            selected$n_variants <- tbl_leads %>% 
                filter(cluster == selected$clusters[length(selected$clusters)]) %>% 
                select(cluster, id, lead_variantId) %>%
                nrow()
            req(!is.null(selected$rows))
            selected$forest_plot <- plot_forest(tbl_genes, selected$clusters[length(selected$clusters)])
            selected$regional_plot <- plot_regional(tbl_genes, selected$clusters[length(selected$clusters)])
            ## get number of variants on y axis to scale the figure size
            selected$n_breaks <- length(ggplot_build(selected$forest_plot)$layout$panel_scales_y[[1]]$breaks)
            selected$plot_height <- (selected$n_breaks * 0.25) + 0.75
            selected$forest_plot_interactive <- girafe(ggobj = selected$forest_plot,
                                                       height_svg = selected$plot_height)
            selected$regional_plot_interactive <- girafe(ggobj = selected$regional_plot)
        }
    })
    
    output$cluster_graph <- renderVisNetwork({
        if (!is.null(selected$rows)){
            view_cluster(g_all$cn, g_all, selected$clusters)
        }
    })
    
    
    observe({output$forest_plot <- renderGirafe({
        if (length(selected$rows)>0){
            selected$forest_plot_interactive
        }
    })
    })
    
    observe({output$regional_plot <- renderGirafe({
      if (length(selected$rows)>0){
        selected$regional_plot_interactive
      }
    })
    
})

}

# Run the application 
shinyApp(ui = ui, server = server)
