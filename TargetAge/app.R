library(shiny)
library(shinyjs)
library(shinyWidgets)
library(shinydashboard)
library(tidyverse)
library(DT)
library(plotly)
library(reactable)
library(ggtext)
library(ggiraph)
library(visNetwork)
library(igraph)
#library(shinyFiles)   ### For PDB viewer
#library(nglShiny)     ### For PDB viewer

## Prepared data
load("data/prepared_table.Rda")  # Target annotations
load("data/genage_legend.Rda")
load("data/hallmarks_legend.Rda")
load("data/cellage_legend.Rda")
load("data/diseases_with_associations.Rda")
load("data/targetage_geneids.Rda")
load("data/genetics_table.Rda")
load("data/graph_all_morbidities.Rda")

## Select columns to show in the output
tbl_output <- tbl_targets %>% 
  select(targetId,
         symbolSynonyms,
         targetSymbol,
         overallMultimorbidityCount,
   starts_with("icon"), 
         all_of(d)) 

## Plot a heatmap showing the genetic association score for each age-related phenotype for one or more targets
plot_multiheatmap <- function(associations, diseases, sort_nmorb = FALSE, save = FALSE) {
    textcol <- "grey40"
    df <- associations %>%
        dplyr::select(targetSymbol, targetName, all_of(diseases)) %>%
        pivot_longer(cols = all_of(diseases)) 
    df <- df %>%
        mutate(name = factor(name, levels = c(
            sort(unique(df$name[df$name != "longevity"])), "longevity")) ## alphabetical order for morbidities, but with longevity on the end
        )
    if (sort_nmorb == TRUE){
        nsum_order <-
            associations %>%
            arrange(-overallMultimorbidityCount) %>%
            pull(targetSymbol)
        df <- df %>%
            mutate(targetSymbol = factor(targetSymbol, levels = unique(nsum_order)))
    }
    df <- df %>% mutate(targetName = paste0(targetName, " (", targetSymbol,")"))
    heatmap <-
        ggplot(df, aes(x = name, y = targetSymbol, fill = value)) +
        geom_tile(colour = "black", size = 0.3) +
        coord_equal() +
        theme_bw(base_size = 10) +
        scale_y_discrete(expand = c(0, 0)) +
        scale_x_discrete(expand = c(0, 0)) +
        scale_fill_gradient(
            "score",
            limits = c(0,1),
            low = "white",
            high = "steelblue",
            na.value = 'gray90',
            guide = guide_colorbar(frame.colour = "black", frame.linewidth = 0.8)
        ) +
        guides(fill = FALSE) +
        theme(
            axis.text.x = element_text(
                colour = textcol,
                angle = 45,
                hjust = 1
            ),
            axis.text.y = element_text(vjust = 0.2, colour = textcol),
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
            panel.border = element_rect(colour = "black", size = 0.8)
        )
    if (save == TRUE){
        ggsave(heatmap, file = "genage_heatmap.png", width=7, height=4, dpi=900)
    }
    else{
        return(heatmap)
    }
}

#### A function to view visNetwork graph of genetic association cluster
view_cluster <- function(nodes, g, curr_cluster){
  cluster_nodes <- nodes %>% filter(cluster %in% curr_cluster) %>%
    mutate(color = ifelse(color == "black", "#5E4B56", color)) %>% 
    mutate(color.border = ifelse(has_sumstats == TRUE, "black", color),
           borderWidth = 2) %>%
    rename(color.background = color) 
  
  cluster_edges <- bind_rows(g$edges %>% mutate(width = 1), g$coloc_edges %>% mutate(width = 3)) %>% 
    filter((from %in% cluster_nodes$id) & (to %in% cluster_nodes$id)) 
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
    scale_x_continuous(limits = limits, labels = scales::label_number(scale = 1/1000000, suffix = "Mb", accuracy = 0.01)) + 
    scale_y_continuous(limits = c(0, NA)) + 
    labs(x = paste("Position on chromosome", unique(plot_data$lead_chrom)),
         y = expression(-log[10]("p-value"))) +
    theme(panel.grid = element_blank(),
          legend.position = "bottom",
          legend.title = element_blank(),
          legend.box="vertical")
}



### To rotate the headers ----
headerCallback <- c(
    "function(thead, data, start, end, display){",
    "  var $ths = $(thead).find('th');",
    "  $ths.css({'vertical-align': 'bottom', 'white-space': 'nowrap'});",
    "  var betterCells = [];",
    "  $ths.each(function(){",
    "    var cell = $(this);",
    "    var newDiv = $('<div>', {height: 'auto', width: '50px'});",
    "    var newInnerDiv = $('<div>', {text: cell.text()});",
    "    newDiv.css({margin: 'auto'});",
    "    newInnerDiv.css({",
    "      transform: 'rotate(180deg)',",
    "      'writing-mode': 'tb-rl',",
    "      'white-space': 'nowrap'",
    "    });",
    "    newDiv.append(newInnerDiv);",
    "    betterCells.push(newDiv);",
    "  });",
    "  $ths.each(function(i){",
    "    $(this).html(betterCells[i]);",
    "  });",
    "}"
)

## Custom JS aggregate function for reactable that drops NA
uniqueDropNA <- JS("function(values) {
      const isNotNA = function(x) { return x != null && x !== 'NA' }
      let uniqueValues = [...new Set(values)]
      return uniqueValues.filter(isNotNA).join(', ')  
    }")

## Options for GO filtering
go_options <- go_legend %>% 
    unite(label, c(goIcon, goHallmark), sep = " ", remove = FALSE)
    
## Initialise selected rows
s <- NULL

## To do: phenodigm
## Add mouse icon to phenodigm
#phenodigm <- 
#    phenodigm %>% 
#    mutate(phenodigm.mouse_model = paste("<img src=\"noun_Mouse_3551360.png\" height=\"20\" data-toggle=\"tooltip\" data-placement=\"right\" title=\"", human.phenotype, "\"></img>"))


## Colour and values for table colour formatting
brks <- seq(0.05, 1, 0.05)
clrs <- colorRampPalette(c("white", "#6baed6"))(21)

table_explanation1 <- "The table below shows each target's Open Targets genetic association score with each age-related diseases/trait, as well as the multimorbidity count (total number of associated morbidities). Icons indicate any Hallmarks of Ageing that the target is involved in according to its GO annotations." 
table_explanation2 <- "Select a target for more detailed information, and to view a structure if one is available. Much more information about a target and its association with disease, including drug and clinical trial information, can be found via the link to its Open Targets Platform profile on the right."
genage_explanation <- p(a(href = "https://genomics.senescence.info/genes/human.html", "GenAge"), " is an expert-curated database of age-related genes. The GenAge ID number links to the GenAge entry for relevant targets, and icons indicate whether animal model orthologs have experimental evidence of pro-longevity or anti-longevity effects.")
cellage_explanation <- p(a(href = "https://genomics.senescence.info/cells/", "CellAge"), " is a manually-curated database of genes associated with cell senescence through genetic experiments in human cell types. The icons indicate whether the gene was found to induce or inhibit senescence.")
all_explanation <- "All targets with genetic association (>= 0.5) with a cluster containing two or more phenotypes of interest. Adjust filters using the dropdown options above."
targetage_explanation <- "Targets implicated in a genetic association cluster involving two or more age-related diseases/traits."
#lof_explanation <- "Targets with predicted loss-of-function (LOF) variants that are protective against at least one morbidity of interest." 
longevity_explanation <- "Targets genetically associated with longevity and at least one morbidity."

# Define UI for dataset viewer app ----
ui <- dashboardPage(
    # App title ----
    dashboardHeader(title = "TargetAge"),
    dashboardSidebar(disable = T),
    
    # Dashboard ----
    dashboardBody(useShinyjs(),
                  # Row containing target table and details ----
                  fluidRow(
                      # Main panel column ----
                      column(width = 9,
                             # Box for table output
                             box(
                                 width = 12,
                                 
                                 ## Dropdown filtering options ----
                                 dropdownButton(
                                     # Toggle AND or OR ----
                                     switchInput(
                                         inputId = "and_or",
                                         onLabel = "AND",
                                         offLabel = "OR",
                                         label = "Filter type",
                                         labelWidth = "80px"
                                     ),

                                     ### Filtering options ----
                                     
                                     # Input: morbidities to show ----
                                     multiInput(
                                         inputId = "show_morbs",
                                         label = "Morbidities:",
                                         choices = NULL,
                                         choiceNames = sort(d),
                                         choiceValues = sort(d),
                                         width = "100%"
                                     ),
                                     
                                     # Input: hallmarks to show ----
                                     multiInput(
                                         inputId = "show_hallmarks",
                                         label = "Hallmarks:",
                                         choices = NULL,
                                         choiceNames = go_options$label,
                                         choiceValues = go_options$goHallmark
                                     ),
                                     
                                     # Input: Numeric entry for minimum multimorbidity count ----
                                     numericInput(
                                         inputId = "minmmcount",
                                         label = "Minimum multimorbidity count:",
                                         value = 2
                                     ),
                                     
                                     # # Input: Only targets without PDB structure ----
                                     # radioGroupButtons(
                                     #     inputId = "pdbfilter",
                                     #     label = "Filter by solved crystal structures",
                                     #     choices = c(
                                     #         "All" = "all",
                                     #         "With structures" = "haspdb",
                                     #         "No structures" = "nopdb"
                                     #     ),
                                     #     selected = "all"
                                     # ),
                                    
                                     # Button to show filter settings ----
                                     circle = TRUE,
                                     status = "danger",
                                     icon = icon("sliders"),
                                     width = "50%",
                                     tooltip = tooltipOptions(title = "Filter targets and morbidities")
                                 ),
                                 #----
                                 
                                 br(),
                                 tags$div(tags$p(table_explanation1),
                                          tags$p(table_explanation2)),
                                 br(),
                                 
                                 tabsetPanel(id = "targetlist",
                                             selected = "all",
                                             # Tab 1 ----
                                             tabPanel("All targets",
                                                      tags$div(
                                                          br()
                                                      ), 
                                                      # Output: HTML table with associations for selected diseases
                                                      DT::dataTableOutput("view"),
                                                      value = "all"),
                                             # Tab 2 ----
                                              tabPanel("TargetAge",
                                                       tags$div(br()),
                                                       # Output: HTML table with associations for targetage set
                                                       DT::dataTableOutput("targetage"),
                                                       value = "targetage")
                                             # Tab 6 ----
                                             # tabPanel("All triaged",
                                             #          tags$div(
                                             #              br()),
                                             #          #    tags$p("Triage notes compiled by Dundee and MDC for targets of interest.")), 
                                             #          # Output: Triage information)
                                             #          br(),
                                             #          DT::dataTableOutput("triagetbl"),
                                             #          value = "triagetbl")
                                             
                                             # End of tabs ----
                                 ), # end of tabsetPanel
                             ), # end of box
                             ## Box for plots -----
                           #  box(
                           #      width = 12,
                            #     tabsetPanel(
                                     # Output: heatmap of selected target(s)
                                    # tabPanel("Heatmap", plotOutput("heatmap", height = "700px"))
                                     # # Output: PCA
                                     # tabPanel("PCA", plotlyOutput("pcaplot", height = "700px")),
                             #    ) # end of tabsetPanel
                            # ) #end of box
                      ), ## end of column
                      column(width = 3,
                             box(
                                 width = 12,
                                 title = "Target information",
                                 "Select a target to see more detailed information.",
                                 tabsetPanel(
                                     tabPanel("Details",
                                              tableOutput("targetdetails")),
                                      ### View top hit PDB structure
                                      tabPanel(
                                          "Structure",
                                          br(),
                                     #     h3(textOutput("pdb_code")),
                                     #     br(),
                                     #     #  actionButton("fitButton", "Reset position"),
                                     #     actionButton('pdb_go', 'Load PDB'),
                                     #     hr(),
                                     #     nglShinyOutput('nglShiny')
                                      )
                                 ), # end of tabset panel
                                 footer = htmlOutput("open_ot"),
                                 style = "overflow-y:scroll; max-height: 130vh"
                             ), # end of box
                             # box(
                             #     width = 12,
                             #     title = "Triage Legend",
                             #     DT::dataTableOutput("triagelegend")
                             # ),
                             box(
                                 width = 12,
                                 title = "Hallmarks Icons",
                                 tableOutput("icon_legend")
                             ),
                             box(
                                 width = 12,
                                 title = "GenAge & CellAge Icons",
                                 genage_explanation,
                                 tableOutput("genage_legend"),
                                 cellage_explanation,
                                 tableOutput("cellage_legend")
                             )
                      ), # end of fluidRow
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
                  ) # end of column 2 
    ) # end of dashboard body
) # end of UI

# Define server logic to summarize and view selected dataset ----
server <- function(input, output, session) {
    min_gen_score = 0.5
    
    # some global datatable options 
    dt.options <- list(
        paging = TRUE,
        pageLength = 25,
        scrollX = TRUE,
        scrollY = TRUE,
        autoWidth = TRUE,
        server = FALSE,
        dom = "Bfrtip",
        search = list(regex = TRUE),
        headerCallback = JS(headerCallback),
        #                     fixedColumns = list(leftColumns = 1),
        columnDefs = list(list(
            width = '0px',
            targets = "_all",
            className = "dt-center"
        ),
        list(targets = c(0,1), 
             visible = FALSE)))
    
    ## selected morbidities to show
    morbs <- reactive({
        if (length(input$show_morbs) == 0) {
            d
        } else {
            input$show_morbs
        }
    })
    
    # test_pdb_structure <- function(n.pdbs, pdbchoice) {
    #     if (pdbchoice == "haspdb") {
    #         n.pdbs > 0 
    #     } else if (pdbchoice == "nopdb") {
    #         n.pdbs == 0
    #     } else {
    #         TRUE
    #     }
    # }
    
    min_gen_score = 0.05 # make this a slider
    
    ## update list according to settings
    all_dat <- reactive({
        selected_morbidities <- morbs()
        selected_dat <- tbl_output %>%
           filter(overallMultimorbidityCount >= input$minmmcount) %>%
           arrange(-overallMultimorbidityCount)
       # Filter hallmarks if any are selected
       if (length(input$show_hallmarks) > 0) {
           selected_dat <- selected_dat %>%
               filter(str_detect(icon.goHallmarks, 
                                 paste(input$show_hallmarks, 
                                       collapse = "|")))
       }
       if (input$and_or) {
           ## if AND targets must be associated with all morbidities
           selected_dat %>% filter_at(selected_morbidities, all_vars(. > min_gen_score))
       } else {
           ## if OR then targets must be associated with any selected morbidity
           selected_dat %>% filter_at(selected_morbidities, any_vars(. > min_gen_score))
       }
    })
    
    ## Table of targets and genetic associations with all morbidities
    output$view <- DT::renderDataTable({
      all_dat() %>%
            datatable(
                options = dt.options,
                class = 'row-border',
                escape = FALSE,
                rownames = FALSE,
                filter = "none",
                caption = all_explanation
            ) %>%
            formatStyle(d, backgroundColor = styleInterval(brks, clrs)) 
    })
    
    output$targetage <- DT::renderDataTable({
      all_dat() %>%
        filter(targetId %in% targetage) %>% 
        datatable(
          options = dt.options,
          class = 'row-border',
          escape = FALSE,
          rownames = FALSE,
          filter = "none",
          caption = all_explanation
        ) %>%
        formatStyle(d, backgroundColor = styleInterval(brks, clrs)) 
    })
      
    
    ## Genetics 
    sticky_style <- list(position = "sticky", left = 0, background = "#fff", zIndex = 1,
                         borderRight = "1px solid #eee")
    
    get_selected_gene_tbl <- function(df, selected_targets){
        curr <- tbl_genes %>% 
          filter(str_detect(targetIds, 
                        paste(selected_targets, 
                              collapse = "|")))
    }
    
    output$cluster_table <- reactable::renderReactable({
      selected$tbl_genes %>%  
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
                    ## Aggregate number of communities in the cluster subgraph 
                    communities = colDef(
                      aggregate = "max"
                    ),
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
                  defaultSorted = list(communities = "asc", n_morbidities = "desc"),
                  theme = reactableTheme(
                    rowSelectedStyle = list(backgroundColor = "#eee", boxShadow = "inset 2px 0 0 0 #ffa62d")
                  ))
    })
    
    
    selected <- reactiveValues()
    observe({
      selected$plot_height = 1000
      selected$rows <- getReactableState("cluster_table", "selected")
      selected$targets <- selected_targets()
      selected$tbl_genes <- get_selected_gene_tbl(tbl_genes, selected$targets)
       if (!is.null(selected$rows)){
         selected$clusters <- tail(unique(selected$tbl_genes$cluster[selected$rows]),1)
         selected$n_variants <- selected$tbl_genes %>% 
           filter(cluster == selected$clusters[length(selected$clusters)]) %>% 
           select(cluster, id, lead_variantId) %>%
           nrow()
         req(!is.null(selected$rows))
         selected$forest_plot <- plot_forest(selected$tbl_genes, selected$clusters[length(selected$clusters)])
         selected$regional_plot <- plot_regional(selected$tbl_genes, selected$clusters[length(selected$clusters)])
      #   ## get number of variants on y axis to scale the figure size
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
    
    # output$triagetbl <- renderDataTable(
    #     triage %>%
    #         select(Colour, Meaning, everything()) %>% 
    #         rename(target = symbol) %>% 
    #         datatable(
    #             options = list(
    #                 scrollX = TRUE,
    #                 scrollY = TRUE,
    #                 autoWidth = TRUE,
    #                 server = FALSE,
    #                 #          dom = "Bfrtip",
    #                 #          buttons = c('csv', 'excel'),
    #                 search = list(regex = TRUE),
    #                 columnDefs = list(
    #                     list(
    #                         #            targets = c(2, 3, 16, 17, 18),
    #                         targets = c(2, 3, 15, 16, 17),
    #                         className = "dt-center"),
    #                     list(targets = c(0,1),
    #                          visible = FALSE))
    #             ),
    #             caption = triage_explanation, 
    #             rownames = FALSE,
    #             escape = FALSE)
    #     %>% 
    #         formatTriage()
    # )
    # 
    # 
    
    ## Selected rows
    selected_rows = reactive({
        if (input$targetlist == "all"){
            input$view_rows_selected
        } else if (input$targetlist == "targetage") {
            input$targetage_rows_selected 
        } else {
            input$view_rows_selected
        }
    })
    
    ## First 10 rows
    default_rows = reactive({
        if (input$targetlist == "all"){
            all_dat() %>%
                slice(1:10) %>%
                pull(targetId)
        } else if (input$targetlist == "targetage") {
            all_dat() %>%
                filter(targetId %in% targetage) %>%
                slice(1:10) %>%
                pull(targetId)
        } else {
            all_dat() %>%
                slice(1:10) %>%
                pull(targetId)      
          }
    })
    
    ## All selected targets
    selected_targets <- reactive({
        s <- selected_rows()
        if (length(s)) {
            if (input$targetlist == "targetage") {
                all_dat() %>%
                    filter(targetId %in% targetage) %>%
                    slice(s) %>%
                    pull(targetId)
            } else {
                all_dat() %>%
                    slice(s) %>%
                    pull(targetId)
            }
        } else {
            NULL
        }
    })
    
    ## OT URL for selected target
    ot_url <- reactive({
        st <- selected_targets()
        if (!is.null(st)){
          paste0("https://platform.opentargets.org/target/", st[length(st)])
        } else {
          "https://platform.opentargets.org/"
        }
    })
    
    output$open_ot <- renderUI({
        d_url <- ot_url()
        tags$a(href = d_url, "See more at Open Targets")
    })
    
    ## Get target details for most recently selected target
     selected_targetdetails <- reactive({
         s <- selected_rows()
         if (length(s)) {
             st <- selected_targets()
    #         # filter for the full details
             selected_target <-
                 tbl_targets %>%
                 filter(targetId == st[length(st)]) %>%
                 slice(1)
         }
     })
    
     output$targetdetails <- renderTable ({
         details <- selected_targetdetails()
         if (!is.null(details)) {
           details$details_output %>% t()
         }
     },
     spacing = "s",
     rownames = TRUE,
     colnames = FALSE,
     sanitize.text.function = function(x) x)
    
    output$heatmap <- renderPlot({
        s <- selected_targets()
         if (length(s) < 1) {
            s <- default_rows()
         }
            tbl_targets %>% filter(targetId %in% s) %>%
                plot_multiheatmap(diseases = d)
    })
    

    # output$pcaplot <- renderPlotly({
    #     s <- selected_targets()
    #     s_genes <- a %>% filter(targetId %in% s) %>% pull(targetSymbol)
    #     ggplot(pca.df,
    #            aes(
    #                x = PC1,
    #                y = PC2,
    #                colour = as.character(cluster),
    #                fill = as.character(cluster),
    #                label = target
    #            )) +
    #         geom_point() +
    #         geom_point(data = subset(pca.df, target %in% s_genes), aes(x = PC1, y = PC2), shape = 21, colour = "black") + 
    #         theme_bw() +
    #         theme(panel.grid = element_blank(), legend.position = "none")
    #     ggplotly(tooltip = "label")
    # })
    # 
    # pca_selected <- reactive({
    #   p <- event_data("plotly_selected")
    #   if (!is.null(p)){
    #     p_selected <-
    #       p %>%
    #       rename(PC1 = x, PC2 = y) %>%
    #       inner_join(pca.df) %>%
    #       pull(target)
    #     s <- gen_dat() %>%
    #       pull(target.gene_info.name) %>%
    #       which(. %in% p_selected)
    #     print(s)
    #   }
    # })
    
    ## Hallmarks legend
    output$icon_legend <- renderTable ( 
        go_legend,
        spacing = "s",
        rownames = FALSE,
        colnames = FALSE,
        sanitize.text.function = function(x) x
    )
    
    ## Genage legend
    output$genage_legend <- renderTable ( 
        genage_legend,
        spacing = "s",
        rownames = FALSE,
        colnames = FALSE,
        sanitize.text.function = function(x) x
    )
    
    ## CellAge legend
    output$cellage_legend <- renderTable ( 
        cellage_legend,
        spacing = "s",
        rownames = FALSE,
        colnames = FALSE,
        sanitize.text.function = function(x) x
    )
    
    headerCallbackRemoveHeaderFooter <- c(
        "function(thead, data, start, end, display){",
        "  $('th', thead).css('display', 'none');",
        "}"
    )
    
    # output$triagelegend <- renderDataTable ( 
    #     triage %>% 
    #         select(Colour, Meaning) %>% 
    #         unique() %>%
    #         mutate(Colour = "") %>% 
    #         datatable(
    #             options = list(
    #                 dom = "t",
    #                 ordering=FALSE, 
    #                 paging = FALSE,
    #                 searching = FALSE, 
    #                 #     autoWidth = TRUE,
    #                 headerCallback = JS(headerCallbackRemoveHeaderFooter),
    #                 server = FALSE),
    #             selection = 'none',
    #             callback = JS("$('table.dataTable.no-footer').css('border-bottom', 'none');"),
    #             class = 'row-border',
    #             extensions = "Buttons",
    #             escape = FALSE,
    #             rownames = FALSE,
    #             filter = "none"
    #         ) %>%
    #         formatStyle("Colour", valueColumns = "Meaning", backgroundColor = styleEqual(unique(triage$Meaning), triage_colours))
    # )
    
    
    ### PDB viewer
    #   observeEvent(input$fitButton, ignoreInit = TRUE, {
    #      fit(session)
    #    })
    
    # target_pdbs <- reactive({
    #     st <- selected_targets()
    #     all_target_pdbs %>% 
    #         filter(targetId == st[length(st)])
    # })
    # 
    # output$pdb_code <- renderText({
    #     if (length(selected_rows()) > 0) { 
    #         pdbs <- target_pdbs()
    #         pdb_id <- pdbs$id[1]
    #         if (!is.na(pdb_id)) {
    #             paste("Structure: ", pdb_id)
    #         } else {
    #             "No PDB structure found"
    #         }
    #     } else {
    #         "Select a target to view structure"
    #     }
    # })
    # 
    # observeEvent(input$pdb_go, ignoreInit = TRUE, {
    #     pdbs <- target_pdbs()
    #     pdb_id <- pdbs$id[1]
    #     selection_value <- pdbs$value.chains
    #     residues <- strsplit(selection_value, "=")[[1]][2]
    #     chain <- strsplit(strsplit(selection_value, "=")[[1]][1], "/")[[1]][1]
    #     selection_string <- paste0(":",chain," and ", residues)
    #     session$sendCustomMessage(type = "removeAllRepresentations", message =
    #                                   list())
    #     session$sendCustomMessage(type = "setPDB",
    #                               message = list(
    #                                   uri = sprintf(
    #                                       "https://files.rcsb.org/download/%s.pdb",
    #                                       toupper(pdb_id),
    #                                       selection = selection_string
    #                                   )
    #                               ))
    #     session$sendCustomMessage(type = "showSelection",
    #                               message = list(
    #                                   selection = selection_string,
    #                                   representation = "cartoon"
    #                               ))
    #     fit(session)
    #     #     print(selection_string)
    # })
    # options <- list(pdbID = '')
    # output$nglShiny <- renderNglShiny(nglShiny(options, 5, 5))
}

# Create Shiny app ----
shinyApp(ui = ui, server = server)


