## process genage data
library(dplyr)
library(shiny)

## GenAge Human - human genes thought to be involved in ageing
# (mostly based on evidence from model organisms)
genage_human <-
  read.csv("data/HAGR/genAge/genage_human.csv",
           header = TRUE,
           stringsAsFactors = FALSE
  )

## GenAge model organisms
## model organism genes
genage_model <- read.csv("data/HAGR/genAge/genage_models.csv",
                           header = TRUE,
                           stringsAsFactors = FALSE
) %>%
  rename(GenAge.model.organism.symbol = symbol,
         GenAge.model.organism = organism,
         GenAge.model.id = GenAge.ID,
         GenAge.model.name = name,
         GenAge.model.organism.entrez.id = entrez.gene.id,
         GenAge.max.avg.lifespan.change  = avg.lifespan.change..max.obsv.,
         GenAge.lifespan.effect = lifespan.effect,
         GenAge.longevity.influence = longevity.influence)

## human orthologs of animal model age-related genes
## all of these are already in genage_human
## so we only need it to connect GenAge to the extra details
## on the model organism ortholog
orthologs <- read.csv("data/HAGR/genAge/genage_models_orthologs.csv",
                      stringsAsFactors = FALSE,
                      header = TRUE) %>%
  rename(entrez.gene.id = Entrez.ID,
         GenAge.model.organism = Model.Organism,
         GenAge.model.organism.entrez.id = Model.Organism.Entrez.ID) %>%
  select(-X, -Species, -Model.Organism.Symbol, -Symbol) # we will join on entrez id


## Combine the databases
genage <- genage_human %>%
  left_join(orthologs, by = c("entrez.gene.id")) %>%
  left_join(genage_model, by = c("GenAge.model.organism", "GenAge.model.organism.entrez.id")) %>% 
  rename(targetSymbol = symbol,
         GenAge.why = why)

save(genage, file = "data/analysis/genage.Rda")

cellage <- read.csv("data/HAGR/cellAge/cellAge1.csv", sep = ";") %>%
  rename_all( ~ paste0("CellAge.", .x)) %>% 
  rename(targetSymbol = CellAge.gene_name) 

save(cellage, file = "data/analysis/cellage.Rda")

