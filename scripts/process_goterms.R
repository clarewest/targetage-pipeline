library(GO.db)
library(tidyverse)

high_go <- read.csv("data/goterm_list.csv", col.names = c("name", "goHallmark", "goHallmarkTerm"))

go_offspring <- as.list(GOBPOFFSPRING) %>% enframe() %>% unnest(value) 

go_hallmarks <- high_go %>% 
  left_join(go_offspring, by = "name") %>%  
  rowwise() %>% 
  mutate(term = ifelse(is.na(value), NA, Term(GOTERM[[value]]))) %>%
  ungroup() %>%
  dplyr::rename(goHallmarkId = `name`, goId = value, goTerm = term) %>%
  bind_rows(high_go %>% 
              dplyr::rename(goHallmarkId = `name`) %>% 
              mutate(goId = goHallmarkId,
                     goTerm = goHallmarkTerm)) %>%
  filter(!is.na(goId)) # terms with no children have NA from the join

go_hallmarks_genes <- AnnotationDbi::select(org.Hs.eg.db,
                                    keys = go_hallmarks$goId %>% unique(),
                                    columns = c("ENTREZID", "SYMBOL", "ENSEMBL"),
                                    keytype = "GO") %>%
  as_tibble #%>%
 # filter(!duplicated(ENTREZID))

write.csv(go_hallmarks, file = "data/full_goterm_list.csv", row.names = FALSE)
write.csv(go_hallmarks_genes, file = "data/full_goterm_genes.csv", row.names = FALSE)
