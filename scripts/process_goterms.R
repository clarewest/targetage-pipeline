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
              rename(goHallmarkId = name) %>% 
              mutate(goId = goHallmarkId,
                     goTerm = goHallmarkTerm)) %>%
  filter(!is.na(goId)) # terms with no children have NA from the join

write.csv(go_hallmarks, file = "data/full_goterm_list.csv", quote = FALSE, row.names = FALSE)

