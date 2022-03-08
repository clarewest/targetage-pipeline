### Triage
library(readxl)

## Map Charlotte's pathway annotations to full names and icons
triage_map <- data.frame(
  abbreviation = c("\\bS\\b",
                   "\\bM\\b",
                   "\\bT\\b",
                   "\\bINF\\b",
                   "\\bP\\b",
                   "\\bDDR\\b",
                   "\\bUPR\\b"),
  pathway = c("cellular senescence",
              "metabolism",
              "transcription",
              "inflammation",
              "proteostasis",
              "DNA damage repair",
              "unfolded protein response (UPR)"),
  icon = c("<i class=\"fab fa-canadian-maple-leaf\"",
           "<i class=\"fas fa-recycle\"",
           "<i class=\"fas fa-pen-alt\"",
           "<i class=\"fas fa-fire-alt\"",
           "<i class=\"fa fa-cut\"",
           "<i class=\"fa fa-dna\"",
           "<i class=\"fa fa-cut\"")
) %>%
  mutate(icon.plain = paste0(icon, "></i>"),
         icon.tooltip = paste0(icon, " data-toggle=\"tooltip\" data-placement=\"right\" title=\"", pathway, "\"></i>"))


triage_mapper = tibble::deframe(triage_map[c(1,5)])

triage <- read_excel("~/work/shinytargetage/data/combined_triage.xlsx") %>%
  rename(symbol = `Gene name`) %>%
  mutate(Pathways = stringr::str_replace_all(string = Pathways,
                                             pattern= triage_mapper)) %>%
  mutate(Pathways = stringr::str_remove_all(Pathways, ",")) %>%
  select(-Biomarker) %>% 
  rename(Biomarker = `Biomarker notes`) %>%
  select(symbol, 
         Colour,
         Meaning,
         Pathways, 
         Indication,
         `Mode of action`,
         `Connection to ageing`,
         `Toxicity issues`,
         `Toxicity notes`,
         Considerations,
         Recommendation,
         `Link to hypothesis`,
         `Validated link to ageing models/pathways`,
         `Reagent and assay availability`,
         `Competitor landscape`,
         Biomarker, 
         Structure, `Chemical matter`)
save(triage, file="~/work/shinytargetage/data/triage.Rda")
