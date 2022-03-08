# studies + variants
studies_variants <- g_all$cn %>% filter(n_morbidities > 1) %>% ungroup() %>% select(studyId, lead_variantId) %>% unique()

if(0) {
  library(ghql)
  library(jsonlite)
  
  initialise_api <- function(){
    cli <- GraphqlClient$new(
      url = "https://genetics-api.opentargets.io/graphql"
    )
    return(cli)
  }
  
  initialise_queries <- function(){
    qry <- Query$new()
    qry$query('l2g_query', 'query l2gQuery($studyId: String!, $variantId: String!){
  studyLocus2GeneTable(studyId: $studyId, variantId: $variantId){
    rows {
      gene {
        id
        symbol
      }
      hasColoc
      yProbaModel
      yProbaDistance
      yProbaInteraction
      yProbaMolecularQTL
      yProbaPathogenicity
      distanceToLocus
    }
  }
  variantInfo(variantId: $variantId){
      mostSevereConsequence
  }
}')
    return(qry)
  }
  
  ## Do the API call
  fetch_l2g <- function(df, variables){
    result <- fromJSON(cli$exec(qry$queries$l2g_query, variables, flatten = TRUE))$data
    
    l2g_result <- result$studyLocus2GeneTable %>% bind_cols(result$variantInfo) %>% bind_cols(df) 
    return(l2g_result)
  }
  
  #########
  
  cli <- initialise_api()
  qry <- initialise_queries()
  
  ## split the data frame into smaller chunks (1000 rows)
  ## I don't really know if we need to do this but just in case
  ## I'm scared of breaking the server with too many successive calls
  
  n <- 100
  nr <- nrow(studies_variants)
  studies_variants_split <-
    studies_variants %>% 
    split(., rep(1:ceiling(nr / n), each = n, length.out = nr))
  
  # Somewhere to hold the results
  l2g_all <- vector(mode = "list", length = length(studies_variants_split))
  
  ## Do the first chunk on its own to check 
  l2g_all[[1]] <-  studies_variants_split[[1]]  %>% 
    group_by(studyId, lead_variantId) %>% 
    group_split() %>% 
    ## API call for each studyID + variantID
    purrr::map(~fetch_l2g(df = ., variables = list(studyId = .$studyId, variantId = .$lead_variantId))) %>%
    bind_rows() 
  
  ## Do all the other chunks
  for (i in seq(2, length(l2g_all))){
    print(i)
    l2g_all[[i]] <- studies_variants_split[[i]] %>% 
      group_by(studyId, lead_variantId) %>% 
      group_split() %>% 
      ## API call for each studyID + variantID
      purrr::map(~fetch_l2g(df = ., variables = list(studyId = .$studyId, variantId = .$lead_variantId))) %>%
      bind_rows() 
    Sys.sleep(3)
  }
  
  l2g_all_joined <- l2g_all %>% bind_rows() %>% jsonlite::flatten() %>% filter(yProbaModel >= 0.05)
  save(l2g_all_joined, file="ltg_all.Rda")
  
} else {
  load("ltg_all.Rda")
}
