choose_enriched_genes <- function(genemat){
  findgenes <- genemat %>%
    dplyr::select(NAME, Blood.neutrophil, V.neutrophil) %>%
    mutate(difference = Blood.neutrophil - V.neutrophil)
  
  enriched.blood <- findgenes %>%
    arrange(desc(difference))
  
  enriched.blood <- enriched.blood[1:10, ]
  
  enriched.at <- findgenes %>%
    arrange(difference)
  
  enriched.at<- enriched.at[1:10, ]
  
  
  mostenriched <- bind_rows(enriched.blood, enriched.at)
  mostenriched$celltype <- ifelse(mostenriched$difference > 0, "blood", "at")
  return(mostenriched)
}