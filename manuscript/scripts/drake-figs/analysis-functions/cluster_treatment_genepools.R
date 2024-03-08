cluster_treatment_genepools <- function(comptreat){
  comptreat <- comptreat[,-6]
  colnames(comptreat) <- c("Gene", "Airspace", "VAT", "Endotoxin", "Exercise", 
                           # "LungCancer",
                           "Sepsis", "Synovial",  "Active TB")
  
  fordist <- comptreat %>%
    column_to_rownames(var = "Gene") %>%
    as.matrix() %>%
    t()
  
  hold <- hclust(dist(fordist))
  
  hold <- dendro_data(hold)
  
  upsetdat <- comptreat %>%
    select(as.character(hold$labels$label))
  
  return(upsetdat)
}
