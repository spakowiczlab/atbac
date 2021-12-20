format_pathway_dat <- function(exprs, meta, generoles){
  ctcoms <- max(str_count(generoles$Genes, ",")) + 1
  lgenes.l <- generoles %>%
    separate(Genes, sep = ",", into = paste0("gene", 1:ctcoms)) %>%
    gather(-Role, key = "gcount", value = "Gene") %>%
    filter(!is.na(Gene)) %>%
    dplyr::select(-gcount)
  
  exprs.l <- exprs %>%
    gather(-Gene, key = "sample", value = "counts")
  
  boxin <- lgenes.l %>%
    left_join(exprs.l) %>%
    left_join(meta)
  
  return(boxin)
}