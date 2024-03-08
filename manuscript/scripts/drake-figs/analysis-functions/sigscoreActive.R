sigscoreActive <- function(){
  neutmetas <- pull_iso_metas()
  neutcounts <- pull_iso_counts(neutmetas)
  labelled.genes = read.csv("../../data/from-alecia/pathgenes.csv", stringsAsFactors = F)
  
  pathsigs <- str_split(labelled.genes$Genes, ",")
  names(pathsigs) <-labelled.genes$Role
  
  calculateAvgZScore <-
    function (gene_matrix, genes) 
    {
      gene_missing <- tmesig::checkGenes(gene_matrix$Gene.Symbol, score = NULL, 
                                         expected.genes = genes)
      if (length(gene_missing) > 0) {
        mis.genes.short <- paste(gene_missing, collapse = ",")
        warn.message <- paste("Genes missing from set:", mis.genes.short)
        print(warn.message)
      }
      avg_z_score <- gene_matrix %>% tidyr::pivot_longer(-Gene.Symbol, 
                                                         names_to = "sample", values_to = "counts") %>% dplyr::group_by(Gene.Symbol) %>% 
        dplyr::mutate(Average = mean(counts), std_dev = sd(counts)) %>% 
        dplyr::ungroup() %>%
        dplyr::mutate(z_score = (counts - 
                                   Average)/std_dev) %>%
        dplyr::filter(Gene.Symbol %in% genes & std_dev != 0 ) %>%
        dplyr::group_by(sample) %>% dplyr::summarize(avg_z_score = mean(z_score))
      return(avg_z_score)
    }
  
  calcManyZScores <- function(inputdat){
    tmp <- lapply(names(pathsigs), function(x) calculateAvgZScore(inputdat, pathsigs[[x]]) %>%
                    rename(!!x := "avg_z_score"))
    
    res.df <- purrr::reduce(tmp, full_join)
  }
  
  neutcounts.form <- lapply(neutcounts, function(x) 
    x %>%
      rename("Gene.Symbol" = Gene) %>%
      mutate_if(is.integer, as.double))
  
  zdf.ls <- lapply(neutcounts.form, calcManyZScores)
  
  celllabs <- as.data.frame(cbind(
    Neutrophil= c("bldvat", "endotox", "exercise", "sepsis", "tb", "syno", "airspace"),
    neutlabs = c("VAT", "Endotoxin", "Exercise", "Sepsis", "Active TB", "Synovial", "Airspace")
  ))

  treatment.scores <- lapply(names(neutmetas), function(x)
    neutmetas[[x]] %>%
      filter(source == "treatment") %>%
      select(sample) %>%
      mutate(Neutrophil = x) %>%
      left_join(zdf.ls[[x]])) %>%
    bind_rows() %>%
    left_join(celllabs)
  
  return(treatment.scores)
  
  # This section used to check significant differences between each neutrophil and VAT
  # 
  # toKWtest <- function(nname){
  #   tmp <- treatment.scores %>%
  #     filter(neutlabs %in% c("VAT", nname))
  #   
  #   kw.res <- lapply(names(pathsigs), function(x) kruskal.test(tmp[[x]], tmp[["neutlabs"]]) %>%
  #                      broom::tidy() %>%
  #                      mutate(pathway = x)) %>%
  #     bind_rows()
  #   
  #   return(kw.res)
  # }
  # 
  # neuts <- unique(treatment.scores$neutlabs)[-1]
  # 
  # kws <- lapply(neuts, toKWtest)
  # names(kws) <- neuts
}