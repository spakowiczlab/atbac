nmds_neuts_pathways <- function(davgroups, reactkey){
  path.pvals <- lapply(names(davgroups), function(x)  pathwayDescriptionAssign(davgroups[[x]], reactkey) %>% 
                         dplyr::select(Category, Term, Description, PValue) %>%
                         mutate(source = x,
                                FDR = p.adjust(PValue, method = "fdr"))) %>%
    bind_rows()
  
  
  path.bi.in <- path.pvals %>%
    dplyr::select(source, FDR, Description) %>%
    filter(!is.na(Description)) %>%
    filter(Description != "") %>%
    mutate(FDR = 1-FDR) %>%
    spread(key = "Description", value = 'FDR') %>%
    gather(-source, key = "Description", value = "FDR") %>%
    mutate(FDR = ifelse(is.na(FDR), 0, FDR)) %>%
    spread(key = "Description", value = 'FDR')
  
  path.novar <- lapply(path.bi.in[,-1], function(x) var(x))
  path.novar <- names(subset(path.novar, path.novar == 0))
  
  path.nmds.in <- path.bi.in %>%
    dplyr::select(-path.novar)%>%
    column_to_rownames(var = "source")
  
  set.seed(112358)
  path.nmds <- metaMDS(path.nmds.in)
  
  path.nmds.points <- path.nmds$points %>%
    as.data.frame() %>%
    rownames_to_column(var = "pointname") %>%
    mutate(pointtype = "source")
  
  path.nmds.species <- path.nmds$species %>%
    as.data.frame() %>%
    rownames_to_column(var = "pointname") %>%
    mutate(pointtype = "pathway")
  
  
  path.nmds.ggdat <- bind_rows(path.nmds.points, path.nmds.species) %>%
    mutate(pointlab = ifelse(pointtype == "source", pointname,
                             "")
           # ifelse(MDS1 > 0.5, pointname,
           #        ifelse(MDS2 > 1, pointname, "")))
    )
  
  return(path.nmds.ggdat)
}

pathwayDescriptionAssign <- function(davidoutput, reactkey){
  db <- c('REACTOME_PATHWAY', 'KEGG_PATHWAY', 'BIOCARTA', 'BBID')
  tmp.ls <- list()
  for(d in db){
    tmp.ls[[d]] <- davidoutput %>%
      filter(Category == d)
  }
  
  tmp.ls$REACTOME_PATHWAY <- tmp.ls$REACTOME_PATHWAY %>%
    mutate(Term = gsub("(*.):.*","\\1", Term)) %>%
    left_join(reactkey)
  
  tmp.ls$KEGG_PATHWAY <- tmp.ls$KEGG_PATHWAY %>%
    separate(Term, into = c("Term", "Description"), sep = ":")
  
  tmp.ls$BIOCARTA <- tmp.ls$BIOCARTA %>%
    separate(Term, into = c("Term", "Description"), sep = ":")
  
  tmp.ls$BBID <- tmp.ls$BBID %>%
    separate(Term, into = c("Term", "Description"), sep = ".")
  
  tmp <- bind_rows(tmp.ls) %>%
    arrange(PValue)
  return(tmp)
}
