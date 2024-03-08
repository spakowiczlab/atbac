nmds_neuts_counts <- function(countls, metals){
  datsource <- c("bldvat", "endotox", "exercise", "sepsis", "tb", "syno", 
                 # "lcanc",
                 "airspace")
  metals.treat <- lapply(datsource, function(x) metals[[x]] %>% mutate(neutgroup = x) %>%
                           filter(source == "treatment"))
  names(metals.treat) <- datsource
  metadf <- bind_rows(metals.treat) %>%
    dplyr::select(neutgroup, sample) %>%
    rename("pointname" = "sample")
  
  countls.red <- lapply(datsource, function(x) countls[[x]] %>% 
                          dplyr::select("Gene", metals.treat[[x]]$sample))
  
  countdf <- reduce(countls.red, inner_join) %>%
    column_to_rownames(var = "Gene") %>%
    t() %>%
    as.data.frame %>%
    rownames_to_column(var = "sample")

  # countdf[is.na(countdf)] <- 0
  counts.novar <- lapply(countdf[,-1], function(x) var(x))
  counts.novar <- names(subset(counts.novar,counts.novar == 0))
  
  counts.nmds.in <- countdf %>%
    dplyr::select(-all_of(counts.novar))%>%
    column_to_rownames(var = "sample")

  set.seed(112358)
  counts.nmds <- metaMDS(counts.nmds.in)
  
  counts.nmds.points <- counts.nmds$points %>%
    as.data.frame() %>%
    rownames_to_column(var = "pointname") %>%
    mutate(pointtype = "sample")
  
  counts.nmds.species <- counts.nmds$species %>%
    as.data.frame() %>%
    rownames_to_column(var = "pointname") %>%
    mutate(pointtype = "gene")
  
  
  counts.nmds.ggdat <- bind_rows(counts.nmds.points, counts.nmds.species) %>%
    left_join(metadf)

  
  return(counts.nmds.ggdat)
}

