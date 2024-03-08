correlate_neut_types <- function(desres){
  cor.df <- lapply(names(desres), function(x) desres[[x]] %>%
                     mutate(source = x) %>% 
                     dplyr::select(Gene, source, log2FoldChange)) %>%
    bind_rows() %>%
    spread(key = "source", value = "log2FoldChange") %>%
    drop_na()
  
  neutgroups <- colnames(cor.df)[-1]
  
  inres <- list()
  outres <- list()
  
  for(i in neutgroups){
    for(j in neutgroups){
      inres[[j]] <- cor.test(cor.df[[i]], cor.df[[j]], method = "spearman") %>%
        tidy() %>%
        dplyr::select(estimate, p.value) %>%
        mutate(ingroup = j,
               outgroup = i,
               pairing = paste(sort(c(i,j)), collapse = ","))
    }
    outres[[i]] <- inres
    inres <- list()
  }
  
  cor.res <- bind_rows(lapply(outres, function(x) bind_rows(x)))
  
  cor.dat <- cor.res %>%
    dplyr::select("estimate", "p.value", "pairing") %>%
    separate("pairing", into = c("neutrgroup1", "neutgroup2"), sep = ",")
  cor.dat <- cor.dat[!duplicated(cor.dat),]
}
