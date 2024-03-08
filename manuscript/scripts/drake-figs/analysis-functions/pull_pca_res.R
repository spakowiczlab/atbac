pull_pca_res <- function(exprs, meta){
  exprs.2 <- exprs %>%
    remove_rownames() %>%
    column_to_rownames(var = "Gene") %>%
    t() %>%
    as.data.frame()%>%
    rownames_to_column(var = "sample") %>%
    arrange(sample)
  
  meta <- exprs.2 %>%
    select(sample) %>%
    left_join(meta)
  
  selcol <- bind_cols(lapply(exprs.2[,-1], function(x) var(x))) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(var = "Gene") %>%
    filter(V1 > 0)
  
  exprs.4 <- exprs.2 %>%
    dplyr::select(selcol$Gene)
  
  x <- prcomp(exprs.4)
  
  datpoints <- x$x[, c(1,2)] %>%
    as.data.frame() %>%
    mutate(sample = meta$sample) %>%
    left_join(meta)
  
  output <- list(datpoints, x)
  names(output) <- c("plotdat", "pca")
  return(output)
}
