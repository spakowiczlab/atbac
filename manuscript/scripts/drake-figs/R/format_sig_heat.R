format_sig_heat <- function(genemat){
  newlabs <- c("NAME", "Naive B Cells", "Memory B Cells", 
               "Plasma Cells", "CD8 T Cells", "Naive CD4 T Cells", 
               "Resting Memory CD4 T Cells", "Activated Memory CD4 T Cells", "Follicular Helper T Cells",
               "Gamma Delta T Cells", "Regulatory T Cells", "Resting NK Cells", 
               "Activated NK Cells", "Monocytes", "M0 Macrophages",
               "M1 Macrophages", "M2 Macrophages", "Resting Dendritic Cells", 
               "Activated Dendritic Cells", "Resting Mast Cells", "Activated Mast Cells", 
               "Eosonophils", "Peripheral Blood Neutrophils", "Visceral AT Neutrophils")
  
  colnames(genemat) <- newlabs
  
  genenames <- genemat$NAME
  
  genemat2 <- as.matrix(genemat[,-1])
  rownames(genemat2) <- genenames
  choosegenes <- rev(sort(rowSums(genemat2)))[1:100]
  choosegenes <- names(choosegenes)
  genemat2 <- t(genemat2) %>%
    as.data.frame() %>%
    select(choosegenes) %>%
    as.matrix()
  
  geneord <- genemat2 %>%
    t() %>%
    dist() %>%
    hclust() %>%
    dendro_data()
  geneord <- as.character(geneord$labels$label)
  
  dend.dat <- genemat2 %>%
    dist() %>%
    hclust() %>%
    dendro_data()
  cellord <- as.character(dend.dat$labels$label)
  
  heatdat <- genemat2 %>%
    as.data.frame() %>%
    rownames_to_column(var = "celltype") %>%
    gather(-celltype, key = "gene", value = "counts") %>% 
    mutate(gene = fct_relevel(gene, geneord),
           celltype = fct_relevel(celltype, cellord))
  
  return(list(heatdat, dend.dat))
}
