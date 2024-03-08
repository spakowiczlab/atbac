sigscoreSCRNA <- function(){
  scrna <- readRDS("../../data/scrna-neutrophils.rds")
  sig.genes <- read.delim("../../data/CIBERSORT/sig-genes/custom-sig-genes_bldat.txt")
  
  neut.discrim <- sig.genes %>%
    mutate(V.enriched = V.neutrophil - Blood.neutrophil) %>%
    arrange(desc(V.enriched))
  
  VAT.genes <- neut.discrim$NAME[1:10]
  
  scrna.frmt <- scrna %>%
    select("Sample", all_of(VAT.genes)) %>%
    column_to_rownames(var = "Sample") %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(var = "Gene.Symbol")
  
  scrna.sig.score <- calculateAvgZScore(scrna.frmt, VAT.genes) %>%
    rename("sig_score" = "avg_z_score", "Sample" = "sample") %>%
    left_join(scrna[c('Sample', 'Cell Type', 'umap_1', 'umap_2')], by = "Sample")
  
  return(scrna.sig.score)
}