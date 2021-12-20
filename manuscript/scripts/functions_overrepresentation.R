# Rebecca Hoyd's function to avoid clicking buttons in the gene ontology website
# This document provides code to build a database from HGNC symbols to PANTHER pathway mappings and perform an overrepresentation test using it.


# The code chunk commented out here displays how the db was produced.
# library(PANTHER.db)
# library(org.Hs.eg.db)
# 
# pthOrganisms(PANTHER.db) <- "HUMAN"
# 
# entrez.keys <- keys(PANTHER.db, keytype = "ENTREZ")
# hugo.trans <- AnnotationDbi::select(org.Hs.eg.db, keys = entrez.keys,
#                                     keytype = "ENTREZID", columns = c("SYMBOL", "ENTREZID")) %>%
#   filter(!is.na(SYMBOL))
# 
# 
# grab.mappings <- AnnotationDbi::select(PANTHER.db, keys = entrez.keys,
#                                        keytype = "ENTREZ", columns = c("ENTREZ", "PATHWAY_TERM")) %>%
#   rename("ENTREZID" = "ENTREZ") %>%
#   left_join(hugo.trans) %>%
#   filter(!is.na(SYMBOL))
# 
# all.syms <- unique(grab.mappings$SYMBOL)
# n.syms <- length(all.syms)
# 
# formatted.table <- grab.mappings %>%
#   add_count(PATHWAY_TERM) %>%
#   mutate(perc.genes = n/n.syms) %>%
#   dplyr::select(-n)
# # 
# write.table(formatted.table, "../data/pathway-db/overrep_custom-panther-paths.txt", sep = "\t",
#             quote = F,  row.names = F)

pathdb <- read.table("../data/pathway-db/overrep_custom-panther-paths.txt", sep = "\t", header = T, stringsAsFactors = F)

library(dplyr)
pathOverrepTest <- function(intgenes){
  ngenes = length(intgenes) 
  tmp.true <- as.data.frame(cbind(SYMBOL = intgenes, present = 1))
  
  calcdat <- pathdb %>%
    left_join(tmp.true) %>%
    filter(present == 1) %>%
    add_count(PATHWAY_TERM, name = "n.int") %>%
    group_by(PATHWAY_TERM) %>%
    summarise(exp.percent = unique(perc.genes),
              n.int = unique(n.int),
              genes.present = paste(SYMBOL, collapse = ", "))
    # mutate(pval = binom.test(x = n.int, n = ngenes, p = exp.percent)$p.value,
    #        padj = p.adjust(pval, method = "fdr"))
    
pvals <- unlist(lapply(1:nrow(calcdat), function(y) binom.test(x = calcdat$n.int[y], n = ngenes, p = calcdat$exp.percent[y])$p.value))

calcdat$obs.percent <- calcdat$n.int/ngenes
calcdat$pval <- pvals
calcdat$padj <- p.adjust(calcdat$pval, method = "fdr")

return(calcdat)
}
