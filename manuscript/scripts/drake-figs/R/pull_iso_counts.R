pull_iso_counts <- function(metals){
  bldvat <- read.table("../../data/ncbi_processed/iso_bldvat.txt", 
                       sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  sepsis <- read.table("../../data/ncbi_processed/iso_sepsis.txt",
                       sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  tb <- read.table("../../data/ncbi_processed/iso_tb.txt",
                   sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  exercise <- read.table("../../data/ncbi_processed/iso_exercise.txt",
                         sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  endotox <- read.table("../../data/ncbi_processed/iso_endotox.txt",
                        sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  syno <- read.table("../../data/ncbi_processed/iso_syno_GSE116899.txt",
                     sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  # lcanc <- read.table("../../data/ncbi_processed/iso_lung-canc_GSE68795.txt", 
  #                     sep = "\t", header = T, stringsAsFactors = F)
  # For this to work you're going to need to mess with the packages. But we are excluding lcanc anyway.
  # lcanc <- lcanc %>%
  #   column_to_rownames(var = "Entrez")
  # lcanc <- collapseRows(lcanc[,-1], rowGroup = lcanc$Symbol, rowID = rownames(lcanc))
  # lcanc <- lcanc$datETcollapsed %>%
  #   as.data.frame() %>%
  #   dplyr::select(lcanc.meta$sample)
  
  airspace <- endotox %>%
    gather(-Gene, key = "sample", value = "count") %>%
    mutate(count = round(count)) %>%
    spread(key = "sample", value = "count")
  
  expls <- list(bldvat, endotox, exercise, sepsis, tb, syno,
                # lcanc,
                airspace)
  datsource <- c("bldvat", "endotox", "exercise", "sepsis", "tb", "syno", 
                 # "lcanc",
                 "airspace")
  names(expls) <- datsource
  
  expls.rightsamps <- lapply(datsource, function(x) expls[[x]] %>% select(c("Gene", metals[[x]]$sample)))
  names(expls.rightsamps) <- datsource
  return(expls.rightsamps)
  
}
