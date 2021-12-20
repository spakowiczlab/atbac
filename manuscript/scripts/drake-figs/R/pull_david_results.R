pull_david_results <- function(){
  davgroups <- c("endotox", "exercise", "sepsis", "tb", "bldvat", "syno", 
                 # "lcanc",
                 "airspace")
  davgroups <- lapply(davgroups, function(x) read.table(paste0("../../data/DAVID/output/fig2_top200_", x, ".txt"), 
                                                        sep = "\t", header = T, stringsAsFactors = F))
  names(davgroups) <- c("Endotoxin", "Exercise", "Sepsis", "Active.TB", "VAT", "Synovial", 
                        # "Lung.Cancer", 
                        "Airspace")
  
  return(davgroups)
}