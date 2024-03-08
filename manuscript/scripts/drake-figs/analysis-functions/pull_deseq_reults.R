pull_deseq_results <- function(){
  bldvat <- read.csv("../../data/DESeq2/bldvat.csv")
  endo <- read.csv("../../data/DESeq2/endotox.csv")
  exer <- read.csv("../../data/DESeq2/exercise.csv")
  seps <- read.csv("../../data/DESeq2/sepsis.csv")
  tb <- read.csv("../../data/DESeq2/tb.csv")
  syno <- read.csv("../../data/DESeq2/syno.csv")
  lcanc <- read.csv("../../data/DESeq2/lcanc.csv")
  airspace <- read.csv("../../data/DESeq2/airspace.csv")
  
  des.res <- list(bldvat, endo, exer, seps, tb, syno, lcanc, airspace)
  names(des.res) <- c("VAT", "Endotoxin", "Exercise", "Sepsis", "Active.TB", "Synovial", "Lung.Cancer",
                      "Airspace")
  return(des.res)
}
