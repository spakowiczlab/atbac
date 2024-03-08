pull_GTEx <- function(){
  files_CUS <- list.files(path = "../../data/CIBERSORT/CIBERSORTx-deconvolved/GTEx", pattern = "*_CUS.csv", full.names = T)
  ciber_CUS <- sapply(files_CUS, read_csv, simplify=FALSE) %>% bind_rows(.id = "id")
  
  VAT <- ciber_CUS %>%
    mutate(Mixture = gsub("(GTEX)\\.(.*)\\.(\\d+)\\.(SM)\\.(.*)", "\\1\\-\\2\\-\\3\\-\\4\\-\\5",ciber_CUS$Mixture))%>%
    mutate(Sample = gsub("(GTEX)\\.(.*)\\.(\\d+)\\.(.*)\\.(SM)\\.(.*)", "\\1\\-\\2\\-\\3\\-\\4\\-\\5\\-\\6",Mixture)) %>%
    select(Sample,`V neutrophil`)%>%
    rename("VAT Neutrophil"="V neutrophil")
  
  
  files_LM <- list.files(path = "../../data/CIBERSORT/CIBERSORTx-deconvolved/GTEx", pattern = "*_LM22.csv", full.names = T)
  ciber_LM <- sapply(files_LM, read_csv, simplify=FALSE) %>% bind_rows(.id = "id")
  Blood <- ciber_LM %>% 
    mutate(Mixture = gsub("(GTEX)\\.(.*)\\.(\\d+)\\.(.*)\\.(SM)\\.(.*)", "\\1\\-\\2\\-\\3\\-\\4\\-\\5\\-\\6",ciber_LM$Mixture))%>%
    mutate(Sample = gsub("(GTEX)\\.(.*)\\.(\\d+)\\.(SM)\\.(.*)", "\\1\\-\\2\\-\\3\\-\\4\\-\\5",Mixture)) %>%
    select(Sample,`Neutrophils`) %>%
    rename("Blood Neutrophil"="Neutrophils")
  
  
  GTEx_dat <- left_join(VAT,Blood)
  
  return(GTEx_dat)
}
