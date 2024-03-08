TCGA_box <- function(TCGA_CUS,TCGA_LM,cancer){
  if(is.null(cancer)) {cancer == ".*"}
  TCGA_dat <- TCGA_LM %>%
    left_join(TCGA_CUS)%>%
    rename("VAT Neutrophil" = "V.neutrophil", "Blood Neutrophil" = "Neutrophils", "C.B" = "Blood.neutrophil", "Tissue" = "Cancer") %>%
    mutate(file_name.expression = paste0(Mixture,".FPKM.txt.gz"))
  dat_plot <- TCGA_dat %>% 
    gather("cell_type", "relative_abundance", 3:5) %>%
    mutate(Tissue2nd = paste0(Tissue, cell_type)) %>%
    filter(!cell_type == "C.B")%>%
    filter(Tissue == ifelse(is.null(cancer), ".*",cancer))%>%
    mutate(patnum = as.numeric(as.factor(Tissue)))
  
  
  # categorize based on blood and VAT
  tiskey <- TCGA_dat %>% 
    group_by(Tissue) %>% 
    dplyr::summarise(sum_blood = sum(`Blood Neutrophil`), sum_vat = sum(`VAT Neutrophil`)) %>% 
    mutate(type = ifelse(sum_blood > sum_vat, "Blood Neutrophil", "VAT Neutrophil"))
  colnames(tiskey) <- c("Tissue", "sum_VAT","sum_blood","type")
  
  
  # create object to fill the plot background based on category
  numtissuekey <- dat_plot %>%
    select(Tissue, patnum) %>%
    distinct()%>%
    left_join(tiskey)%>%
    mutate(ylow = patnum - .45,
           yhigh = patnum + .45,
           xmin = 0,
           xmax = 0.75) %>%
    gather(xmin, xmax, key = "end", value = "x")
  
  list.output <- list(TCGA_dat,dat_plot, numtissuekey)
  names(list.output) <- c("TCGA_dat", "dat_plot", "numtissuekey")
  return(list.output)
}

## change the naming 