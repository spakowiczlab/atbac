GTEx_box <- function(GTEx_Portal,GTEx_dat){
  ncib <- GTEx_Portal %>%
    select(SMTS, SAMPID) %>%
    rename("Tissue"="SMTS", "Sample"="SAMPID") %>%
    left_join(GTEx_dat, )

  #order tissue by VAT
  GTEx_tiskey <- ncib %>%
    group_by(Tissue) %>%
    summarise(med_VAT = median(`VAT Neutrophil`), sum_VAT = sum(`VAT Neutrophil`),
              med_blood = median(`Blood Neutrophil`), sum_blood = sum(`Blood Neutrophil`))%>%
    arrange(med_VAT)


  #reorder the level of dat_plot$Tissue to change the y axis (change based on which order y axis should be)
  GTEx_plot <- ncib %>%
    tidyr::gather("cell_type", "relative_abundance", 3:4)
  GTEx_plot$Tissue <- factor(GTEx_plot$Tissue, levels = GTEx_tiskey$Tissue)
  GTEx_plot <- GTEx_plot [!duplicated(GTEx_plot ),] %>%
    mutate(patnum = as.numeric(as.factor(Tissue)))

#create object for background filling
  numpatkey <- GTEx_plot %>%
    select(Tissue, patnum)%>%
    distinct()%>%
    left_join(GTEx_tiskey) %>%
    mutate(type = ifelse(sum_VAT > sum_blood,"VAT Neutrophil","Blood Neutrophil"),
           ylow = patnum - .45,
           yhigh = patnum + .45,
           xmin = 0,
           xmax = 0.75) %>%
    tidyr::gather(xmin, xmax, key = "end", value = "x")

  list.output <- list(GTEx_plot,numpatkey)
  names(list.output) <- c("GTEx_plot","numpatkey")
  return(list.output)
}
