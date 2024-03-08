TCGA_SURV <- function(TCGA_dat,ALL_CLI,ALL_MAN,cancer){
  # SURV_CR <- TCGA_dat %>%
  #   filter(Tissue %in% c("COAD","READ"))%>%
  #   left_join(ALL_MAN, by = "file_name.expression")%>%
  #   left_join(ALL_CLI)%>%
  #   select(file_id.BAM,days,vital_status,Tissue,VAT_Neut)
  # SURV_REST <- TCGA_dat %>%
  #   filter(!Tissue %in% c("COAD","READ"))%>%
  #   left_join(ALL_MAN, by = "Mixture")%>%
  #   left_join(ALL_CLI)%>%
  #   select(file_id.BAM,days,vital_status,Tissue,VAT_Neut)
  # SURV_dat <- bind_rows(SURV_CR,SURV_REST)
  # 
  # ## plot using facet wrap
  # SURV_d <- SURV_dat %>%
  #   mutate(VAT_type = ifelse(VAT_Neut > sort(VAT_Neut)[0.9*length(VAT_Neut)],"More VAT","less VAT"))%>%
  #   ### filter to only COAD
  #   filter(Tissue == ifelse(is.null(cancer), ".*",cancer))
  # surv_plot <- survfit(Surv(days, vital_status) ~ VAT_type, data = SURV_d)
  
  
  SURV_CR <- TCGA_dat %>%
    filter(Tissue %in% c("COAD","READ"))%>%
    left_join(ALL_MAN, by = "file_name.expression")%>%
    left_join(ALL_CLI)%>%
    select(file_id.BAM,days,vital_status,Tissue,`VAT Neutrophil`,`Blood Neutrophil`)
  SURV_REST <- TCGA_dat %>%
    filter(!Tissue %in% c("COAD","READ"))%>%
    left_join(ALL_MAN, by = "Mixture")%>%
    left_join(ALL_CLI)%>%
    select(file_id.BAM,days,vital_status,Tissue,`VAT Neutrophil`,`Blood Neutrophil`)
  SURV_dat <- bind_rows(SURV_CR,SURV_REST)
  
  ## plot using facet wrap (using the threshold)
  SURV_d <- SURV_dat %>%
    ### filter to only COAD
    filter(Tissue == ifelse(is.null(cancer), ".*",cancer)) %>%
    mutate(VAT_type = ifelse(`VAT Neutrophil` > quantile(`VAT Neutrophil`, .9),"More VAT Neutrophil","Less VAT Neutrophil"),
           Blood_type = ifelse(`Blood Neutrophil` > sort(`Blood Neutrophil`)[0.9*length(`Blood Neutrophil`)],"More Blood Neutrophil","Less Blood Neutrophil"))
    
    
  
  surv_plot_VAT <- survfit(Surv(days, vital_status) ~ VAT_type, data = SURV_d)
  surv_plot_Blo <- survfit(Surv(days, vital_status) ~ Blood_type, data = SURV_d)
  list.output <- list(SURV_d, surv_plot_VAT, surv_plot_Blo)
  names(list.output) <- c("SURV_d","surv_plot_VAT", "surv_plot_Blo")
  return(list.output)
}