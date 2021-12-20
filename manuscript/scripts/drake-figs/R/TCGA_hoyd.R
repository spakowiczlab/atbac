TCGA_hoyd <- function(SURV_d,cancer){
  cancers <- unique(SURV_d$Tissue)
  percentiles<- seq(0.01,0.99,0.01)
  ## creating a dataframe for saving the p values ##
  ptable <- list()
  percents.list <- list()
  canc.list <- list()
  for(e in cancers){
    ##create subset for each cancer 
    canc_set <- SURV_d[SURV_d$Tissue== e,]
    for (p in percentiles) {  ## loop through every percentile
      #find the threshold based on the percentile
      cut <- quantile(canc_set$`VAT Neutrophil`, p)
      canc_set$type <- ifelse(canc_set$`VAT Neutrophil` > cut, 1, 0)
      survtemp <- coxph(Surv(days, vital_status) ~ type, canc_set)
      ##put wanted value in dataframe format
      dat_temp <- merge(data.frame(summary(survtemp)[["conf.int"]]), data.frame(summary(survtemp)[["coefficients"]]))
      ## extract value from dataframe
      ptable[[e]] <- dat_temp %>%
        mutate(hazard.ratio = exp..coef.,
               low.bound = lower..95,
               upper.bound = upper..95,
               pval = log(Pr...z..),
               Percentiles = p,
               Tissue = e) %>%
        select(hazard.ratio, low.bound, upper.bound, Percentiles, pval, Tissue)
      percents.list[[as.character(p)]] <- ptable
      ptable <- list()
    }
    canc.list[[e]] <- bind_rows(lapply(percents.list, function(x) bind_rows(x)))
    percents.list <- list()
  }
  ## decide whether the line goes up or down
  canc.df <- bind_rows(canc.list) %>%
    mutate(`hazard ratio` = ifelse(hazard.ratio >= 1, ">=1", "<1"))%>%
    filter(Tissue == ifelse(is.null(cancer), ".*",cancer))
  
  return(canc.df)
}