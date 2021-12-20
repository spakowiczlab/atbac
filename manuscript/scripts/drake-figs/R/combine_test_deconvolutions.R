combine_test_deconvolutions <- function(){
  Blood1 <- read.csv("../../data/CIBERSORT/deconvolved-samples/Blood_GSE103889.csv")
  Blood2 <- read.csv("../../data/CIBERSORT/deconvolved-samples/Blood_GSE28064-GPL8432.csv")
  Blood3 <- read.csv("../../data/CIBERSORT/deconvolved-samples/Blood_GSE14771-GPL96.csv")
  Blood4 <- read.csv("../../data/CIBERSORT/deconvolved-samples/validate-blood.csv")
  VAT <- read.csv("../../data/CIBERSORT/deconvolved-samples/VAT_GDS4276.csv")
  VAT2 <- read.csv("../../data/CIBERSORT/deconvolved-samples/validate-at.csv")
  
  descells <- c("Blood.neutrophil", "V.neutrophil") #Combines and form vector of blood and vat cell types
  
  blood.frc1 <- Blood1 %>% 
    select(-P.value, -Pearson.Correlation, -RMSE) %>% #Leave out the P.calue Rearson.Correlation and RMSE
    gather(-Input.Sample, key = "celltype", value = "fraction") %>% #gather the colomns and collapse them into Celltypes and fraction
    mutate(celltype = ifelse(celltype %in% descells, celltype, "Other")) %>% #Put in Other if it doesn't fit the celltype
    group_by(Input.Sample, celltype) %>%
    summarise(fraction = sum(fraction)) %>%
    spread(key = "celltype", value = "fraction") %>%
    mutate(source = "blood")
  
  blood.frc2 <- Blood2 %>%
    select(-P.value, -Pearson.Correlation, -RMSE) %>%
    gather(-Input.Sample, key = "celltype", value = "fraction") %>%
    mutate(celltype = ifelse(celltype %in% descells, celltype, "Other")) %>%
    group_by(Input.Sample, celltype) %>%
    summarise(fraction = sum(fraction)) %>%
    spread(key = "celltype", value = "fraction") %>%
    mutate(source = "blood")
  
  
  blood.frc3 <- Blood3 %>%
    select(-P.value, -Pearson.Correlation, -RMSE) %>%
    gather(-Input.Sample, key = "celltype", value = "fraction") %>%
    mutate(celltype = ifelse(celltype %in% descells, celltype, "Other")) %>%
    group_by(Input.Sample, celltype) %>%
    summarise(fraction = sum(fraction)) %>%
    spread(key = "celltype", value = "fraction") %>%
    mutate(source = "blood")
  
  blood.frc4 <- Blood4 %>%
    select(-P.value, -Pearson.Correlation, -RMSE) %>%
    gather(-Input.Sample, key = "celltype", value = "fraction") %>%
    mutate(celltype = ifelse(celltype %in% descells, celltype, "Other")) %>%
    group_by(Input.Sample, celltype) %>%
    summarise(fraction = sum(fraction)) %>%
    spread(key = "celltype", value = "fraction") %>%
    mutate(source = "blood")
  
  at.frc <- VAT %>%
    select(-P.value, -Pearson.Correlation, -RMSE) %>%
    gather(-Input.Sample, key = "celltype", value = "fraction") %>%
    mutate(celltype = ifelse(celltype %in% descells, celltype, "Other")) %>%
    group_by(Input.Sample, celltype) %>%
    summarise(fraction = sum(fraction)) %>%
    spread(key = "celltype", value = "fraction") %>%
    mutate(source = "VAT")
  
  at.frc2 <- VAT2  %>%
    select(-P.value, -Pearson.Correlation, -RMSE, Blood.neutrophil) %>%
    gather(-Input.Sample, key = "celltype", value = "fraction") %>%
    mutate(celltype = ifelse(celltype %in% descells, celltype, "Other")) %>%
    group_by(Input.Sample, celltype) %>%
    summarise(fraction = sum(fraction)) %>%
    spread(key = "celltype", value = "fraction") %>%
    mutate(source = "VAT")
  
  
  frc <- bind_rows(blood.frc1,blood.frc2, blood.frc3, blood.frc4, at.frc, at.frc2) %>%
    arrange(desc(Input.Sample))
  
  all.frc <- as.data.frame(frc)
  all.frc[is.na(all.frc)] <- 0
  
  bardat <- all.frc %>%
    gather(-Input.Sample, -source, key = "celltype", value = "predfrac") 
  
  return(bardat)
}
