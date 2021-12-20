pull_TCGA_deconvolved <- function(){
  #load files to match id
  TCGA_id_link <- read_csv("../../data/TCGA-link-bam-and-expression-files.csv")
  TCGA.COAD_id <- read.delim("../../data/gdc_manifest.2020-05-28.txt")
  TCGA.READ_id <- read.delim("../../data/gdc_manifest.2020-05-28 _READ.txt")
  ## add cancer type to the data frame
  TCGA.COAD_id <- data.frame(TCGA.COAD_id$id, rep("COAD", length(TCGA.COAD_id$id)))
  TCGA.READ_id <- data.frame(TCGA.READ_id$id, rep("READ", length(TCGA.READ_id$id)))
  colnames(TCGA.COAD_id) <- c("file_id.BAM", "Tissue")
  colnames(TCGA.READ_id) <- c("file_id.BAM", "Tissue")
  TCGA_tissue_id <- rbind(TCGA.COAD_id, TCGA.READ_id)
  # load cibersortx file
  TCGA_CUS <- read.csv("../../data/CIBERSORT/CIBERSORTx-deconvolved/COAD-READ_cus.csv", header = TRUE)
  TCGA_LM <- read.csv("../../data/CIBERSORT/CIBERSORTx-deconvolved/COAD-READ_lm.csv", header = TRUE)
  BRCA_C <- read.csv("../../data/CIBERSORT/CIBERSORTx-deconvolved/BRCA_CUS.csv", header = TRUE)%>%
    mutate(Cancer = "BRCA")
  SARC_C <- read.csv("../../data/CIBERSORT/CIBERSORTx-deconvolved/SARC_CUS.csv", header = TRUE)%>%
    mutate(Cancer = "SARC")
  LUSC_C <- read.csv("../../data/CIBERSORT/CIBERSORTx-deconvolved/LUSC_CUS.csv", header = TRUE)%>%
    mutate(Cancer = "LUSC")
  LUAD_C <- read.csv("../../data/CIBERSORT/CIBERSORTx-deconvolved/LUAD_CUS.csv", header = TRUE)%>%
    mutate(Cancer = "LUAD")
  BRCA_L <- read.csv("../../data/CIBERSORT/CIBERSORTx-deconvolved/BRCA_LM22.csv", header = TRUE)%>%
    mutate(Cancer = "BRCA")
  SARC_L <- read.csv("../../data/CIBERSORT/CIBERSORTx-deconvolved/SARC_LM22.csv", header = TRUE)%>%
    mutate(Cancer = "SARC")
  LUSC_L <- read.csv("../../data/CIBERSORT/CIBERSORTx-deconvolved/LUSC_LM22.csv", header = TRUE)%>%
    mutate(Cancer = "LUSC")
  LUAD_L <- read.csv("../../data/CIBERSORT/CIBERSORTx-deconvolved/LUAD_LM22.csv", header = TRUE)%>%
    mutate(Cancer = "LUAD")
  KIRC_C <- read.csv("../../data/CIBERSORT/CIBERSORTx-deconvolved/KIRC_CUS.csv", header = TRUE)%>%
    mutate(Cancer = "KIRC")
  KIRC_L <- read.csv("../../data/CIBERSORT/CIBERSORTx-deconvolved/KIRC_LM22.csv", header = TRUE)%>%
    mutate(Cancer = "KIRC")
  
  ##combine tissue and deconvolved CUS
  # TCGA_tissue_id <- left_join(TCGA_id_link, TCGA_tissue_id)
  # colnames(TCGA_tissue_id) <- c("file_id.BAM","Mixture","Cancer")
  # TCGA_CUS <- left_join(TCGA_CUS, TCGA_tissue_id)
  # TCGA_CUS <- TCGA_CUS %>% select(Cancer, everything())
  TCGA_tissue_id <- TCGA_tissue_id %>%
    left_join(TCGA_id_link) %>%
    rename("Cancer"="Tissue", "Mixture"="file_id.expression")
  
  
  ##combine tissue and deconvolved LM
  # TCGA_LM <- left_join(TCGA_LM, TCGA_tissue_id)
  # TCGA_LM <- TCGA_LM %>% select(Cancer, everything())
  # TCGA_CUS <- bind_rows(TCGA_CUS,SARC_C,LUAD_C,LUSC_C,KIRC_C,BRCA_C)%>%
  #   dplyr::select(Mixture, Cancer, V.neutrophil, Blood.neutrophil)
  # TCGA_LM <- bind_rows(TCGA_LM,SARC_L,LUAD_L,LUSC_L,KIRC_L, BRCA_L)%>%
  #   dplyr::select(Mixture, Cancer, Neutrophils)
  TCGA_CUS <- TCGA_CUS %>%
    left_join(TCGA_tissue_id)%>%
    select(Cancer, everything())
  TCGA_LM <- TCGA_LM %>%
    left_join(TCGA_tissue_id)%>%
    select(Cancer, everything())
  TCGA_CUS <- bind_rows(TCGA_CUS,SARC_C,LUAD_C,LUSC_C,KIRC_C,BRCA_C)%>%
    dplyr::select(Mixture, Cancer, V.neutrophil, Blood.neutrophil)
  TCGA_LM <- bind_rows(TCGA_LM,SARC_L,LUAD_L,LUSC_L,KIRC_L, BRCA_L)%>%
    dplyr::select(Mixture, Cancer, Neutrophils)
    # rename("file_id.expression" = "Mixture")
  
  #return(list(TCGA_CUS, TCGA_LM))
  list.output <- list(TCGA_CUS, TCGA_LM)
  names(list.output) <- c("TCGA_CUS","TCGA_LM")
  return(list.output)
}