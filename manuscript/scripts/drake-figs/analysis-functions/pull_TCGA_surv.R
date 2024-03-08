pull_TCGA_surv <- function(){
  # BRCA_CLI <- read.csv("../../data/survival info/BRCA_CLI.csv", header = TRUE)
  # LUAD_CLI <- read.csv("../../data/survival info/LUAD_CLI.csv", header = TRUE)
  # LUSC_CLI <- read.csv("../../data/survival info/LUSC_CLI.csv", header = TRUE)
  # SARC_CLI <- read.csv("../../data/survival info/SARC_CLI.csv", header = TRUE)
  # COAD_CLI <- read.csv("../../data/survival info/COAD_CLI.csv", header = TRUE)
  # READ_CLI <- read.csv("../../data/survival info/READ_CLI.csv", header = TRUE)
  # KIRC_CLI <- read.csv("../../data/survival info/KIRC_CLI.csv", header = TRUE)
  # ALL_CLI <- bind_rows(LUAD_CLI,LUSC_CLI,READ_CLI,COAD_CLI,KIRC_CLI,SARC_CLI,BRCA_CLI)%>%
  #   select(file_id.BAM, days_to_last_follow_up,days_to_death,vital_status)%>%
  #   mutate(vital_status = ifelse(vital_status=="Alive",0,1),
  #          days = ifelse(is.na(days_to_death),days_to_last_follow_up,days_to_death))
  # 
  # ## combine manifest files
  # BRCA_MAN <- read.csv("../../data/survival info/BRCA_MAN.csv", header = TRUE)
  # LUAD_MAN <- read.csv("../../data/survival info/LUAD_MAN.csv", header = TRUE)
  # LUSC_MAN <- read.csv("../../data/survival info/LUSC_MAN.csv", header = TRUE)
  # SARC_MAN <- read.csv("../../data/survival info/SARC_MAN.csv", header = TRUE)
  # COAD_MAN <- read.csv("../../data/survival info/COAD_MAN.csv", header = TRUE)
  # READ_MAN <- read.csv("../../data/survival info/READ_MAN.csv", header = TRUE)
  # KIRC_MAN <- read.csv("../../data/survival info/KIRC_MAN.csv", header = TRUE)
  # ALL_MAN <- bind_rows(LUAD_MAN,LUSC_MAN,COAD_MAN,KIRC_MAN,SARC_MAN,READ_MAN,BRCA_MAN) %>%
  #   select(file_name.expression,file_id.BAM,file_id.expression)%>%
  #   rename("Mixture" = "file_id.expression")
  
  files_CLI <- list.files(path = "../../data/survival info", pattern = "*_CLI.csv", full.names = T)
  ALL_CLI <- sapply(files_CLI, read_csv, simplify=FALSE) %>%
    bind_rows(.id = "id")%>%
    select(file_id.BAM, days_to_last_follow_up,days_to_death,vital_status)%>%
    mutate(vital_status = ifelse(vital_status=="Alive",0,1),
           days = ifelse(is.na(days_to_death),days_to_last_follow_up,days_to_death))
  
  files_MAN <- list.files(path = "../../data/survival info", pattern = "*_MAN.csv", full.names = T)
  ALL_MAN <- sapply(files_MAN, read_csv, simplify=FALSE) %>%
    bind_rows(.id = "id")%>%
    select(file_name.expression,file_id.BAM,file_id.expression)%>%
    rename("Mixture" = "file_id.expression")
  
  list.output <- list(ALL_CLI, ALL_MAN)
  names(list.output) <- c("ALL_CLI","ALL_MAN")
  return(list.output)
}
