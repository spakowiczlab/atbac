pull_iso_metas <- function(){
  bldvat.meta <- readRDS("../../data/ncbi_processed/meta_bldvat.RDS")
  sepsis.meta <- readRDS("../../data/ncbi_processed/meta_sepsis.RDS")
  tb.meta <- readRDS("../../data/ncbi_processed/meta_tb.RDS")
  exercise.meta <- readRDS("../../data/ncbi_processed/meta_exercise.RDS")
  endotox.meta <- readRDS("../../data/ncbi_processed/meta_endotox.RDS")
  syno.meta <- readRDS("../../data/ncbi_processed/meta_syno.RDS")
  lcanc.meta <- read.csv("../../data/ncbi_processed/meta_lung-canc.csv", stringsAsFactors = F)
  
  bldvat.meta <- bldvat.meta %>%
    mutate(sample = as.character(sample),
           source = ifelse(source == "VAT", "treatment", "control")) %>%
    arrange(sample)
  sepsis.meta <- sepsis.meta %>%
    mutate(source = ifelse(`sample_type:ch1` == "Patient", "treatment", "control")) %>%
    arrange(geo_accession) %>%
    mutate(sample = geo_accession)
  tb.meta <- tb.meta %>%
    filter(grepl("Neut", title)) %>%
    mutate(source = ifelse(grepl("Control", characteristics_ch1.3), "control", "treatment")) %>% 
    arrange(geo_accession) %>%
    mutate(sample = geo_accession)
  airspace.meta <- endotox.meta %>%
    filter(cell.type != "in vitro") %>%
    mutate(source = ifelse(cell.type == "circulating", "control", "treatment"),
           sample = as.character(sample)) %>%
    arrange(sample)
  endotox.meta <- endotox.meta %>%
    filter(cell.type == "circulating") %>%
    mutate(source = ifelse(agent == "endotoxin", "treatment", "control"),
           sample = as.character(sample)) %>%
    arrange(sample)
  exercise.meta <- exercise.meta %>%
    mutate(source = ifelse(protocol == "after exercise", "treatment", "control"),
           sample = as.character(sample)) %>%
    arrange(sample)
  syno.meta <- syno.meta %>%
    mutate(sample = as.character(sampnames),
           source = ifelse(tissue == "synovial", "treatment", "control")) %>%
    arrange(sample)
  # lcanc.meta <- lcanc.meta %>%
  #   mutate(sample = as.character(samps),
  #          source = ifelse(status == "tumor", "treatment", "control")) %>%
  #   arrange(sample)
  
  metals <- list(bldvat.meta, endotox.meta, exercise.meta, sepsis.meta, tb.meta, syno.meta, 
                 # lcanc.meta,
                 airspace.meta)
  names(metals) <- c("bldvat", "endotox", "exercise", "sepsis", "tb", "syno",
                     # "lcanc", 
                     "airspace")
  return(metals)
}
