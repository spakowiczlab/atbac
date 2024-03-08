mouseDistances <- function(){
  seqtab <- read.csv("../../data/mouse-data/second-experiment/ASV_Abundance_Table.csv")
  metadata <- seqtab %>%
    rename("sample" = "seqs") %>%
    select(sample) %>% 
    # Pre-format the combined infor
    mutate(sample.form = gsub("Week\\.", "Week", 
                              gsub("VATO", "VAT.O", 
                                   gsub("\\.$","",
                                        gsub("\\.\\.",".", 
                                             gsub("long\\.","",
                                                  gsub("Ly6G\\.","",
                                                       gsub("End\\.","",
                                                            gsub("Pellet\\.","",
                                                                 gsub("T2..T6821","T6821(T2)",sample))))))))),
           sample.form = ifelse(grepl("Pellet", sample), paste0("Pellet.", sample.form), sample.form)
    ) %>%
    separate(sample.form,remove = F,
             into = c("tissue", "bmi.status", "space", "timepoint", "diet", "mouse", "patient", "patient.sex"), sep = "\\.") %>%
    # Fix pellet issues
    mutate(patient.sex = ifelse(tissue == "Pellet", patient, patient.sex),
           patient = ifelse(tissue == "Pellet",mouse, patient),
           mouse = ifelse(tissue == "Pellet", NA, mouse)) %>%
    # Fix obese, lean gavage
    mutate(gavage = grepl("slurry", sample),
           bmi.status = ifelse(gavage == T, tissue, bmi.status),
           tissue = ifelse(gavage == F, tissue, NA),
           patient.sex = ifelse(gavage == T, diet, patient.sex),
           patient = ifelse(gavage ==T, timepoint, patient),
           diet = ifelse(gavage == F, diet, NA),
           timepoint = ifelse(gavage == F, timepoint, NA)) %>%
    # Fix saline gavages
    mutate(gavage = gavage == T | tissue == "Saline",
           bmi.status = gsub("\\d", "Saline", bmi.status),
           tissue = ifelse(tissue == "Saline", NA, tissue)) %>%
    # Misc indicators
    mutate(Ly6G = grepl("Ly6G", sample),
           timepoint = gsub("Week1s", "Week1", timepoint),
           long.exp = grepl("long", sample)) %>%
    # Fix patient IDs with infor from Dharti
    mutate(patient = ifelse(long.exp == T & bmi.status == "Obese", "T6821", patient),
           patient= ifelse(long.exp == F & bmi.status == "Obese" & tissue != "Pellet" & is.na(patient), "T6824", patient))
  
  # Calculate the distances
  bc.dist <- seqtab %>%
    column_to_rownames(var = "seqs") %>%
    vegan::vegdist(method = "bray") 
  
  bc.form <- bc.dist %>%
    as.matrix() %>%
    as.data.frame() %>%
    rownames_to_column(var = "sample") %>%
    pivot_longer(-sample, names_to = "sample2", values_to = "distance")
  
  compsamps <- metadata %>%
    filter(timepoint == "Week1" & tissue %in% c("VAT", "Lung", "Pellet") & bmi.status == "Obese") %>%
    select(sample,tissue)
  
  VAT.selected <- compsamps %>%
    filter(tissue == "VAT")
  other.select <- compsamps %>%
    filter(tissue != "VAT") %>%
    rename("sample2" = "sample")
  
  bc.plot <- bc.form %>%
    filter(sample %in% VAT.selected$sample) %>%
    inner_join(other.select)
  
  # Check species overlap
  specseqs <- seqtab %>%
    pivot_longer(-seqs, names_to = "sequence", values_to = "counts") %>%
    rename("sample" = "seqs") %>%
    filter(counts > 0) %>%
    inner_join(compsamps)
  
  VATseqs <- filter(specseqs, tissue == "VAT")$sequence
  lungseqs <- filter(specseqs, tissue == "Lung")$sequence
  stoolseqs <- filter(specseqs, tissue == "Pellet")$sequence
  
  listseqs <- list(VAT = unique(VATseqs),
                   Lung = unique(lungseqs),
                   Stool = unique(stoolseqs))
  
  euldat <- eulerr::euler(combinations = listseqs, shape = "ellipse")
  
  # Use this to get numbers in each group
  # compsamps %>%
  #   group_by(tissue) %>%
  #   tally()
  
  output <- list(bc.plot, euldat)
  names(output) <- c("bc.plot", "euldat")
  return(output)
}