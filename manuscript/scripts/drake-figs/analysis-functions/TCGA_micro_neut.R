TCGA_micro_neut <- function(clinical, ra_COAD,gram_neg_list, dat_plot){
  TCGA_id_link <- read_csv("../../data/TCGA-link-bam-and-expression-files.csv")
  g <- gram_neg_list$gram.negative
  id <- select(clinical, sample.microbe,sample.expression)%>%
    rename("Mixture"="sample.expression","sample"="sample.microbe")
  TCGA_id <- TCGA_id_link %>%
    rename("Mixture"="file_id.expression", "sample"="file_id.BAM")%>%
    rbind(id)
  
  
  micro <- ra_COAD %>%
    mutate(gg = gsub("g__(.*)(-s.*)","\\1",genus),
           gg = gsub("g__(.*)","\\1",gg))%>%
    filter(gg %in% g)%>%
    select(sample, exo.ra, microbe)%>%
    spread(microbe, exo.ra)%>%
    mutate(sum.ra = rowSums(.[,-1]))%>%
    select(sample, sum.ra)
  
  neut <- TCGA_id %>%
    left_join(dat_plot)%>%
    left_join(micro)%>%
    select(-Tissue2nd,-patnum,-file_name.expression) %>%
    filter(!cell_type == "C.B",
           !is.na(sum.ra)) %>%
    mutate(relative_abundance = relative_abundance,
           sum.ra = log(sum.ra))
  
  return(neut)
  
}