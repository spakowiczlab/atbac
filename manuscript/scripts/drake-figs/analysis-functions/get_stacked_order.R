get_stacked_order <- function(taxcounts, bmi){
  sampord <- bmi %>%
    mutate(sample = paste0(sequence.id, "-I")) %>%
    left_join(taxcounts) %>%
    group_by(sequence.id) %>%
    mutate(totco = sum(count)) %>%
    ungroup() %>%
    filter(Phylum == "Firmicutes") %>%
    mutate(phylra = count/totco) %>%
    arrange(BMI, phylra) 
  
  return(sampord)
}