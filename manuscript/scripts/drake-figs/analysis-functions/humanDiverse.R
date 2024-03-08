humanDiverse <- function(){
  seqtab <- readRDS("../../../exploratory/data/seqtabNoCf_update.RDS")
  tax <- readRDS("../../../exploratory/data/taxf_update.RDS")
  
  bmi.stat <- read.csv("../../data/atbac/samples_obesity-status.csv", stringsAsFactors = F) %>%
    select(-sample.name) %>%
    mutate(label = paste0(sequence.id, "-I")) %>%
    filter(!is.na(BMI))
  
  spec.counts <- seqtab %>%
    as.data.frame() %>%
    rownames_to_column(var = "label") %>%
    pivot_longer(-label, names_to = "sequence", values_to = "counts") %>%
    filter(counts != 0) %>%
    group_by(label) %>%
    tally() %>%
    inner_join(bmi.stat)
  
  # Check normality, then significance of number of species in BMI groups
  # shapiro.test(spec.counts$n)
  # t.test(spec.counts$n ~ spec.counts$BMI)
  # kruskal.test(spec.counts$n ~ spec.counts$BMI)
  
  div.calc <- vegan::diversity(seqtab, index = "simpson") %>%
    as.data.frame() %>%
    rownames_to_column(var = "label") %>%
    inner_join(bmi.stat)
  
  # Check significance of diversity differences in BMI groups
  # kruskal.test(div.calc$. ~ div.calc$BMI)
  
  output <- list(spec.counts, div.calc)
  names(output) <- c("spec.counts", "div.calc")
  return(output)
}