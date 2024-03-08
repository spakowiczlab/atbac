mouseVulcRange <- function(){
  load("../../data/plotdat_mouse-exp_vulcanrange.rda")
  mickey <- readRDS("../../data/microbes_found-in-tissues.rds") 
  colnames(mickey)[1] <- "term"
  rangedat.form2 <- rangedat.form %>%
    drop_na(direction) %>%
    left_join(mickey) %>%
    mutate(tissue = replace_na(tissue, "Other"))
  
  output <- list(rangedat.form, rangedat.form2, patkey)
  names(output) <- c("rangedat.form", "rangedat.form2", "patkey")
  return(output)
}