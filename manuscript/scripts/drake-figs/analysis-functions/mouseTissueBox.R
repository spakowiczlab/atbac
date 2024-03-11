mouseTissueBox <- function(){
  x <- read_csv("../../data/mouse-data/first-experiment/absolute.abundance.csv") %>%
    drop_na(sample_id)
  
  tissue <- 
    x %>%
    mutate(source = if_else(grepl("Pre|Post|End|^Leangavage$|^Obesegavage$", customer_label),
                            true = "stool",
                            false = "tissue"),
           tissue_type = if_else(grepl("VAT", customer_label),
                                 true = "VAT",
                                 false = if_else(grepl("Liver", customer_label),
                                                 true = "Liver",
                                                 false = if_else(grepl("Brain", customer_label),
                                                                 true = "Brain",
                                                                 false = "stool"))),
           Gavage = if_else(grepl("Lean", customer_label),
                            true = "Lean",
                            false = if_else(grepl("Obese", customer_label),
                                            true = "Obese",
                                            false = "Saline")),
           Diet = if_else(grepl("HFD", customer_label),
                          true = "HFD",
                          false = if_else(grepl("Chow", customer_label),
                                          true = "Chow",
                                          false = "Other")))
  
  return(tissue)
}
