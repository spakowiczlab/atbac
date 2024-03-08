zymoFauxDes <- function(){
  zymo.seqtab <- read.csv("../../data/mouse-data/first-experiment/ASV_Abundance_Table.csv")
  zymo.tax <- read.table("../../data/mouse-data/first-experiment/ASV_tax_assignments.txt", sep = "\t")
  
  formatTax <- function(){
    tmp <- zymo.tax %>%
      select(V1, V2) %>%
      separate(V2, into = c("kingdom", "phylum", "class", "order",
                            "family", "genus", "species"),
               sep = ";") %>%
      mutate(species = gsub("s__", "", species),
             species = paste0("s__", gsub("g__", "", genus),".", species))
    tmp$species <- make.names(tmp$species)
    return(tmp)
    
  }
  
  exoRAtowide <- function(data, taxlev){
    tmp <- data
    tmp$Taxa <- tmp[[taxlev]]
    tmp.wide<- tmp %>%
      filter(!is.na(Taxa)) %>%
      group_by(seqs, Taxa) %>%
      summarise(ra = sum(counts, na.rm = T)) %>%
      spread(key = "Taxa", value = "ra")
    
    return(tmp.wide)
  }
  
  exoToDF <- function(taxalevels = c("kingdom", "phylum", "class", 
                                     "order", "family", "genus", "species"), 
                      data){
    w.ls <- lapply(taxalevels, function(x) exoRAtowide(data, x))
    w.df <- reduce(w.ls, function(x,y) left_join(x,y)) 
    return(w.df)
  }
  
  formatModelData <- function(){
    tax <- formatTax()
    dessamps <- zymo.seqtab %>%
      filter(grepl("HFD.VAT", seqs) & !grepl("Saline", seqs))
    all0 <- names(subset(colSums(dessamps[,-1]), colSums(dessamps[,-1]) == 0))
    
    desseqs <- dessamps %>%
      select(-all0) %>%
      pivot_longer(-seqs, names_to = "V1", values_to = "counts") %>%
      left_join(tax)
    
    desseqs.w <- exoToDF(data = desseqs)
    return(desseqs.w)
  }
  
  capture.models.univ <- function(outcome, lfun, modin, mics){
    mods.list <- lapply(mics, function(x) try({glm(as.formula(paste0(outcome, " ~ `", x, "`")), family = Gamma, data = modin) %>%
        broom::tidy()})
    )
    
    mods.list.clean <- rlist::list.clean(mods.list, function(x) is.null(x))
    mods.df <- bind_rows(mods.list.clean)
    return(mods.df)
  }
  
  calculatePValues <- function(){
    moddat <- formatModelData()
    microbes <- colnames(moddat[,-1])
    
    moddat.lab <- moddat %>%
      mutate(ob.lean = ifelse(grepl("Lean", seqs), 1,2))
    
    mod.pvals <- capture.models.univ("ob.lean", "gamma", moddat.lab, microbes) %>%
      filter(term != "(Intercept)") %>%
      select(term, p.value)
    
    return(mod.pvals)
    
  }
  
  calculateFoldChange <- function(){
    moddat <- formatModelData()
    
    foldcalc <- moddat %>%
      pivot_longer(-seqs, names_to = "microbes", values_to = "counts") %>%
      mutate(gavage = gsub("gavage.*", "", seqs)) %>%
      group_by(gavage, microbes) %>%
      summarise(gmean = mean(counts)) %>%
      pivot_wider(names_from = "gavage", values_from = gmean) %>%
      mutate(log2foldchange = log(Obese/Lean, base = 2)) %>%
      rename("term" = "microbes") %>%
      select(term, log2foldchange)
    
    return(foldcalc)
  }
  
  volcanoPlotDat <- function(){
    pvals <- calculatePValues()
    fold <- calculateFoldChange()
    
    plotin <- pvals %>%
      left_join(fold) %>%
      mutate(taxlev = ifelse(p.value < 0.1, str_sub(term, 1, 1), NA),
             textlab = ifelse(p.value < 0.1, term, NA),
             inf.rescale = ifelse(log2foldchange == Inf, 10,
                                  ifelse(log2foldchange == -Inf, -10, log2foldchange)))
    
    return(plotin)
  }
  
  volplot <- volcanoPlotDat()
  
  pval.table <- volplot %>%
    mutate(term = gsub("`", "", term)) %>%
    select(-taxlev, -textlab, -inf.rescale) %>%
    rename("microbe" = "term")
  
  miccounts <- formatModelData()
  micmeans <- miccounts %>%
    pivot_longer(-seqs, names_to = "microbe", values_to = "counts") %>%
    mutate(gavage = paste0("mean.", gsub("gavage.*", "", seqs))) %>%
    group_by(gavage, microbe) %>%
    summarise(gmean = mean(counts)) %>%
    pivot_wider(names_from = "gavage", values_from = gmean) %>%
    left_join(pval.table)
  
  output <- list(volplot, micmeans)
  names(output) <- c("volplot", "micmeans") 
  return(output)
}