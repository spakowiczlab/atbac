plan <- drake_plan(
  reactkey = read.delim("../../data/pathway-db/ReactomePathways.txt", sep = "\t", header = F, stringsAsFactors = F) %>%
    rename("Term" = "V1",
           "Description" = "V2",
           "Species" = "V3") %>%
    select(-Species),
  des.res = pull_deseq_results(),
  david.res = pull_david_results(),
  neutmetas = pull_iso_metas(),
  neutcounts = pull_iso_counts(neutmetas),
  sourceres = read.table("../../data/atbac/mixing_proportions.txt",
                          sep = "\t", header = T),
  bmi.stat = read.csv("../../data/atbac/samples_obesity-status.csv", stringsAsFactors = F) %>%
    select(-sample.name) %>%
    mutate(label = paste0(sequence.id, "_int")) %>%
    filter(!is.na(BMI)),
  
  bacteria.desres = read.csv("../../data/atbac/DESeq2_all-levs_DNA.csv", stringsAsFactors = F),
  labelled.genes = read.csv("../../data/from-alecia/pathgenes.csv", stringsAsFactors = F),
  exprs = read.table("../../data/ncbi_processed/iso_bldvat.txt", sep = "\t", header = TRUE, stringsAsFactors = F,
                      check.names = F),
  exprs.meta = readRDS("../../data/ncbi_processed/meta_bldvat.RDS") %>%
    mutate(combcat = paste(source, obes.stat, sep = ".")),
  
  samples.bmi.firmicutes = get_stacked_order(taxcounts.nocontam.phyl, bmi.stat),
  taxcounts.nocontam.phyl = readRDS("../../data/atbac/taxcounts_interface-DNA_no-neg-contams.RDS") %>%
    mutate(Phylum = as.character(Phylum),
           Phylum = ifelse(is.na(Phylum), "Unclassified", Phylum)) %>%
    group_by(sample, Phylum) %>%
    summarise(count = sum(counts)), 
  
  siggenes = read.table("../../data/CIBERSORT/sig-genes/custom-sig-genes_bldat.txt", sep = "\t", header = TRUE, stringsAsFactors = F),
  siggenes.heatdend = format_sig_heat(siggenes),
  siggenes.long = siggenes.heatdend[[1]],
  siggenes.celldend = siggenes.heatdend[[2]],
  siggenes.enriched = choose_enriched_genes(siggenes),
  
  GenePath.meta = format_pathway_dat(exprs, exprs.meta, labelled.genes),
  pca.bldvat.exprs = pull_pca_res(exprs, exprs.meta),
  treatment.genepools = cluster_treatment_genepools(read.csv("../../data/DESeq2/compare-genes_treat.csv")),
  correlate.neut.sigs = correlate_neut_types(des.res),
  # nmds.neuts.pathways = nmds_neuts_pathways(david.res, reactkey), No longer using this version.
  nmds.neuts.counts = nmds_neuts_counts(neutcounts, neutmetas),
  deconvolved.fractions = combine_test_deconvolutions(),
  
  
  ############################
  TCGA_deconv = pull_TCGA_deconvolved(),
  TCGA_box_dat = TCGA_box(TCGA_deconv$TCGA_CUS, TCGA_deconv$TCGA_LM,"COAD"),
  TCGA_SURV_info = pull_TCGA_surv(),
  TCGA_surv_dat = TCGA_SURV(TCGA_box_dat$TCGA_dat, TCGA_SURV_info$ALL_CLI, TCGA_SURV_info$ALL_MAN,"COAD"),
  canc.df = TCGA_hoyd(TCGA_surv_dat$SURV_d,"COAD"),
  clinical = read.csv("../../data/clinical.csv",stringsAsFactors = F),
  ra_COAD = read.csv(file.path(paths$ra_COAD, "ra_COAD.csv"),stringsAsFactors = F),
  # collected from website (https://globalrph.com/bacteriacat/gram-negative-bacteria/)
  gram_neg_list = read.csv("../../data/gram negative list.csv")%>%
    mutate(gram.negative = trimws(gram.negative, which = c("both"))),
  neut = TCGA_micro_neut(clinical, ra_COAD,gram_neg_list,TCGA_box_dat$dat_plot),
  GTEx_Portal = read_csv("../../data/GTEX_v8_sample-location-key.csv"),
  GTEx_dat = pull_GTEx(),
  GTEx_box_dat = GTEx_box(GTEx_Portal,GTEx_dat)
)
