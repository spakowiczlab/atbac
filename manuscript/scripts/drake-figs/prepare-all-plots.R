func.files <- list.files("analysis-functions/", full.names = T)
lapply(func.files, source)

# Analyses related to human data

reactkey = read.delim("../../data/pathway-db/ReactomePathways.txt", sep = "\t", header = F, stringsAsFactors = F) %>%
  rename("Term" = "V1",
         "Description" = "V2",
         "Species" = "V3") %>%
  select(-Species)
des.res = pull_deseq_results()
david.res = pull_david_results()
neutmetas = pull_iso_metas()
neutcounts = pull_iso_counts(neutmetas)
sourceres = read.table("../../data/atbac/mixing_proportions.txt",
                       sep = "\t", header = T)
bmi.stat = read.csv("../../data/atbac/samples_obesity-status.csv", stringsAsFactors = F) %>%
  select(-sample.name) %>%
  mutate(label = paste0(sequence.id, "_int")) %>%
  filter(!is.na(BMI))

bacteria.desres = read.csv("../../data/atbac/DESeq2_all-levs_DNA.csv", stringsAsFactors = F)
labelled.genes = read.csv("../../data/from-alecia/pathgenes.csv", stringsAsFactors = F)
exprs = read.table("../../data/ncbi_processed/iso_bldvat.txt", sep = "\t", header = TRUE, stringsAsFactors = F,
                   check.names = F)
exprs.meta = readRDS("../../data/ncbi_processed/meta_bldvat.RDS") %>%
  mutate(combcat = paste(source, obes.stat, sep = "."))

taxcounts.nocontam.phyl = readRDS("../../data/atbac/taxcounts_interface-DNA_no-neg-contams.RDS") %>%
  mutate(Phylum = as.character(Phylum),
         Phylum = ifelse(is.na(Phylum), "Unclassified", Phylum)) %>%
  group_by(sample, Phylum) %>%
  summarise(count = sum(counts))
samples.bmi.firmicutes = get_stacked_order(taxcounts.nocontam.phyl, bmi.stat)

human.div <- humanDiverse()
treatment.scores <- sigscoreActive()
scrna.sig.score <- sigscoreSCRNA()

siggenes = read.table("../../data/CIBERSORT/sig-genes/custom-sig-genes_bldat.txt", sep = "\t", header = TRUE, stringsAsFactors = F)
siggenes.heatdend = format_sig_heat(siggenes)
siggenes.long = siggenes.heatdend[[1]]
siggenes.celldend = siggenes.heatdend[[2]]
siggenes.enriched = choose_enriched_genes(siggenes)

GenePath.meta = format_pathway_dat(exprs, exprs.meta, labelled.genes)
pca.bldvat.exprs = pull_pca_res(exprs, exprs.meta)
treatment.genepools = cluster_treatment_genepools(read.csv("../../data/DESeq2/compare-genes_treat.csv"))
correlate.neut.sigs = correlate_neut_types(des.res)
# nmds.neuts.pathways = nmds_neuts_pathways(david.res, reactkey), No longer using this version.
nmds.neuts.counts = nmds_neuts_counts(neutcounts, neutmetas)
deconvolved.fractions = combine_test_deconvolutions()

save(samples.bmi.firmicutes, file = "prepared-data/sample-bmi-firmicutes.rda")
save(taxcounts.nocontam.phyl, file = "prepared-data/taxcounts-nocontam-phyl.rda")
save(sourceres, file = "prepared-data/sourceres.rda")
save(bacteria.desres, file = "prepared-data/bacteria-desres.rda")
save(GenePath.meta, file = "prepared-data/GenePath-meta.rda")
save(pca.bldvat.exprs, file = "prepared-data/pca-bldvat-exprs.rda")
save(treatment.genepools, file = "prepared-data/treatment-genepools.rda")
save(correlate.neut.sigs, file = "prepared-data/correlate-neut-sigs.rda")
save(nmds.neuts.counts, file = "prepared-data/nmds-neut-counts.rda")
save(des.res, file = "prepared-data/des-res.rda")
save(human.div, file = "prepared-data/human-div.rda")
save(treatment.scores, file = "prepared-data/score-activation.rda")
save(scrna.sig.score, file = "prepared-data/scrna-sig-score.rda")
save(siggenes.long, file = "prepared-data/siggenes-long.rda")
save(siggenes.celldend, file = "prepared-data/siggenes-celldend.rda")
save(siggenes.enriched, file = "prepared-data/siggenes-enriched.rda")
save(deconvolved.fractions, file = "prepared-data/deconvolved-fractions.rda")

############################ Analyses related to TCGA and GTex data
TCGA_deconv = pull_TCGA_deconvolved()
TCGA_box_dat = TCGA_box(TCGA_deconv$TCGA_CUS, TCGA_deconv$TCGA_LM,"COAD")
TCGA_SURV_info = pull_TCGA_surv()
TCGA_surv_dat = TCGA_SURV(TCGA_box_dat$TCGA_dat, TCGA_SURV_info$ALL_CLI, TCGA_SURV_info$ALL_MAN,"COAD")
canc.df = TCGA_hoyd(TCGA_surv_dat$SURV_d,"COAD")

# These analyses are not included in the final version of the paper.
# clinical = read.csv("../../data/clinical.csv",stringsAsFactors = F)
# ra_COAD = read.csv(file.path(paths$ra_COAD, "ra_COAD.csv"),stringsAsFactors = F)
# collected from website (https://globalrph.com/bacteriacat/gram-negative-bacteria/)
# gram_neg_list = read.csv("../../data/gram negative list.csv")%>%
#   mutate(gram.negative = trimws(gram.negative, which = c("both")))
# neut = TCGA_micro_neut(clinical, ra_COAD,gram_neg_list,TCGA_box_dat$dat_plot)

GTEx_Portal = read_csv("../../data/GTEX_v8_sample-location-key.csv")
GTEx_dat = pull_GTEx()
GTEx_box_dat = GTEx_box(GTEx_Portal,GTEx_dat)

save(TCGA_box_dat, file = "prepared-data/TCGA_box_data.rda")
save(TCGA_surv_dat, file = "prepared-data/TCGA_surv_dat.rda")
save(canc.df, file = "prepared-data/canc-df.rda")
# save(neut, file = "prepared-data/neut.rda")
save(GTEx_box_dat, file = "prepared-data/GTEx_box_dat.rda")


########################## Analyses related to mouse data

zymo.old <- zymoFauxDes()
volplot <- zymo.old$volplot
micmeans <- zymo.old$micmeans
mousedist <- mouseDistances()
bc.plot <- mousedist$bc.plot
euldat <- mousedist$euldat
vulcrange <- mouseVulcRange()
zymo.tissuedat <- mouseTissueBox()

save(volplot, micmeans, file = "prepared-data/zymo-old.rda")
save(mousedist, euldat, file = "prepared-data/mouse-distances.rda")
save(vulcrange, file = "prepared-data/mouse_vulc-range.rda")
save(zymo.tissuedat, file = "prepared-data/mouse_tissuebox.rda")

