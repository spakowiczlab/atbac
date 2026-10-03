<p align="right">
<img src="images/atbac-hex.png" alt="atbac hex sticker" width="180" />
</p>

# atbac [![DOI](https://zenodo.org/badge/271302970.svg)](https://zenodo.org/badge/latestdoi/271302970)

Scripts and data for the analyses in:

> Shantaram D\*, Hoyd R\*, Blaszczak AM\*, Antwi L, Jalilvand A, Wright VP, Liu J, Smith AJ, Bradley D, Lafuse W, Liu Y, Williams NF, Snyder O, Wheeler C, Needleman B, Brethauer S, Noria S, Renton D, Perry KA, Nagareddy P, Wozniak D, Mahajan S, Rana PSJB, Pietrzak M, Schlesinger LS, Spakowicz DJ, Hsueh WA. Obesity-associated microbiomes instigate visceral adipose tissue inflammation by recruitment of distinct neutrophils. *Nature Communications* **15**, 5434 (2024). [https://doi.org/10.1038/s41467-024-48935-5](https://doi.org/10.1038/s41467-024-48935-5). PMID: [38937454](https://pubmed.ncbi.nlm.nih.gov/38937454/).
>
> \*These authors contributed equally.

![Graphical abstract. Visceral adipose tissue from people with obesity carries more neutrophils and a distinct bacterial community. In mice, only a high-fat diet plus stool from a donor with obesity recruits those neutrophils, which is followed by more Th1 cells and fewer regulatory T cells. The VAT neutrophil signature is distinct from blood neutrophils, is widespread across human tissues, and separates survival in colon cancer.](images/graphical-abstract.png)

Sequencing data generated for the paper are in NCBI BioProject [PRJNA766535](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA766535). A snapshot of this repository is archived on Zenodo via the DOI badge above.

## Reproducing the analyses

The notebooks under `manuscript/scripts/` are the published analysis. `exploratory/` and `grants/` are earlier or unrelated work. Flow cytometry, ELISA, qRT-PCR, histology, and the mouse immune-cell time courses (Figures 1, 3C–E, 3H–I, 4, and 6B) were analyzed in FlowJo and GraphPad Prism and are not in this repository. Cell-fraction estimates were produced on the [CIBERSORTx](https://cibersortx.stanford.edu/) website; the signature matrix and the result tables are committed here.

Knit processing notebooks from `manuscript/scripts/`. Knit figure notebooks from `manuscript/scripts/drake-figs/`. Paths inside the notebooks are relative to those directories.

### Environment

The paper was run in R 4.1.0, with some early processing in R 4.0.2. Packages named in the code-availability statement:

| Package | Version | Role |
| --- | --- | --- |
| dada2 | 1.12.1 | ASV inference |
| decontam | 1.4.0 | Negative-control filtering |
| HMP16SData | 1.8.2 | Human Microbiome Project source communities |
| MGnifyR | 0.1.0 | Built-environment source communities |
| DESeq2 | — | Taxon and gene differential abundance |
| tidyverse, readr, broom | 2.0.0 / 2.1.4 / 1.0.4 | Data wrangling |
| vegan | 2.5.7 | Diversity and distances |
| survival, survminer | 3.3.1 / 0.4.9 | Kaplan–Meier curves |
| Seurat | 4.4.0 | Single-cell neutrophil subset |
| tmesig | 0.1.0 | VIN-like score (`devtools::install_github("spakowiczlab/tmesig")`) |
| GenomicDataCommons | 1.15.0 | TCGA RNA-seq download |
| ggdendro, RColorBrewer, viridis, ggrepel, ggfortify, ggforce, flextable | as in the paper | Figures and supplementary tables |

SourceTracker 2.0.1 was run outside R. CIBERSORTx replaced the original CIBERSORT website. `MGnifyR` was installed from `beadyallen/MGnifyR` at the time of the analysis (the call is commented at the top of `processing_sourcetracker.Rmd`).

### Redraw the published panels from committed results

`manuscript/scripts/drake-figs/prepared-data/` already contains the objects the figure notebooks load. From `manuscript/scripts/drake-figs/`, knit:

| Notebook | Paper panels |
| --- | --- |
| `figure_2.Rmd` | 2B phylum bars, 2C SourceTracker bars, 2D taxon volcano |
| `figure_3.Rmd` | 3B tissue bacterial load, 3F–G mouse VAT taxon models |
| `figure_5.Rmd` | 5A–E neutrophil expression |
| `figure_6.Rmd` | 6A heatmap, 6C validation deconvolution, 6D GTEx, 6E–F TCGA-COAD, 6G single-cell score |
| `supplementary.Rmd` | Supplementary Table 2 (public neutrophil cohorts) |
| `supplementary_diversity_human-16s.Rmd` | Human VAT richness and Simpson diversity |
| `supplementary_saline-microbes.Rmd` | Saline-gavage VAT taxa |
| `bray-curtis_mouse-tissues.Rmd` | Mouse community distances and engraftment |

To rebuild those `.rda` files from the committed tables, run `source("prepare-all-plots.R")` in that same directory, then knit. The script sources every file in `analysis-functions/`.

### Figure 2B–D: bacteria in human VAT

1. Download the 16S reads from BioProject PRJNA766535.
2. In `processing_dada2.Rmd`, set `path` to the directory of demultiplexed fastqs. Reads are quality- and length-filtered to 340–440 bp, denoised with dada2 1.12.1, and assigned taxonomy against SILVA v123. The chimera-filtered table committed in the clone is `manuscript/data/atbac/seqtabNoC_update.RDS`. The length filter at the end of that notebook writes `seqtabNoCf_update.RDS`. The copy used by later steps is `exploratory/data/seqtabNoCf_update.RDS`. `processing_decontam_test-for-contaminants.Rmd` and `DESeq2-and-fisher_microbiome.Rmd` look for it at `manuscript/data/seqtabNoCf_update.RDS`, so copy it there before re-running them. Taxonomy is `manuscript/data/atbac/taxf_update.RDS`.
3. `processing_decontam_test-for-contaminants.Rmd` runs decontam, then drops any genus seen in the sequenced lysis blanks. If a negative-control ASV has no genus, sample ASVs at least 97% identical to it are removed. The notebook also reads `manuscript/data/Sample list from Nyelia_10092020.xlsx`, which is not in this clone. The interface-layer DNA table after that filter is already committed as `manuscript/data/atbac/taxcounts_interface-DNA_no-neg-contams.RDS`, and BMI labels are in `samples_obesity-status.csv`.
4. `processing_sourcetracker.Rmd` builds source communities from HMP16SData `V35()` body sites and the MGnify built-environment biome, then calls SourceTracker2. Mixing proportions are committed as `manuscript/data/atbac/mixing_proportions.txt` (Figure 2C: gastrointestinal tract is the dominant source).
5. `DESeq2-and-fisher_microbiome.Rmd` tests obese versus lean VAT at several taxonomic ranks. The result used in Figure 2D is `manuscript/data/atbac/DESeq2_all-levs_DNA.csv`. In that table, Streptococcaceae and *Ruminococcaceae_UCG-014* are higher in VAT from donors with obesity, and the order Bacillales and the genus *Marvinbryantia* are higher in lean VAT.
6. `prepare-all-plots.R` orders samples and saves `prepared-data/sample-bmi-firmicutes.rda`, `taxcounts-nocontam-phyl.rda`, `sourceres.rda`, and `bacteria-desres.rda`. `figure_2.Rmd` draws the three panels.

### Figure 3B, 3F–G and the mouse 16S supplements

Zymo-processed ASV tables for the avatar-mouse experiments are under `manuscript/data/mouse-data/first-experiment/` and `second-experiment/` (abundance, taxonomy, and absolute copies). The raw delivery directory `zr8258.16S_220923.zymo` is gitignored.

`processing-zymo-new-for-vulc.Rmd` fits the univariate gamma generalized linear model (Figure 3F, one obese and one lean donor) and the logistic models that compare each obese donor with lean donors (Figure 3G). `prepare-all-plots.R` then calls `mouseTissueBox()`, `zymoFauxDes()`, `mouseVulcRange()`, and `mouseDistances()`, and `figure_3.Rmd` draws bacterial load by tissue (Figure 3B) and the taxon plots. Bray–Curtis distances, saline-gavage taxa, and diversity are in `bray-curtis_mouse-tissues.Rmd` and `supplementary_saline-microbes.Rmd`.

### Figure 5: VAT neutrophils versus other activation states

Human Ampliseq counts for paired VAT and blood neutrophils are `manuscript/data/ncbi_processed/iso_bldvat.txt` with `meta_bldvat.RDS`. Public comparator cohorts were unpacked with GEOquery in `processing_unpack-geo-objects.Rmd` and `processing_unpack-syno.Rmd` from GSE2322, GSE19443, GSE64457, GSE8668, and GSE116899 (sepsis, exercise, endotoxin, active tuberculosis, synovial fluid, and airspace; Supplementary Table 2). Processed matrices and sample metadata are under `manuscript/data/ncbi_processed/` and `manuscript/data/ncbi_raw-files/`.

`DESeq2_expressions.Rmd` runs one DESeq2 contrast per condition, with Benjamini–Hochberg adjustment, and writes `manuscript/data/DESeq2/*.csv`. Pathway membership used for the boxplots is `manuscript/data/from-alecia/pathgenes.csv`. DAVID over-representation inputs and outputs are in `manuscript/data/DAVID/`; DAVID itself is the web tool, not an R step.

`prepare-all-plots.R` builds the PCA, pathway, shared-DEG, NMDS, and signature-score objects. `figure_5.Rmd` draws them: global PCA (5A), activation-pathway counts (5B), shared differentially expressed genes and their Spearman correlation (5C), ordination of the isolated-neutrophil datasets (5D), and pathway signature scores (5E). `supplementary.Rmd` formats the cohort table.

### Figure 6A and 6C: VIN signature and deconvolution check

The custom signature replaces the LM22 neutrophil column with the blood and VAT neutrophils from this study. The signature gene matrix is `manuscript/data/CIBERSORT/sig-genes/custom-sig-genes_bldat.txt`. Bulk blood and adipose cohorts listed in Supplementary Table 4 were deconvolved on CIBERSORTx. Those result tables are under `manuscript/data/CIBERSORT/CIBERSORTx-deconvolved/` and `manuscript/data/CIBERSORT/deconvolved-samples/`.

`prepare-all-plots.R` clusters the signature genes (Euclidean distance, complete linkage) and combines the validation fractions. `figure_6.Rmd` draws the heatmap and dendrogram (6A) and the blood-versus-VAT stacked bars (6C). VIN-type fractions are confined to adipose samples in that check.

### Figure 6D: GTEx

GTEx v8 tissue labels are `manuscript/data/GTEX_v8_sample-location-key.csv`. Per-batch CIBERSORTx output is `manuscript/data/CIBERSORT/CIBERSORTx-deconvolved/GTEx/` (`*_CUS.csv` for the custom signature and `*_LM22.csv` for the LM22 neutrophil). The expression matrices themselves are not in the clone.

`pull_GTEx()` and `GTEx_box()` in `prepare-all-plots.R` join those tables to tissue and mark whether the VIN-type or blood signature is higher. `figure_6.Rmd` section D draws the boxplots. `GTEx_neut_analysis.Rmd` is an earlier version of the same join.

### Figure 6E–F: TCGA colon adenocarcinoma

Tumor RNA-seq was retrieved with GenomicDataCommons, filtered to colon adenocarcinoma (and, for the negative checks, other obesity-associated cancers). Those matrices are not committed. CIBERSORTx results are `manuscript/data/CIBERSORT/CIBERSORTx-deconvolved/COAD-READ_cus.csv` and `COAD-READ_lm.csv`, plus BRCA, SARC, LUSC, LUAD, and KIRC tables in the same directory. File identifiers are linked with `manuscript/data/TCGA-link-bam-and-expression-files.csv` and the GDC manifests `gdc_manifest.2020-05-28.txt` and `gdc_manifest.2020-05-28 _READ.txt`. Clinical follow-up is in `manuscript/data/survival info/`.

Two equivalent entry points:

- `Fig6_B-C_TCGA-COAD_abundance-and-survival.Rmd`, knitted from `manuscript/scripts/` (the filename uses the panel letters from an earlier draft).
- `prepare-all-plots.R`, which calls `pull_TCGA_deconvolved()`, `TCGA_box()`, `pull_TCGA_surv()`, and `TCGA_SURV()` for COAD, then `figure_6.Rmd` sections E and F.

Samples are split into high and low VIN-type abundance for the Kaplan–Meier curve. The same split on blood-neutrophil abundance does not separate survival. The same contrast was run for other cancers and is not a result in the paper. A join of TCGA microbes to neutrophil fractions is left commented out in `prepare-all-plots.R` because it was not included in the final manuscript. `ONLY TCGA neut analysis.Rmd` is an earlier COAD-only pass.

### Figure 6G: single-cell lung tumors

The public object is the lung-tumor immune-cell atlas of Prazanowska and Lim ([figshare](https://doi.org/10.6084/m9.figshare.c.6222221.v3)). `processing_scrna-seurat.Rmd` subsets cells labeled as neutrophils. That notebook reads a Seurat file from an Ohio Supercomputer Center path that is not in the clone. The neutrophil table the figure uses is committed as `manuscript/data/scrna-neutrophils.rds`.

`sigscoreSCRNA()` takes the ten signature genes with the largest VAT-minus-blood difference and scores each cell with `tmesig::calculateAvgZScore`. `figure_6.Rmd` section G plots cells with score above 0.2 as VIN-like and the remainder as blood-like, alongside all neutrophils.

### What is already enough to redraw a panel

If a `prepared-data/*.rda` file listed above is present, knit the matching figure notebook and skip the upstream download. Re-running dada2, SourceTracker, DESeq2, or CIBERSORTx is only required when you want to regenerate those intermediates from reads or from expression matrices.
