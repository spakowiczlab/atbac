library(dada2)
library(ShortRead)

# Define locations, set up
path <- "/fs/scratch/PAS1479/atbac/atbac-fastq/" # CHANGE ME to the directory containing the fastq files after unzipping.
fns <- list.files(path, full.names = T)
fns <- unlist(lapply(fns, function(x) list.files(x, full.names = T)))

fastqs <- fns[grepl(".fastq.gz$", fns)]
fastqs <- sort(fastqs) # Sort ensures forward/reverse reads are in same order
fnFs <- fastqs[grepl("_R1", fastqs)] # Just the forward read files
fnRs <- fastqs[grepl("_R2", fastqs)] # Just the reverse read files
# Get sample names, assuming files named as so: SAMPLENAME_XXX.fastq
sample.names <- gsub(".*C_SpakowiczD_(.*)_V1C_.*", "\\1", fnFs)

filt_path <- file.path(path, "filtered")
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

# Dereplicate the fastqs
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Get error rates
dadaFs.lrn <- dada(derepFs, err=NULL, selfConsist = TRUE, multithread=TRUE)
errF <- dadaFs.lrn[[1]]$err_out

dadaRs.lrn <- dada(derepRs, err=NULL, selfConsist = TRUE, multithread=TRUE)
errR <- dadaRs.lrn[[1]]$err_out

# Save error rate for sanity plot
saveRDS(dadaFs.lrn, "~/Documents/repos/atbac/data/dada2_errs.RDS")

# Sample inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE) #merge paired reads

seqtab <- makeSequenceTable(mergers) #get seqtab
seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=TRUE) # No chimeras!
saveRDS(seqtab.nochim, "~/Documents/repos/atbac/data/seqtabNoC_update.RDS")

sum(seqtab.nochim)/sum(seqtab)