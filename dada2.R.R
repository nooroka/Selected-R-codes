library(dada2); packageVersion("dada2")
#setwd("~/16S clean reads-20230526T195949Z-001/16S clean reads")
path = "."
fnFs <- sort(list.files(path, pattern="_L1_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_L1_2.fastq", full.names = TRUE))
sample.names1 <-basename(fnFs)
sample.names2 <-basename(fnRs)
filtFs <- file.path(path, "filtered", paste0(sample.names1, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names2, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names1
names(filtRs) <- sample.names2
plotQualityProfile(fnFs)
plotQualityProfile(fnRs)
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(250,250),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
derepF1 <- derepFastq(filtFs, verbose=TRUE)
derepR1 <- derepFastq(filtRs, verbose=TRUE)
errF <- learnErrors(derepF1, multithread=TRUE)
errR <- learnErrors(derepR1, multithread=TRUE)
dadaF1 <- dada(derepF1, err=errF, multithread=FALSE)
dadaR1 <- dada(derepR1, err=errR, multithread=FALSE)
merger1 <- mergePairs(dadaF1, derepF1, dadaR1, derepR1, verbose=TRUE)
merger1.nochim <- removeBimeraDenovo(merger1, multithread=FALSE, verbose=TRUE)
seqtab <- makeSequenceTable(merger1)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE) #оч долго работает
write.csv(seqtab.nochim, "seqtab_nochim.csv", row.names=FALSE)
saveRDS(seqtab.nochim, "seqtab_nochim.rds")
#getN <- function(x) sum(getUniques(x))
#track <- cbind(out, sapply(dadaF1, getN), sapply(dadaR1, getN), sapply(merger1, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
#colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
#sample.names <- sapply(strsplit(basename( sample.names1), "_"), `[`, 2)
#rownames(track) <- sample.names
#head(track)
taxa <- assignTaxonomy(seqtab.nochim, "tax/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "tax/silva_species_assignment_v138.1.fa.gz")
write.csv(taxa, "taxa.csv", row.names=FALSE)
saveRDS(taxa, "taxa.rds")
