rm(list=ls())

R.version

BiocManager::install("dada2", version ="3.14")
library(dada2); packageVersion("dada2")

path <- "/analysis/hwlee/monomer"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE))
fnFs
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply( strsplit( basename(fnFs), "_1.f" ), `[`, 1 )
# as.list(sample.names)

plotQualityProfile(fnFs[c(22:24,25:27,10:12)]) #blanks
plotQualityProfile(fnFs[1:9]) #2nd 4A,5A,6A
plotQualityProfile(fnFs[13:21]) #3rd 4A, 5A, 6A

# as.list(fnRs)
plotQualityProfile(fnRs[c(22:24,25:27,10:12)]) #blanks
plotQualityProfile(fnRs[1:9]) #2nd 4A,5A,6A
plotQualityProfile(fnRs[13:21]) #3rd 4A, 5A, 6A

# Place filtered files in filtered/ subdirectory 
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
# View(filtFs)
# class(filtFs)

# ?filterAndTrim
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(300, 210), maxN=0, maxEE=c(2, 2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)
# On Windows set multithread=FALSE 
head(out)

#에러체크
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)



dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

# dadaFs[[1]]

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample 
head(mergers[[1]])



seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths 
table(nchar(getSequences(seqtab)))


seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
rownames(seqtab.nochim)[c(25:27)] <- c("1st_Blank_1", "1st_Blank_2", "1st_Blank_3")

getN <- function(x) sum(getUniques(x)) 
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs) 
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim") 
rownames(track) <- sample.names 
head(track)
track

?assignTaxonomy
taxa <- assignTaxonomy(seqtab.nochim, "/analysis/hwlee/tax/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)

taxa.print <- taxa
# Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)



library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")

theme_set(theme_bw())

samples.out <- rownames(seqtab.nochim)
# samples.out[c(25:27)] <- c("1st_Blank_1", "1st_Blank_2", "1st_Blank_3")
# subject <- sapply(strsplit(samples.out, "A_"), `[`, 1)
# strsplit(samples.out, "_")[[2]][2]
# class(strsplit(samples.out, "_"))
# strsplit(samples.out, "_")
# subject <- substr(samples.out, 1, nchar(samples.out)-2)
# times <- substr(samples.out, 1, 3)
sample_name <- substr(samples.out, 5, nchar(samples.out)-2)
# repetition <- substr(samples.out,nchar(samples.out),nchar(samples.out))
time_and_sample <- substr(samples.out, 1, nchar(samples.out)-2)
samdf <- data.frame(Sample_name=sample_name, Time_and_Sample=time_and_sample)
# samdf$When <- "Early"
# samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out


ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
# ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample
?prune_samples

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

# all(tax_table(ps) %in% colnames(seqtab.nochim)[-1])
# ?all

plot_richness(ps, x="Time_and_Sample", measures=c("InvSimpson","Shannon", "Simpson", "ACE", "Fisher"), color="Sample_name") +
  scale_x_discrete(limits=c("1st_Blank", "2nd_Blank", "3rd_Blank",
                            "2nd_4A", "3rd_4A",
                            "2nd_5A", "3rd_5A",
                            "2nd_6A", "3rd_6A"))
?plot_richness


top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:2000]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Time_and_Sample", fill="Phylum") +
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

plot_bar(ps.top20, x="Time_and_Sample", fill="Class") +
  geom_bar(aes(color=Class, fill=Class), stat="identity", position="stack")

plot_bar(ps, x="Time_and_Sample", fill="Phylum")+
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")+
  scale_x_discrete(limits=c("1st_Blank", "2nd_Blank", "3rd_Blank",
                            "1st_Blank", "2nd_4A", "3rd_4A",
                            "1st_Blank", "2nd_5A", "3rd_5A",
                            "1st_Blank", "2nd_6A", "3rd_6A"))

top10 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:10]
ent10 <- prune_taxa(top20, ps)
plot_bar(ent10, "Time_and_Sample", fill="Phylum")+
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")
