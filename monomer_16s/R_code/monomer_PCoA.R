#DADA2 에서 이어짐
rm(list = ls())
getwd()

## downloading DECIPHER-formatted SILVA v138 reference
download.file(url="http://www2.decipher.codes/Classification/TrainingSets/SILVA_SSU_r138_2019.RData", destfile="SILVA_SSU_r138_2019.RData")

## loading reference taxonomy object
load("SILVA_SSU_r138_2019.RData")


## loading DECIPHER
library(DECIPHER)
packageVersion("DECIPHER")

## creating DNAStringSet object of our ASVs
dna <- DNAStringSet(getSequences(seqtab.nochim))

## and classifying
tax_info <- IdTaxa(test = dna, trainingSet=trainingSet, strand = "both", processors = NULL)

########

# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep = "_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
# setwd("/analysis/hwlee/monomer/ASV_to_beta")
# getwd()
write(asv_fasta, "ASVs.fa")

#count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.tsv", sep = "\t", quote = F, col.names = NA)

#tax table:
# creating table of taxonomy and setting any that are unclassified as "NA"
ranks <- c("domain","phylum","class","order","family","genus","species")
asv_tax <- t(sapply(tax_info, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(asv_tax) <- ranks
rownames(asv_tax) <- gsub(pattern = ">", replacement = "", x=asv_headers)

write.table(asv_tax, "ASVs_taxonomy.tsv", sep = "\t", quote=F, col.names = NA)

###

# BiocManager::install("decontam")
library(decontam)
#########
colnames(asv_tab)
vector_for_decontam <- c(rep(TRUE, 4), rep(FALSE, 16))

contam_df <- isContaminant(t(asv_tab), neg=vector_for_decontam)

table(contam_df$contaminant)

# getting vector holding the identified contaminant IDs
contam_asvs <- row.names(contam_df[contam_df$contaminant == TRUE, ])

asv_tax[row.names(asv_tax) %in% contam_asvs, ]

#################

# making new fasta file
contam_indices <- which(asv_fasta %in% paste0(">", contam_asvs))
dont_want <- sort(c(contam_indices, contam_indices + 1))
asv_fasta_no_contam <- asv_fasta[- dont_want]

# making new count table
asv_tab_no_contam <- asv_tab[!row.names(asv_tab) %in% contam_asvs, ]

# making new taxonomy table
asv_tax_no_contam <- asv_tax[!row.names(asv_tax) %in% contam_asvs, ]

## and now writing them out to files
write(asv_fasta_no_contam, "ASVs-no-contam.fa")
write.table(asv_tab_no_contam, "ASVs_counts-no-contam.tsv",
            sep="\t", quote=F, col.names=NA)
write.table(asv_tax_no_contam, "ASVs_taxonomy-no-contam.tsv",
            sep="\t", quote=F, col.names=NA)

##########################################################################

library(tidyverse) ; packageVersion("tidyverse") # 1.3.1
library(phyloseq) ; packageVersion("phyloseq") # 1.38.0
library(vegan) ; packageVersion("vegan") # 2.5.7
# BiocManager::install("DESeq2")
library(DESeq2) ; packageVersion("DESeq2") # 1.34.0
# BiocManager::install("dendextend")
library(dendextend) ; packageVersion("dendextend") # 1.15.2
library(viridis) ; packageVersion("viridis") # 0.6.1

# rm(list = ls())


count_tab <- read.table("ASVs_counts-no-contam.tsv", header=T, row.names=1,
                        check.names=F, sep="\t")

tax_tab <- as.matrix(read.table("ASVs_taxonomy-no-contam.tsv", header=T,
                                row.names=1, check.names=F, sep="\t"))

# sample_info_tab <- read.table("sample_info.tsv", header=T, row.names=1,
                              # check.names=F, sep="\t")
# rm(sample.out)
samples.out <- colnames(count_tab)
# subject <- substr(samples.out, 1, nchar(samples.out)-2)
times <- as.factor(substr(samples.out, 1, 3))
sample_name <- as.factor(substr(samples.out, 5, nchar(samples.out)-2))
repetition <- as.factor(substr(samples.out,nchar(samples.out),nchar(samples.out)))
time_and_sample <- as.factor(substr(samples.out, 1, nchar(samples.out)-2))
color <- c("red", "red", "red", "blue", "blue", "blue", "forestgreen", "forestgreen", "forestgreen", "black", "black", "black",
           "red", "red", "red", "blue", "blue", "blue", "forestgreen", "forestgreen", "forestgreen", "black", "black", "black",
           "black", "black", "black")
sample_info_tab <- data.frame(Times=times, Sample_name=sample_name, Repetition=repetition, Time_and_Sample=time_and_sample, color=color)

str(sample_info_tab)


# and setting the color column to be of type "character", which helps later
sample_info_tab$color <- as.character(sample_info_tab$color)

sample_info_tab # to take a peek


# first we need to make a DESeq2 object
deseq_counts <- DESeqDataSetFromMatrix(count_tab, colData = sample_info_tab, design = ~Time_and_Sample)
# we have to include the "colData" and "design" arguments because they are 
# required, as they are needed for further downstream processing by DESeq2, 
# but for our purposes of simply transforming the data right now, they don't 
# matter

deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
# NOTE: If you get this error here with your dataset: "Error in
# estimateSizeFactorsForMatrix(counts(object), locfunc =locfunc, : every
# gene contains at least one zero, cannot compute log geometric means", that
# can be because the count table is sparse with many zeroes, which is common
# with marker-gene surveys. In that case you'd need to use a specific
# function first that is equipped to deal with that. You could run:
# deseq_counts <- estimateSizeFactors(deseq_counts, type = "poscounts")
# now followed by the transformation function:
# deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)

# and here is pulling out our transformed table
vst_trans_count_tab <- assay(deseq_counts_vst)

# and calculating our Euclidean distance matrix
euc_dist <- dist(t(vst_trans_count_tab))

###################

euc_clust <- hclust(euc_dist, method="ward.D2")

# hclust objects like this can be plotted with the generic plot() function
plot(euc_clust) 
# but i like to change them to dendrograms for two reasons:
# 1) it's easier to color the dendrogram plot by groups
# 2) if wanted you can rotate clusters with the rotate() 
#    function of the dendextend package

euc_dend <- as.dendrogram(euc_clust, hang=0.1)
dend_cols <- as.character(sample_info_tab$color[order.dendrogram(euc_dend)])
labels_colors(euc_dend) <- dend_cols

plot(euc_dend, ylab="VST Euc. dist.")

################################################
#####
rownames(sample_info_tab) <- colnames(count_tab)
#####
vst_count_phy <- otu_table(vst_trans_count_tab, taxa_are_rows=T)
sample_info_tab_phy <- sample_data(sample_info_tab)
vst_physeq <- phyloseq(vst_count_phy, sample_info_tab_phy)
# sample_names(vst_count_phy)
# sample_names(sample_info_tab_phy)

# generating and visualizing the PCoA with phyloseq
vst_pcoa <- ordinate(vst_physeq, method="MDS", distance="euclidean")
eigen_vals <- vst_pcoa$values$Eigenvalues # allows us to scale the axes according to their magnitude of separating apart the samples

plot_ordination(vst_physeq, vst_pcoa, color="color") + 
  geom_point(size=1) + labs(col="type") + 
  geom_text(aes(label=rownames(sample_info_tab), hjust=0.3, vjust=-0.4)) + 
  coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle("PCoA") + 
  scale_color_manual(values=unique(sample_info_tab$color[order(sample_info_tab$color)])) + 
  theme_bw() + theme(legend.position="none")

#####################################################

rarecurve(t(count_tab), step=100, col=sample_info_tab$color, lwd=2, ylab="ASVs", label=F)

# and adding a vertical line at the fewest seqs in any sample
abline(v=(min(rowSums(t(count_tab)))))
