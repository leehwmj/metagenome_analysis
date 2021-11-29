# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.14")

# BiocManager::install("phyloseq")

library("phyloseq")
library("ggplot2")

?otu_table

otumat = matrix(sample(1:100, 100, replace = TRUE), nrow = 10, ncol = 10)
?sample
rownames(otumat) <- paste0("OTU", 1:nrow(otumat))
colnames(otumat) <- paste0("Sample", 1:ncol(otumat))

taxmat = matrix(sample(letters, 70, replace = TRUE), nrow = nrow(otumat), ncol = 7)
rownames(taxmat) <- rownames(otumat)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")


OTU = otu_table(otumat, taxa_are_rows = TRUE)
OTU
TAX = tax_table(taxmat)
TAX

physeq = phyloseq(OTU, TAX)
physeq

plot_bar(physeq, fill = "Family")

#########################################

sampledata = sample_data(data.frame(
  Location = sample(LETTERS[1:4], size=nsamples(physeq), replace=TRUE),
  Depth = sample(50:1000, size=nsamples(physeq), replace=TRUE),
  row.names=sample_names(physeq),
  stringsAsFactors=FALSE
))
sampledata

####################################

library("ape")
random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
plot(random_tree)


?import_biom
getwd()
# S7 <- system.file("extdata", "./singlem/evalue_05/oil_car_biom/oil_carbonate_otu_table.S1.10.ribosomal_protein_S7.biom", package = "phyloseq")
# import_biom(S7)

library(dplyr)
# oil_carbonate_otu <- read.csv("./singlem/evalue_20/oil_carbonate_e20_otu_table.csv")
oil_carbonate_otu <- read.table("./singlem/evalue_20/oil_carbonate_e20_otu_table.csv", sep = "\t", header = TRUE)

class(oil_carbonate_otu)

#L2_rplB만 검색
oil_carbonate_otu_L2 <- filter(oil_carbonate_otu, grepl("S1.1.ribosomal_protein_L2_rplB", gene))

#글자 찾아 바꾸기. 샘플 이름 뒤에 "_1" 삭제
oil_carbonate_otu_L2$sample <- gsub("_1", "", oil_carbonate_otu_L2$sample)
View(oil_carbonate_otu_L2$taxonomy)

#sample명, tax명 기준으로 num_hits 합
str(oil_carbonate_otu_L2)
oil_carbonate_otu_L2$num_hits <- as.numeric(oil_carbonate_otu_L2$num_hits)
oil_carbonate_otu_L2_sum <- oil_carbonate_otu_L2 |> group_by(sample, taxonomy) |> summarise(sum_num_hit = sum(num_hits))

install.packages("splitstackshape")
library(splitstackshape)
#Root만 있는 tax 제거
oil_carbonate_otu_L2_sum <- oil_carbonate_otu_L2_sum |> filter(!grepl("^Root$", taxonomy))

library(reshape2)
?cast
oil_carbo_otu_L2_cast <- dcast(oil_carbonate_otu_L2_sum, taxonomy~sample)
rownames(oil_carbo_otu_L2_cast) <- paste0("OTU", 1:length(oil_carbo_otu_L2_cast[,1]))

#otu hit만 골라서 otu table
oil_carbo_L2_otu_table <- oil_carbo_otu_L2_cast[,c(2,3,4)]
class(oil_carbo_L2_otu_table)
oil_carbo_L2_otu_table[is.na(oil_carbo_L2_otu_table)] <- 0
oil_carbo_L2_otu_table

#otu table을 matrix로
class(oil_carbo_L2_otu_table)
oil_L2_otu_mat <- as.matrix(oil_carbo_L2_otu_table)
str(oil_L2_otu_mat)
class(oil_L2_otu_mat)
oil_L2_otu_mat

#otu taxonomy table
class(oil_carbo_otu_L2_cast)
oil_carbo_otu_L2_cast_tax_split <- cSplit(oil_carbo_otu_L2_cast, splitCols = "taxonomy", sep = ";")
class(oil_carbo_otu_L2_cast_tax_split)
oil_carbo_otu_L2_cast_tax_split<- as.data.frame(oil_carbo_otu_L2_cast_tax_split)
class(oil_carbo_otu_L2_cast_tax_split)
rownames(oil_carbo_otu_L2_cast_tax_split) <- rownames(oil_carbo_L2_otu_table)

oil_carbo_L2_tax_table <- oil_carbo_otu_L2_cast_tax_split[,c(5,6,7,8,9,10)]
class(oil_carbo_L2_tax_table)
colnames(oil_carbo_L2_tax_table) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
oil_L2_tax_mat <- as.matrix(oil_carbo_L2_tax_table)
class(oil_L2_tax_mat)

oil_L2_OTU <- otu_table(oil_L2_otu_mat, taxa_are_rows = TRUE)
oil_L2_TAX <- tax_table(oil_L2_tax_mat)
oil_L2_OTU

oil_L2_physeq <- phyloseq(oil_L2_OTU, oil_L2_TAX)
oil_L2_physeq

plot_bar(oil_L2_physeq, fill = "Class")

ig <- make_network(oil_L2_physeq, max.dist = 0.4)
plot_network(ig, oil_L2_physeq)

data(enterotype)
