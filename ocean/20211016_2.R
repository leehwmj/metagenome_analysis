library(devtools)
library(igraph)
library(qgraph)
library(vegan)
library(MCL)
library(SpiecEasi)

data("amgut1.filt")

colnames(amgut1.filt) <- sapply(1:ncol(amgut1.filt),
                                function(x) paste("OTU", x, sep = "_"))
rownames(amgut1.filt) <- sapply(1:nrow(amgut1.filt),
                                function(x) paste("Sample", x, sep = "_"))

#OTU read
#user.table <- read.table(infile, header = T, row.names = 1, sep = ",")

# L2_otu_convert <- t(L2_otu)
# str(L2_otu_convert)

L2_otu_matrix<- as.matrix(L2_otu)
rownames(L2_otu_matrix) <- L2_otu_matrix[,1]
L2_otu_matrix <- L2_otu_matrix[,c(2:7)]
L2_otu_matrix[is.na(L2_otu_matrix)] <- 0
L2_otu_matrix
class(L2_otu_matrix) <- "numeric"
L2_otu_matrix <- L2_otu_matrix[c(1:429),]
L2_otu_convert <- t(L2_otu_matrix)
str(L2_otu_convert)

#otu ralative
otu.relative <- amgut1.filt / rowSums(amgut1.filt)
L2_otu.relative <- L2_otu_convert / rowSums(L2_otu_convert)

# Create dissimilarity matrix
distances <- vegdist(t(otu.relative), method = "bray")
L2_distances <- vegdist(t(L2_otu.relative), method = "bray")
