install.packages("remotes")
remotes::install_github("MadsAlbertsen/ampvis2")

install.packages("BiocManager")
BiocManager::install("biomformat")

library(ampvis2)

data("MiDAS")
MiDAS$metadata

getwd()
# if(TRUE) {
#   biom_otuable <- amp_import_biom("./singlem/evalue_05/oil_car_biom/oil_carbonate_otu_table.S1.10.ribosomal_protein_S7.biom")
#   otu_test <- amp_load(biom_otuable)
# }
# ?amp_import_biom()
# help("Defunct")
# biom_otuable <- amp_import_biom("./singlem/evalue_05/oil_car_biom/oil_carbonate_otu_table.S1.10.ribosomal_protein_S7.biom")
otu_test <- amp_load(otutable = "./singlem/evalue_05/oil_car_biom/oil_carbonate_otu_table.S1.10.ribosomal_protein_S7.biom")
otu_test$metadata[,2]<-c("Blank","PPC 4.2","PPC 1.2")
otu_test$metadata[,1]<-c("ECHO-B","ECHO-C","ECHO-F")
colnames(otu_test$metadata) <- c("sampleID","substance")
rownames( otu_test$metadata)<- c(1,2,3)
otu_test$metadata
otu_test
View(otu_test$tax)
View(otu_test$abund)
# View(MiDAS$tax)

otu_test$abund[1,]

rownames(otu_test$abund)[2405]
rownames(otu_test$tax)[2405]

paste("OTU", c(1,2,3), sep = "_")

OTU_names <- paste( "OTU", c( 1:length(rownames(otu_test$abund)) ), sep = "_" )
rownames(otu_test$abund) <- OTU_names
rownames(otu_test$tax) <- OTU_names
otu_test$tax[,8] <- OTU_names

# otu_subset <- amp_subset_samples(otu_test, substance %in% c("Blank", "PPC 1.2"))

amp_heatmap(otu_test, tax_aggregate="Class", tax_show = 25)
amp_heatmap(otu_test, tax_aggregate="Order", tax_show = 25)
amp_heatmap(otu_test, tax_aggregate="Family", tax_show = 25)
amp_heatmap(otu_test, tax_aggregate="Genus", tax_show = 25)
amp_heatmap(otu_test, tax_aggregate="Species", tax_show = 25)
?amp_heatmap()




