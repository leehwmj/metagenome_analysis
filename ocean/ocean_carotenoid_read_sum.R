library(dplyr)

ocean_orfs_carotenoid <- ocean_orfs[grep('Carotenoid',ocean_orfs$KEGGPATH),]

colnames(ocean_orfs_carotenoid) <- c("Contig_ID", "Molecule", "Method",
                                "Length_NT", "Length_AA", "GC_perc", "Gene_name", "Tax", "KEGG_ID",
                                "KEGGFUN", "KEGGPATH", "COG_ID", "COGFUN", "COGPATH", "PFAM", "TPM_16depth",
                                "TPM_25depth", "Coverage_16depth", "Coverage_25depth", "Raw_read_count_16depth",
                                "Raw_read_count_25depth", "Raw_base_count_16depth", "Raw_base_count_25depth", "Hits")


ocean_orfs_carotenoid_16sum <-ddply(ocean_carotenoid, .(Gene_name, KEGG_ID),summarize, Raw_read_count_16depth =sum(Raw_read_count_16depth))
ocean_orfs_carotenoid_25sum <-ddply(ocean_carotenoid, .(Gene_name),summarize, Raw_read_count_25depth =sum(Raw_read_count_25depth))
?ddply

ocean_carotenoid_sum <- full_join(ocean_orfs_carotenoid_16sum, ocean_orfs_carotenoid_25sum, by="Gene_name")

