library(readr)
antarctic_ocean_otu_table_single <- read_delim("R_project/ocean/antarctic_ocean_otu_table_single.csv", 
                                               delim = "\t", escape_double = FALSE, 
                                               trim_ws = TRUE)


View(antarctic_ocean_otu_table_single)

antarctic_ocean_otu_table <- read_delim("R_project/ocean/antarctic_ocean_otu_table.csv", 
                                               delim = "\t", escape_double = FALSE, 
                                               trim_ws = TRUE)

# L2 
L2_otu <- antarctic_ocean_otu_table |> filter(gene == "S1.1.ribosomal_protein_L2_rplB")
L2_16_1_L1_otu <- L2_otu |> filter(sample == "16-1_TTAGGC_L001_R1_001")
L2_16_2_L1_otu <- L2_otu |> filter(sample == "16-2_ACTTGA_L001_R1_001")
L2_16_1_L2_otu <- L2_otu |> filter(sample == "16-1_TTAGGC_L002_R1_001")
L2_16_2_L2_otu <- L2_otu |> filter(sample == "16-2_ACTTGA_L002_R1_001")
L2_25_1_L1_otu <- L2_otu |> filter(sample == "25_GGCTAC_L001_R1_001")
L2_25_1_L2_otu <- L2_otu |> filter(sample == "25_GGCTAC_L002_R1_001")

# L2_25_otu <- full_join(L2_25_1_L1_otu[,c(6,4)], L2_25_1_L2_otu[,c(6,4)], by = "taxonomy")
# L2_16_1_otu <- full_join(L2_16_1_L1_otu[,c(6,4)], L2_16_1_L2_otu[,c(6,4)], by = "taxonomy")
# L2_16_2_otu <- full_join(L2_16_2_L1_otu[,c(6,4)], L2_16_2_L2_otu[,c(6,4)], by = "taxonomy")
# a <- L2_25_1_L1_otu[,c(6,4)]
# b <- a |> group_by(taxonomy) |> summarise(sum_hits = sum(num_hits))

L2_25_1_L1 <- L2_25_1_L1_otu[,c(6,4)] |> group_by(taxonomy) |> summarise(num_hits = sum(num_hits))
L2_25_1_L2 <- L2_25_1_L2_otu[,c(6,4)] |> group_by(taxonomy) |> summarise(num_hits = sum(num_hits))
L2_16_1_L1 <- L2_16_1_L1_otu[,c(6,4)] |> group_by(taxonomy) |> summarise(num_hits = sum(num_hits))
L2_16_1_L2 <- L2_16_1_L2_otu[,c(6,4)] |> group_by(taxonomy) |> summarise(num_hits = sum(num_hits))
L2_16_2_L1 <- L2_16_2_L1_otu[,c(6,4)] |> group_by(taxonomy) |> summarise(num_hits = sum(num_hits))
L2_16_2_L2 <- L2_16_2_L2_otu[,c(6,4)] |> group_by(taxonomy) |> summarise(num_hits = sum(num_hits))

L2_25_otu <- full_join(L2_25_1_L1, L2_25_1_L2, by = "taxonomy")
L2_16_1_otu <- full_join(L2_16_1_L1, L2_16_1_L2, by = "taxonomy")
L2_16_2_otu <- full_join(L2_16_2_L1, L2_16_2_L2, by = "taxonomy")


#
colnames(L2_16_1_otu) <- c("taxonomy","L2_16_1_L1", "L2_16_1_L2")
colnames(L2_16_2_otu) <- c("taxonomy","L2_16_2_L1", "L2_16_2_L2")
colnames(L2_25_otu) <- c("taxonomy","L2_25_L1", "L2_25_L2")

L2_16_otu <- full_join(L2_16_1_otu, L2_16_2_otu, by="taxonomy")
L2_otu <- full_join(L2_16_otu, L2_25_otu, by="taxonomy")

