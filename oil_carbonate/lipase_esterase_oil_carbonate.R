rm(a)
getwd()
system.time(load("./RDatas/oil_spades.RData"))

oil_spades$functions$KEGG$abund
oil_KO_abund <- oil_spades$functions$KEGG$abund

oil_KO_abund <- as.data.frame(oil_KO_abund)
oil_KO_abund$KO_number <- rownames(oil_KO_abund)
oil_KO_abund<- oil_KO_abund[,c(4,1,2,3)]
rownames(oil_KO_abund) <- c(1:12479)
oil_KO_abund_join <- inner_join(KO_DB, oil_KO_abund, by='KO_number')
total_esterase <- oil_KO_abund_join[grep("esterase", oil_KO_abund_join$definition),]

# blank보다 큰 C, 그리고 C보다 큰 F를 만족하는 행 중에, abund가 100이 넘는 행만 filter
oil_KO_abund_filter <- oil_spades$functions$KEGG$abund |> as.data.frame() |> filter(ECHO_Blank <= ECHO_C) |> 
  filter(ECHO_C <= ECHO_F) |> filter(ECHO_F > 100)
oil_KO_abund_filter2 <- oil_spades$functions$KEGG$abund |> as.data.frame() |> filter(ECHO_Blank <= ECHO_C) |> 
  filter(ECHO_C <= ECHO_F)
oil_KO_abund_filter2["KO_number"] <- rownames(oil_KO_abund_filter2)
oil_KO_abund_filter2 <- oil_KO_abund_filter2[,c(4,1,2,3)]

#unmapped 제거
oil_KO_abund_filter <-oil_KO_abund_filter[!(oil_KO_abund_filter$ECHO_Blank >1000000),]
# class(oil_KO_abund_filter)

KO_DB <- read.csv("./oil_carbonate/ko00001.keg.csv")


oil_KO_abund_filter["KO_number"] <- rownames(oil_KO_abund_filter)
# rownames(oil_KO_abund_filter) = c(1:607)
# oil_KO_abund_filter <- oil_KO_abund_filter[,c(4,1,2,3)]

kegg_abund_oil_carbonate <- inner_join(oil_KO_abund_filter, KO_DB, by='KO_number')

lipase <- kegg_abund_oil_carbonate[grep("lipase", kegg_abund_oil_carbonate$definition),c(4,5,1,2,3,6)]
esterase <- kegg_abund_oil_carbonate[grep("esterase", kegg_abund_oil_carbonate$definition),c(4,5,1,2,3,6)]
polyurethanase <- kegg_abund_oil_carbonate[grep("polyurethanase", kegg_abund_oil_carbonate$definition),c(4,5,1,2,3,6)]

##

write.csv(kegg_abund_oil_carbonate[,c(4,5,1,2,3,6)], 'oil_carbonate_kegg_abund.csv')

fatty_acid_degra_kegg <- c("K00001,K00022,K00121,K00128,K00149,K00232,K00248,K00249,K00252,K00255,K00496,K00529,K00626,K00632,K01692,K01782,K01825,K01897,K01909,K04072,K05297,K05939,K06445,K07422,K07425,K07508,K07509,K07511,K07513,K07514,K07515,K07516,K07517,K08765,K08766,K09478,K09479,K10527,K13238,K13239,K13767,K13951,K13952,K13953,K13954,K13980,K14085,K14338,K15013,K15401,K18857,K18880,K19523,K19524,K20495,K21738,K22567,K22887")
class(fatty_acid_degra_kegg)
fatty_degra_kegg_split <- strsplit(fatty_acid_degra_kegg, split= ",") |> as.data.frame()
colnames(fatty_degra_kegg_split)[1] <- "KO_number"

fatty_degra_in_join <- inner_join(fatty_degra_kegg_split, KO_DB, by='KO_number')
fatty_degra_in_join_abund <- inner_join(fatty_degra_in_join, oil_KO_abund_filter2, by='KO_number')


