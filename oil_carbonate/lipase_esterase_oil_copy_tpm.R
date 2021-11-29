oil_KO_copy <- oil_spades$functions$KEGG$copy_number
oil_KO_copy <- as.data.frame(oil_KO_copy)
oil_KO_copy$KO_number <- rownames(oil_KO_copy)
rownames(oil_KO_copy) <- c(1:length(oil_KO_copy[[1]]))
oil_KO_copy<- oil_KO_copy[,c(4,1,2,3)]
oil_KO_copy_join <- inner_join(KO_DB, oil_KO_copy, by='KO_number')
oil_KO_copy_esterase <- oil_KO_copy_join[grep("esterase", oil_KO_copy_join$definition),] |> filter(ECHO_Blank <= ECHO_C) |> 
  filter(ECHO_C <= ECHO_F)
oil_KO_copy_lipase <- oil_KO_copy_join[grep("lipase", oil_KO_copy_join$definition),] |> filter(ECHO_Blank <= ECHO_C) |> 
  filter(ECHO_C <= ECHO_F)

oil_KO_tpm <- oil_spades$functions$KEGG$tpm
oil_KO_tpm <- as.data.frame(oil_KO_tpm)
oil_KO_tpm$KO_number <- rownames(oil_KO_tpm)
rownames(oil_KO_tpm) <- c(1:length(oil_KO_tpm[[1]]))
oil_KO_tpm<- oil_KO_tpm[,c(4,1,2,3)]
oil_KO_tpm_join <- inner_join(KO_DB, oil_KO_tpm, by='KO_number')
oil_KO_tpm_esterase <- oil_KO_tpm_join[grep("esterase", oil_KO_tpm_join$definition),] |> filter(ECHO_Blank <= ECHO_C) |> 
  filter(ECHO_C <= ECHO_F)
oil_KO_tpm_lipase <- oil_KO_tpm_join[grep("lipase", oil_KO_tpm_join$definition),] |> filter(ECHO_Blank <= ECHO_C) |> 
  filter(ECHO_C <= ECHO_F)


 <- oil_spades$functions$KEGG$tpm
