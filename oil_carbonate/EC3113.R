### KO = KEGG Orthology
oil_KO_tpm <- read.table(file = './HN00149294_RawData_squeezemeta_spades_0623/results/tables/HN00149294_RawData_squeezemeta_spades_0623.KO.tpm.tsv', sep = '\t', header = TRUE)
## KEGG number에 대응되는 name 수집
oil_KO_name <- read.table(file = './HN00149294_RawData_squeezemeta_spades_0623/results/tables/HN00149294_RawData_squeezemeta_spades_0623.KO.names.tsv', sep = '\t', header = TRUE, quote = "", row.names=NULL, stringsAsFactors = FALSE)

# EC4_list <- c("EC:4.1", "EC:4.2", "EC:4.3", "EC:4.4", "EC:4.5", "EC:4.6", "EC:4.7", "EC:4.99")
# EC4_list <- c("EC:4.1.1", "EC:4.1.2", "EC:4.1.3", "EC:4.1.99", "EC:4.1.1")
# EC4_list <- "EC:3.1.1.3"
# length(EC4_list)
# 
# EC4 <- matrix(c(0,0,0), nrow = 1, ncol = 3)
# 
# for(i in 1:length(EC4_list)){
#   print (paste("EC:4.", i, sep = ""))
#   if(i == 1){
#     EC4 <- oil_KO_name |> filter(grepl(EC4_list[i], Name)) |> inner_join(oil_KO_tpm) |> select(ECHO_Blank, ECHO_C, ECHO_F) |> as.matrix() |> colSums() |> as.matrix()
#   }else{
#     EC4 <- cbind(EC4, oil_KO_name |> filter(grepl(EC4_list[i], Name)) |> inner_join(oil_KO_tpm) |> select(ECHO_Blank, ECHO_C, ECHO_F) |> as.matrix() |> colSums() |> as.matrix())
#   }
#   colnames(EC4)[i] <- EC4_list[i]
# }


EC3113 <- data.frame(
  Name <- c("EC:3.1.1.3", "EC:3.1.1.3", "EC:3.1.1.3"),
  variable <- c("B", "C", "F"),
  value <- c(180.5046, 205.2182, 170.2210)
)

ggplot(EC3113, aes(x=variable, y=Name, fill=value))+
  geom_tile(colour="black")+
  scale_y_discrete(limits="EC:3.1.1.3")+
  scale_fill_gradient(low = "#F9E79F", high = "#E67E22") +
  labs(fill="TPM", x="", y="")+
  scale_x_discrete(labels=c("B"="Blank","C"="C_P8.9","F"="F_N33"))
