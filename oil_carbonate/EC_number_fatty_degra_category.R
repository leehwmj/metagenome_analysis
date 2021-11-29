getwd()
setwd("/analysis/hwlee/collect-temp-sequencing/HN00149294_Rawdata")

### KO = KEGG Orthology
oil_KO_tpm <- read.table(file = './HN00149294_RawData_squeezemeta_spades_0623/results/tables/HN00149294_RawData_squeezemeta_spades_0623.KO.tpm.tsv', sep = '\t', header = TRUE)
## KEGG number에 대응되는 name 수집
oil_KO_name <- read.table(file = './HN00149294_RawData_squeezemeta_spades_0623/results/tables/HN00149294_RawData_squeezemeta_spades_0623.KO.names.tsv', sep = '\t', header = TRUE, quote = "", row.names=NULL, stringsAsFactors = FALSE)

fatty_deg_EC_list <- c("EC:1.2.1.48", "EC:1.1.1.192", "EC:6.2.1.3", "EC:2.3.1.21", "EC:1.3.3.6", "EC:1.3.8.7",
                       "EC:1.3.99", "EC:1.3.8.8", "EC:1.3.8.9", "EC:4.2.1.17", "EC:4.2.1.74", "EC:1.1.1.35", 
                       "EC:1.1.1.211", "EC:2.3.1.16", "EC:2.3.1.9", "EC:6.2.1.6", "EC:1.3.8.6")

library(dplyr)

for(i in 1:length(fatty_deg_EC_list)){
  print (fatty_deg_EC_list[i])
  if(i == 1){
    fatty_deg <- oil_KO_name |> filter(grepl(fatty_deg_EC_list[i], Name)) |> inner_join(oil_KO_tpm) |> select(ECHO_Blank, ECHO_C, ECHO_F) |> as.matrix() |> colSums() |> as.matrix()
  }else{
    fatty_deg <- cbind(fatty_deg, oil_KO_name |> filter(grepl(fatty_deg_EC_list[i], Name)) |> inner_join(oil_KO_tpm) |> select(ECHO_Blank, ECHO_C, ECHO_F) |> as.matrix() |> colSums() |> as.matrix())
  }
  colnames(fatty_deg)[i] <- fatty_deg_EC_list[i]
}

fatty_deg_sum <- t(fatty_deg)
fatty_deg_sum <- cbind(fatty_deg_sum, fatty_deg_EC_list)[,c(4,1,2,3)]
rownames(fatty_deg_sum)<-1:length(fatty_deg_EC_list)
colnames(fatty_deg_sum)[1] <- "Name"
str(fatty_deg_sum)
class(fatty_deg_sum)
fatty_deg_sum<- as.data.frame(fatty_deg_sum)
# class(fatty_deg_sum_a)
# str(fatty_deg_sum_a)
# fatty_deg_sum[c(2,3,4)] <- sapply(fatty_deg_sum[c(2,3,4)], as.numeric)
# str(fatty_deg_sum_a)
fatty_deg_sum <- filter(fatty_deg_sum, ECHO_F!=0)

library(reshape2)
# fatty_deg_sum <- as.data.frame(fatty_deg_sum)
fatty_deg_sum_melt <- melt(fatty_deg_sum, id="Name")
str(fatty_deg_sum_melt)
fatty_deg_sum_melt[,3] <- as.double(fatty_deg_sum_melt[,3])
fatty_deg_sum_melt$variable <- factor(fatty_deg_sum_melt$variable, levels = c("ECHO_F","ECHO_C","ECHO_Blank"))


library(ggplot2)
ggplot(fatty_deg_sum_melt, aes(x=Name, y=value, fill=variable))+
  geom_bar(stat = "identity",  colour="black") +
  scale_x_discrete(limits=rev(fatty_deg_sum[,1])) +
  coord_flip() +
  scale_fill_discrete(breaks=c("ECHO_Blank","ECHO_C","ECHO_F"), labels=c("Blank", "C_P8.9", "F_N33")) +
  labs(fill="", x="", y="TPM")


ggplot(fatty_deg_sum_melt, aes(x=variable, y=Name, fill=value))+
  geom_tile(colour="black")+
  # scale_y_discrete(limits=rev(fatty_deg_sum[,1]))+
  scale_fill_gradient(low = "#F9E79F", high = "#E67E22") +
  labs(fill="TPM", x="", y="")+
  scale_x_discrete(labels=c("ECHO_Blank"="Blank","ECHO_C"="C_P8.9","ECHO_F"="F_N33"))


#-------------normalization

class(fatty_deg_sum[,c(2,3,4)])
as.numeric(apply())

fatty_deg_sum_nor <- apply(fatty_deg_sum[,c(2,3,4)], 2, function(x) as.numeric(as.character(x)))
fatty_deg_sum_nor <-  as.data.frame(fatty_deg_sum_nor)
class(fatty_deg_sum_nor)
fatty_deg_sum_nor <- cbind(fatty_deg_sum_nor, fatty_deg_sum[,1])[,c(4,1,2,3)]
colnames(fatty_deg_sum_nor)[1] <- "Name"
fatty_deg_sum_nor
str(fatty_deg_sum_nor)


fatty_deg_sum_nor[1,2]/(fatty_deg_sum_nor[1,2]+fatty_deg_sum_nor[1,3]+fatty_deg_sum_nor[1,4])
fatty_deg_sum_nor[1,3]/(fatty_deg_sum_nor[1,2]+fatty_deg_sum_nor[1,3]+fatty_deg_sum_nor[1,4])
fatty_deg_sum_nor[1,4]/(fatty_deg_sum_nor[1,2]+fatty_deg_sum_nor[1,3]+fatty_deg_sum_nor[1,4])

fatty_deg_sum_normal <- fatty_deg_sum_nor 

for (i in 1:length(fatty_deg_sum_nor[,1])) {
  print(i)
  fatty_deg_sum_normal[i,2] <- fatty_deg_sum_nor[i,2]/(fatty_deg_sum_nor[i,2]+fatty_deg_sum_nor[i,3]+fatty_deg_sum_nor[i,4])
  fatty_deg_sum_normal[i,3] <- fatty_deg_sum_nor[i,3]/(fatty_deg_sum_nor[i,2]+fatty_deg_sum_nor[i,3]+fatty_deg_sum_nor[i,4])
  fatty_deg_sum_normal[i,4] <- fatty_deg_sum_nor[i,4]/(fatty_deg_sum_nor[i,2]+fatty_deg_sum_nor[i,3]+fatty_deg_sum_nor[i,4])
}

fatty_deg_nor_melt <- melt(fatty_deg_sum_normal, id="Name")
str(fatty_deg_nor_melt)
fatty_deg_nor_melt$variable <- factor(fatty_deg_nor_melt$variable, levels = c("ECHO_Blank","ECHO_C","ECHO_F"))


ggplot(fatty_deg_nor_melt, aes(x=variable, y=Name, fill=value))+
  geom_tile(colour="black")+
  scale_fill_gradient(low = "#F9E79F", high = "#E67E22") +
  labs(fill="TPM", x="", y="")+
  scale_x_discrete(labels=c("ECHO_Blank"="Blank","ECHO_C"="C_P8.9","ECHO_F"="F_N33"))
