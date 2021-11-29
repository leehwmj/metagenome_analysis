# EC number로 나누는 작업
## KEGG number에 따른 TPM 수치 먼저 수집
getwd() #"/analysis/hwlee/collect-temp-sequencing/HN00149294_Rawdata"
### KO = KEGG Orthology
oil_KO_tpm <- read.table(file = './HN00149294_RawData_squeezemeta_spades_0623/results/tables/HN00149294_RawData_squeezemeta_spades_0623.KO.tpm.tsv', sep = '\t', header = TRUE)
## KEGG number에 대응되는 name 수집
oil_KO_name <- read.table(file = './HN00149294_RawData_squeezemeta_spades_0623/results/tables/HN00149294_RawData_squeezemeta_spades_0623.KO.names.tsv', sep = '\t', header = TRUE, quote = "", row.names=NULL, stringsAsFactors = FALSE)
# 아래와 같은 에러 발생하여 'quote = "", row.names=NULL, stringsAsFactors = FALSE' 추가함
# Warning messages:
# 1: In scan(file = file, what = what, sep = sep, quote = quote, dec = dec,  :
#  EOF within quoted string
# 2: In scan(file = file, what = what, sep = sep, quote = quote, dec = dec,  :
#  number of items read is not a multiple of the number of columns

library(tidyverse)
# lipase, esterase, protease로 검색했을 때, EC number를 가지고 있지 않은 것들 검색


#루프 만들어서 EC 골라내기
## EC number list 만들기
EC1_list <- paste("EC:1.", c(1:23), sep = "")
EC1_list <- c(EC1_list, paste("EC:1.", c(97:99), sep = ""))


rm(a)

#EC number를 뽑아서 합계를 구하는 코드
for(i in 1:length(EC1_list)){
  print (EC1_list[i])
  if(i == 1){
    a <- oil_KO_name |> filter(grepl(EC1_list[i], Name)) |> inner_join(oil_KO_tpm) |> select(ECHO_Blank, ECHO_C, ECHO_F) |> as.matrix() |> colSums() |> as.matrix()
  }else{
    a <- cbind(a, oil_KO_name |> filter(grepl(EC1_list[i], Name)) |> inner_join(oil_KO_tpm) |> select(ECHO_Blank, ECHO_C, ECHO_F) |> as.matrix() |> colSums() |> as.matrix())
  }
  colnames(a)[i] <- EC1_list[i]
}

EC1_all_sum <- t(a)
EC1_all_sum <- cbind(EC1_all_sum, EC1_list)[,c(4,1,2,3)]
rownames(EC1_all_sum)<-1:length(EC1_list)
colnames(EC1_all_sum)[1] <- "Name"

library(reshape2)
EC1_all_sum <- as.data.frame(EC1_all_sum)
EC1_melt <- melt(EC1_all_sum, id="Name")
str(EC1_melt)
EC1_melt[,3] <- as.double(EC1_melt[,3])
EC1_melt$variable <- factor(EC1_melt$variable, levels = c("ECHO_F","ECHO_C","ECHO_Blank"))
str(EC1_melt)

#
library(ggplot2)
library(plyr)
ggplot(EC1_melt, aes(x=Name, y=value, fill=variable))+
  geom_bar(stat = "identity",  colour="black")+
  scale_x_discrete(limits=rev(EC1_all_sum[,1]))+
  coord_flip() +
  scale_fill_discrete(breaks=c("ECHO_Blank","ECHO_C","ECHO_F"), labels=c("Blank", "C_P8.9", "F_N33")) +
  labs(fill="", x="", y="TPM")

ggplot(EC1_melt, aes(x=Name, y=value, fill=variable))+
  geom_bar(stat = "identity", position = "fill", colour="black")+
  scale_x_discrete(limits=rev(EC1_all_sum[,1]))+
  coord_flip() +
  scale_fill_discrete(breaks=c("ECHO_Blank","ECHO_C","ECHO_F"), labels=c("Blank", "C_P8.9", "F_N33")) +
  labs(fill="", x="", y="Ratio")



