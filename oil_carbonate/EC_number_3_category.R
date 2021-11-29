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
test <- filter(oil_KO_name, grepl('lipase', Name)&!grepl('EC:3', Name))
teat4 <- filter(oil_KO_name, grepl('lipase', Name))
test2 <- filter(oil_KO_name, grepl('esterase', Name)&!grepl('EC:', Name))
test3 <- filter(oil_KO_name, grepl('protease', Name)&!grepl('EC:', Name))

#EC:3.1 골라내서 join하고 샘플별로 합치기
EC31 <- filter(oil_KO_name, grepl('EC:3.1', Name))
EC31_join <- inner_join(EC31[,c(1,2)], oil_KO_tpm)
EC31_sum <- apply(EC31_join[,c(-1,-2)], 2, sum)
View(EC31_sum)
EC31_sum_df<- as.data.frame(EC31_sum)

#EC:3.2 골라내서 join하고 샘플별로 합치기
EC32 <- filter(oil_KO_name, grepl('EC:3.2', Name))
EC32_join <- inner_join(EC32[,c(1,2)], oil_KO_tpm)
# summarise(EC32_join, sum_b=sum(ECHO_Blank), sum_c=sum(ECHO_C), sum_f=(ECHO_F))
summarise(EC32_join, sum_b=sum(ECHO_Blank))
summarise(EC32_join, sum_c=sum(ECHO_C))


# EC32_sum <- apply(EC32_join[,c(-1,-2)], 2, sum)
View(EC32_sum)
EC32_sum_df<- as.data.frame(EC32_sum)


#루프 만들어서 EC 골라내기
## EC number list 만들기
EC3_list <- paste("EC:3.", c(1:13), sep = "")
# a <- matrix(c(0,0,0), nrow = 1, ncol = 3)
# colnames(a) <- c("ECHO_Blank", "ECHO_C", "ECHO_F")

rm(a)

#EC number를 뽑아서 합계를 구하는 코드
for(i in 1:13){
  print (paste("EC:3.", i, sep = ""))
  if(i == 1){
    a <- oil_KO_name |> filter(grepl(EC3_list[i], Name)) |> inner_join(oil_KO_tpm) |> select(ECHO_Blank, ECHO_C, ECHO_F) |> as.matrix() |> colSums() |> as.matrix()
  }else{
    a <- cbind(a, oil_KO_name |> filter(grepl(EC3_list[i], Name)) |> inner_join(oil_KO_tpm) |> select(ECHO_Blank, ECHO_C, ECHO_F) |> as.matrix() |> colSums() |> as.matrix())
  }
  colnames(a)[i] <- paste("EC:3.", i, sep = "")
}

EC3_all_sum <- t(a)
EC3_all_sum <- cbind(EC3_all_sum, paste("EC:3.", 1:13, sep = ""))[,c(4,1,2,3)]
rownames(EC3_all_sum)<-1:13
colnames(EC3_all_sum)[1] <- "Name"

library(reshape2)
EC3_all_sum <- as.data.frame(EC3_all_sum)
EC3_melt <- melt(EC3_all_sum, id="Name")
str(EC3_melt)
EC3_melt[,3] <- as.double(EC3_melt[,3])
EC3_melt$variable <- factor(EC3_melt$variable, levels = c("ECHO_F","ECHO_C","ECHO_Blank"))

#
library(ggplot2)
library(plyr)
ggplot(EC3_melt, aes(x=Name, y=value, fill=variable))+
  geom_bar(stat = "identity",  colour="black")+
  scale_x_discrete(limits=rev(EC3_all_sum[,1]))+
  coord_flip() +
  scale_fill_discrete(breaks=c("ECHO_Blank","ECHO_C","ECHO_F"), labels=c("Blank", "C_P8.9", "F_N33")) +
  labs(fill="", x="", y="TPM")

ggplot(EC3_melt, aes(x=Name, y=value, fill=variable))+
  geom_bar(stat = "identity", position = "fill", colour="black")+
  scale_x_discrete(limits=rev(EC3_all_sum[,1]))+
  coord_flip() +
  scale_fill_discrete(breaks=c("ECHO_Blank","ECHO_C","ECHO_F"), labels=c("Blank", "C_P8.9", "F_N33")) +
  labs(fill="", x="", y="Ratio")



