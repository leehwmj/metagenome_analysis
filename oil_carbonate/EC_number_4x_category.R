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
# EC:3.1로 시작하는 모두 검색
# EC31 <- filter(oil_KO_name, grepl('EC:3.1', Name))
# EC31_join <- inner_join(EC31[,c(1,2)], oil_KO_tpm)
# 
# EC31dash <- filter(EC31_join, grepl('EC:3.1.-', Name))
# EC31dash_sum <- as.matrix(apply(EC31dash[,c(-1,-2)], 2, sum))
# colnames(EC31dash_sum)[1] <- "EC:3.1.-"


for (i in 1:30) {
  print(paste("EC:4.1.", i, sep = ""))
  if (i!=30){
    
    if(i==1){
      EC31x <- EC31_join |> filter(grepl( paste("EC:3.1.", i, sep = ""), Name )) |> select(ECHO_Blank, ECHO_C, ECHO_F) |> apply(2,sum) |>
        as.matrix() |> print()
    }else{
      EC31x <- cbind(EC31x ,EC31_join |> filter(grepl( paste("EC:3.1.", i, sep = ""), Name )) |> select(ECHO_Blank, ECHO_C, ECHO_F) |> apply(2,sum) |>
        as.matrix() |> print())
    }
    colnames(EC31x)[i] <- paste("EC:3.1.", i, sep = "")
    
  }else{
    EC31x <- cbind(EC31x ,EC31_join |> filter(grepl( "EC:3.1.-", Name )) |> select(ECHO_Blank, ECHO_C, ECHO_F) |> apply(2,sum) |>
      as.matrix() |> print())
    colnames(EC31x)[i] <- "EC:3.1.-"
  }
}



# EC:3.4로 시작하는 모두 검색
EC34 <- filter(oil_KO_name, grepl('EC:3.4', Name))
EC34_join <- inner_join(EC34[,c(1,2)], oil_KO_tpm)

# EC31dash <- filter(EC31_join, grepl('EC:3.4.-', Name))
# EC31dash_sum <- as.matrix(apply(EC31dash[,c(-1,-2)], 2, sum))
# colnames(EC31dash_sum)[1] <- "EC:3.1.-"


for (i in 1:30) {
  print(paste("EC:3.4.", i, sep = ""))
  if (i!=30){
    
    if(i==1){
      EC34x <- EC34_join |> filter(grepl( paste("EC:3.4.", i, sep = ""), Name )) |> select(ECHO_Blank, ECHO_C, ECHO_F) |> apply(2,sum) |>
        as.matrix() |> print()
    }else{
      EC34x <- cbind(EC34x ,EC34_join |> filter(grepl( paste("EC:3.4.", i, sep = ""), Name )) |> select(ECHO_Blank, ECHO_C, ECHO_F) |> apply(2,sum) |>
                       as.matrix() |> print())
    }
    colnames(EC34x)[i] <- paste("EC:3.4.", i, sep = "")
    
  }else{
    EC34x <- cbind(EC34x ,EC34_join |> filter(grepl( "EC:3.4.-", Name )) |> select(ECHO_Blank, ECHO_C, ECHO_F) |> apply(2,sum) |>
                     as.matrix() |> print())
    colnames(EC34x)[i] <- "EC:3.4.-"
  }
}



EC31x <- t(EC31x)
EC34x <- t(EC34x)

EC31x <- cbind(EC31x, rownames(EC31x))[,c(4,1,2,3)]
rownames(EC31x)<-1:30
colnames(EC31x)[1] <- "Name"
EC31x<- as.data.frame(EC31x)
class(EC31x)
EC31x <- EC31x[!(EC31x$ECHO_Blank == 0), ]

EC34x <- cbind(EC34x, rownames(EC34x))[,c(4,1,2,3)]
rownames(EC34x)<-1:30
colnames(EC34x)[1] <- "Name"
EC34x<- as.data.frame(EC34x)
class(EC34x)
EC34x <- EC34x[!(EC34x$ECHO_Blank == 0), ]

str(EC31x)

library(reshape2)
EC31_melt <- melt(EC31x, id="Name")
str(EC31_melt)
EC31_melt$value <- as.double(EC31_melt$value)
EC31_melt$variable <- factor(EC31_melt$variable, levels = c("ECHO_Blank","ECHO_C","ECHO_F"))

EC34_melt <- melt(EC34x, id="Name")
EC34_melt$value <- as.double(EC34_melt$value)
# EC34_melt$variable <- factor(EC34_melt$variable, levels = c("ECHO_Blank","ECHO_C","ECHO_F"))


library(ggplot2)
ggplot(EC31_melt, aes(x=variable, y=Name, fill=value))+
  geom_tile(colour="black")+
  scale_y_discrete(limits=rev(EC31x[,1]))+
  scale_fill_gradient(low = "#F9E79F", high = "#E67E22") +
  labs(fill="TPM", x="", y="")+
  scale_x_discrete(labels=c("ECHO_Blank"="Blank","ECHO_C"="C_P8.9","ECHO_F"="F_N33"))

ggplot(EC34_melt, aes(x=variable, y=Name, fill=value))+
  geom_tile(colour="black")+
  scale_y_discrete(limits=rev(EC34x[,1]))+
  scale_fill_gradient(low = "#F9E79F", high = "#E67E22")+
  labs(fill="TPM", x="", y="")+
  scale_x_discrete(labels=c("ECHO_Blank"="Blank","ECHO_C"="C_P8.9","ECHO_F"="F_N33"))
