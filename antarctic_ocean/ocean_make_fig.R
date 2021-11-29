rm(list = ls())
getwd()
#kegg_names에서 carotenoid 관여 gene만 분리
# kegg_names <- read.table("./ocean_trim_sqmreads2tables/antarctic_ocean_SQM_reads.KO.names.tsv",
#                          sep="\t", header = TRUE, check.names = FALSE)
# kegg_names_in_carotenoid <- kegg_names[grep('Carotenoid', kegg_names$Path),]

library(tidyr)
library(tidyverse)
#kegg 사이트에서 carotenoid pathway 연관 gene 추출
carotenoid_ortholog_table <- read.csv("./antarctic_ocean/carotenoid_ortholog_table.csv",
                                       header = TRUE)
carotenoid_orth_sep <- separate(data = carotenoid_ortholog_table, col = Definition, sep = ";", into = c("gene_name", "def"))

#해수샘플에서 sqm reads로 부터 얻은 kegg abund 불러오기
kegg_abund <- read.table("./antarctic_ocean/antarctic_ocean_SQM_reads.KO.abund.tsv",
                         sep="\t", header = TRUE, check.names = FALSE)
colnames(kegg_abund)[1] <- ("kegg_name")

#inner join with kegg_abund and carotenoid_orth_sep
kegg_abund_carotenoid_ocean <- inner_join(kegg_abund, carotenoid_orth_sep, by='kegg_name')

kegg_abund_carotenoid_ocean_2 <- kegg_abund_carotenoid_ocean[,c(4,2,3)]
kegg_abund_carotenoid_ocean_2[25,1] <- "crtD(K20611)"

library(reshape2)

kegg_abund_carotenoid_ocean_melt<- melt(kegg_abund_carotenoid_ocean_2, id="gene_name")
kegg_abund_carotenoid_ocean_melt[,'sum']<- log10(apply(kegg_abund_carotenoid_ocean_2[,c(2,3)], 1, sum))

install.packages("ggthemes")
library(ggthemes)

cp <- coord_polar("y")
cp$is_free <- function() TRUE

ggplot(kegg_abund_carotenoid_ocean_melt, aes(x=sum/2, y = value, fill = variable, width = sum)) +
  geom_bar(position="fill", stat="identity", alpha=0.9) +
  facet_wrap(~ gene_name, ncol=8) +
  theme_minimal() +
  theme(panel.grid.major = element_line(colour = "grey50"), axis.text.x = element_text(size= 0)) +
  scale_fill_manual(values = c('16depth'='#21AFFB','25depth'='#1208E7')) +
  cp
  

