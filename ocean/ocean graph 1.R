install.packages("reshape2")
install.packages("grid")
library(ggplot2)
library(SQMtools)
library(grid)

load(file = "./ocean_trim.RData")

rm(image, image2, image3, image4, image5, image6, kegg_1, kegg_2, pfam_1,cog_1)
rm(crtB,crtP,crtQ,ZDS_crtQ)
rm(sum_16_25depth, sum_16depth, sum_25depth)
#pfam_tpm<- Ocean$functions$PFAM$tpm

plotTaxonomy(Ocean, rank='class', count='percent')

ocean_carotenoid <- Ocean$orfs$table[grep('crt[a-z]|LUT5|GGPS$',ignore.case=T, Ocean$orfs$table$'Gene name'),]

ocean_orfs <- Ocean$orfs$table
write.csv(ocean_orfs, file = "./ocean_orfs.csv")

oceanK10027 <- Ocean$orfs$table[grep('K10027',ignore.case=T, Ocean$orfs$table$'KEGG ID'),]
oceanK15746 <- Ocean$orfs$table[grep('K15746',ignore.case=T, Ocean$orfs$table$'KEGG ID'),]

str(ocean_carotenoid)

colnames(ocean_carotenoid) <- c("Contig_ID", "Molecule", "Method",
                                "Length_NT", "Length_AA", "GC_perc", "Gene_name", "Tax", "KEGG_ID",
                                "KEGGFUN", "KEGGPATH", "COG_ID", "COGFUN", "COGPATH", "PFAM", "TPM_16depth",
                                "TPM_25depth", "Coverage_16depth", "Coverage_25depth", "Raw_read_count_16depth",
                                "Raw_read_count_25depth", "Raw_base_count_16depth", "Raw_base_count_25depth", "Hits")


# crtQ <- subset(ocean_carotenoid, Gene_name=="crtQ")
# ZDS_crtQ <- subset(ocean_carotenoid, Gene_name=="ZDS, crtQ")
# crtB <- subset(ocean_carotenoid, Gene_name=="crtB")
# crtP <- subset(ocean_carotenoid, Gene_name=="crtP")

install.packages("plyr")
library(plyr)

#sum_tpm_16depth<- ddply(ocean_carotenoid, .(Gene_name, Tax), summarise, sum_16=sum(TPM_16depth))
#sum_tpm_25depth<- ddply(ocean_carotenoid, .(Gene_name, Tax), summarise, sum_25=sum(TPM_25depth))

#sum_tpm_16_25depth <- merge(x=sum_tpm_16depth, y=sum_tpm_25depth, by=c('Gene_name','Tax'))

sum_read_16depth<- ddply(ocean_carotenoid, .(Gene_name, Tax), summarise, sum_16depth=sum(Raw_read_count_16depth))
sum_read_25depth<- ddply(ocean_carotenoid, .(Gene_name, Tax), summarise, sum_25depth=sum(Raw_read_count_25depth))
sum_read_16_25depth <- merge(x=sum_read_16depth, y=sum_read_25depth, by=c('Gene_name','Tax'))

write.csv(sum_read_16_25depth, "sum_read_16_25depth.csv")

fixed_sum_reads <- read.csv(file = 'fixed_sum_read_16_25depth.csv', )

function()

library(tidyverse)
library(ggplot2)

sum_read_16_25depth$Gene_name <- 
  factor(sum_read_16_25depth$Gene_name, levels=rev(c("GGPS", "crtB", "crtI", "lcyB, crtL1, crtY",
                                                  "crtD", "crtF", "crtISO, crtH", "crtO", "crtP",
                                                  "crtQ", "crtR", "crtU, cruE", "crtW, BKT", "crtX",
                                                  "LUT5, CYP97A3", "PDS, crtP", "ZDS, crtQ")))


depth16 <- ggplot(sum_read_16_25depth, aes(x= Gene_name, y=sum_16depth, fill=Tax)) +
  geom_bar(stat = "identity") +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank()
        ) +
  coord_flip()
  #scale_x_continuous(breaks = c(500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000), labels = c("500","1000","1500","2000","2500","3000","3500","4000","4500","5000","5500","6000"))
  # labs(x="")
depth16

library(scales)

depth25 <- 
  ggplot(sum_read_16_25depth, aes(x= Gene_name, y=sum_25depth, fill=Tax)) +
  geom_bar(stat = "identity") +
  theme(legend.position = "none",
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(hjust=0.5)
        ) +
  coord_flip()+
  labs(x="")+
  scale_y_continuous(breaks = seq(0,2000,2000))

depth25

grid.newpage()
print(depth25, vp = viewport(x = 0.8, y = 0.5, width = 0.3, height = 1.0))
print(depth16, vp = viewport(x = 0.35, y = 0.5, width = -0.6, height = 1.0))

ggplot(sum_read_16_25depth, aes(x=sum_16depth/(sum_25depth + sum_16depth), 
                                y=sum_25depth/(sum_25depth + sum_16depth), color=Tax))+
  geom_point()+
  theme(legend.position= "none")

sum_read_16_25depth$relative_RC_16depth <- sum_read_16_25depth$sum_16depth/sum(sum_read_16_25depth$sum_16depth)*100
sum_read_16_25depth$relative_RC_25depth <- sum_read_16_25depth$sum_25depth/sum(sum_read_16_25depth$sum_16depth)*100

sum_read_16_25depth |>
  group_by(Gene_name) |>
  summarise(relative_RC_16depth = sum(relative_RC_16depth)) |>
ggplot(aes(x=Gene_name, y=relative_RC_16depth))+
  geom_point()

#==============

depth16 <- 
  sum_read_16_25depth |>
  group_by(Gene_name) |>
  summarise(sum_16depth=sum(sum_16depth)) |>
  ggplot(aes(x= Gene_name, y=sum_16depth)) +
  geom_bar(stat = "identity") +
  scale_y_log10() +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank()
  ) +
  coord_flip()

#scale_x_continuous(breaks = c(500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000), labels = c("500","1000","1500","2000","2500","3000","3500","4000","4500","5000","5500","6000"))
# labs(x="")
depth16

library(scales)

depth25 <-
  sum_read_16_25depth |>
  group_by(Gene_name) |>
  summarise(sum_25depth=sum(sum_25depth)) |>
  ggplot(aes(x= Gene_name, y=sum_25depth)) +
  geom_bar(stat = "identity") +
  theme(legend.position = "none",
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(hjust=0.5)
  ) +
  coord_flip()+
  labs(x="")+
  scale_y_log10()
depth25

grid.newpage()
print(depth25, vp = viewport(x = 0.8, y = 0.5, width = 0.3, height = 1.0))
print(depth16, vp = viewport(x = 0.35, y = 0.5, width = -0.6, height = 1.0))

library(reshape2)

sum_read_16_25depth2 <- melt(sum_read_16_25depth[,c("Gene_name", "sum_16depth", "sum_25depth")])
colnames(sum_read_16_25depth2)

ggplot(sum_read_16_25depth2, aes(x= variable, y= value, fill=variable))+
  geom_bar(stat="identity")+
  facet_wrap(Gene_name~., scales = "free", ncol=6)+
  labs(fill="")+
  annotate(geom = "segment", x = c(-Inf, -Inf), y = c(-Inf, -Inf), xend=c(-Inf, Inf), yend = c(Inf, -Inf))+
  theme(legend.position = "top",
        axis.text.x = element_text( angle = 30 , hjust=1),
        panel.background = element_blank())
