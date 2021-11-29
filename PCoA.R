# install.packages("tidyverse")
getwd()
setwd("/home/hwlee/analysis/oil")
library(tidyverse)
oil_orf <- read.csv("oil_orf.csv")
head(oil_orf)

oil_orf$Gene.name
sum(oil_orf$Coverage.ECHO_Blank)
sum(oil_orf$Coverage.ECHO_C)
sum(oil_orf$Coverage.ECHO_F)

is.numeric(oil_orf$Coverage.ECHO_C)

oil_orf2 <- 
  
oil_orf%>%
  group_by(Gene.name, Tax)%>%
  summarise(Coverage.ECHO_Blank = sum(Coverage.ECHO_Blank), 
            Coverage.ECHO_C = sum(Coverage.ECHO_C), 
            Coverage.ECHO_F = sum(Coverage.ECHO_F))
oil_orf3 <- 
data.frame(
gene_tax = paste(oil_orf2$Gene.name, oil_orf2$Tax, sep = "%"),
coverage_Blank = oil_orf2$Coverage.ECHO_Blank,
coverage_C = oil_orf2$Coverage.ECHO_C,
coverage_F = oil_orf2$Coverage.ECHO_F)


rownames(oil_orf3) <- oil_orf3[, 1]

oil_orf_mat <- t(oil_orf3[,-1])

oil_orf_mat

require(vegan)

dat.dist <- vegdist(oil_orf_mat
  , method= "bray")

dat.dist
# ?cmdscale
dat.pcoa <- cmdscale(dat.dist, eig=TRUE)
dat.pcoa.points <- data.frame(dat.pcoa$points)
dat.pcoa_weight <- c(dat.pcoa$eig[1]/sum(dat.pcoa$eig),
                     dat.pcoa$eig[2]/sum(dat.pcoa$eig))
dat.pcoa_weight

colnames(dat.pcoa.points) <- c("PCoA1", "PCoA2")

dat.pcoa.points$sample_label <- c("Blank", "C", "F")

ggplot(dat.pcoa.points, aes(x=PCoA1, y=PCoA2))+
  geom_point()+
  geom_text(aes(label=sample_label, vjust=2))+
  coord_fixed()+
  xlab(paste0("PCoA1 (", round(dat.pcoa_weight[1]*100, 2), " %)"))+
  ylab(paste0("PCoA2 (", round(dat.pcoa_weight[2]*100, 2), " %)"))+
  ggtitle("PCoA based on Bray-Curtis distance")

Tax <- str_split(oil_orf2$Tax, ";")

kingdom <- c() 
for(i in 1:length(Tax)){
  annotate <- if(Tax[[i]][1] %in% c("", NA)){
    "Not annotated"
  }else{
    Tax[[i]][1]
  }
  kingdom <- c(kingdom, annotate)
}

phylum <- c() 
for(i in 1:length(Tax)){
  annotate <- if(Tax[[i]][2] %in% c("", NA)){
    "Not annotated"
  }else{
    Tax[[i]][2]
  }
  phylum <- c(phylum, annotate)
}

class <- c() 
for(i in 1:length(Tax)){
  annotate <- if(Tax[[i]][3] %in% c("", NA)){
    "Not annotated"
  }else{
    Tax[[i]][3]
  }
  class <- c(class, annotate)
}

order <- c()
for(i in 1:length(Tax)){
  annotate <- if(Tax[[i]][4] %in% c("", NA)){
    "Not annotated"
  }else{
    Tax[[i]][4]
  }
  order <- c(order, annotate)
}

family <- c()
for(i in 1:length(Tax)){
  annotate <- if(Tax[[i]][5] %in% c("", NA)){
    "Not annotated"
  }else{
    Tax[[i]][5]
  }
  order <- c(order, annotate)
}

genus <- c()
for(i in 1:length(Tax)){
  annotate <- if(Tax[[i]][6] %in% c("", NA)){
    "Not annotated"
  }else{
    Tax[[i]][6]
  }
  genus <- c(genus, annotate)
}

species <- c()
for(i in 1:length(Tax)){
  annotate <- if(Tax[[i]][7] %in% c("", NA)){
    "Not annotated"
  }else{
    Tax[[i]][7]
  }
  species <- c(species, annotate)
}

