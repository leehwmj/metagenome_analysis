library()
library(SQMtools)

getwd()

load('./RDatas/antarctic_dupli_untrim.RData')

rm(list = ls())

antarctic_dupli_untrim = loadSQM('/analysis/hwlee/antartic_ocean/antarctic_ocean_dupli_untrim')

plotTaxonomy(antarctic_dupli, rank='phylum', count='percent')
plotTaxonomy(antarctic_dupli_untrim, rank='phylum', count='percent')
