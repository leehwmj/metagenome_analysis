rm(list=ls())
rm(.pbd_env)
rm(.Random.seed)

ls(all.names = TRUE)

getwd()

setwd("~/R project/ocean")
library('SQMtools')
# antarctic_dupli = load('antarctic_dupli.RData')
# antarctic_dupli = load('/analysis/hwlee/antartic_ocean/antarctic_dupli.RData')
# antarctic_dupli = loadSQM('/analysis/hwlee/antartic_ocean/antarctic_ocean_dupli')

plotTaxonomy(antarctic_dupli, rank = 'phylum', count='percent')

str(antarctic_dupli)

