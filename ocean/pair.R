install.packages("devtools")
devtools::install_github("dwinter/pafr")

getwd()
rm(test_alignment)

library(pafr, quietly=TRUE)
test_alignment <- system.file("extdata", "ragtag.paf", package="pafr")
ali <- read_paf(test_alignment)
dotplot(ali)
