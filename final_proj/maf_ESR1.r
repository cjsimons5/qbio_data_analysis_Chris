if(!requireNamespace("maftools"))BiocManager::install(c("maftools"))
library(maftools)
library(TCGAbiolinks)

mutation <- GDCquery_Maf(tumor = "BRCA",save.csv=TRUE, pipeline="mutect2")
maf_dataframe = read.maf(mutation)

jpeg("LolliPlot_ALL.jpeg")
lollipopPlot(maf=maf_dataframe, gene="ESR1")