if(!requireNamespace("maftools"))BiocManager::install(c("maftools"))
library(maftools)
library(TCGAbiolinks)
library(SummarizedExperiment)

mutation <- GDCquery_Maf(tumor = "BRCA",save.csv=TRUE, pipeline="mutect2")
maf_dataframe = read.maf(mutation)

jpeg("LolliPlot_ALL.jpeg")
lollipopPlot(maf=maf_dataframe, gene="ESR1")

query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - Counts")
GDCdownload(query) #only need this line of code ONCE to download the data
sum_exp <- GDCprepare(query)

patient_data <- colData(sum_exp)
patient_ages <- patient_data$paper_age_at_initial_pathologic_diagnosis
#here we are creating a NEW COLUMN in patient_data called "age_category"
#NOTE: This will NOT be added to colData(sum_exp). Instead it will only be added to the patient_data data table.

#create new column with age categories for all patient samples
patient_data$age_category = ifelse(patient_ages < 50, "Young", "Old")

#get shortened patient barcodes so that we can compare with
short_maf <- substr(maf_dataframe@clinical.data$Tumor_Sample_Barcode, 1,12)

#create a new column in maf_dataframe
maf_dataframe@clinical.data$short_barcodes <- short_maf

#age vector has either young or old based on age category at that barcode
maf_ages <- patient_data[short_maf, "age_category"]

#Add new column to the maf dataframe containing age info
maf_dataframe@clinical.data$Ages <- maf_ages

#Extract codes for each age group
young_codes <- maf_dataframe@clinical.data[maf_dataframe@clinical.data$Ages =="Young",]
old_codes <- maf_dataframe@clinical.data[maf_dataframe@clinical.data$Ages =="Old",]

#create maf subsets for each age group
young_maf <- subsetMaf(maf_dataframe, tsb = young_codes$Tumor_Sample_Barcode)
old_maf <- subsetMaf(maf_dataframe, tsb = old_codes$Tumor_Sample_Barcode)

jpeg("LolliPlot2OldYoung.jpeg")
lollipopPlot2(m1=young_maf, m2=old_maf, m1_name = "Young", m2_name = "Old", gene="ESR1")
?lollipopPlot2
