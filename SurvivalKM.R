library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(DT)
library(MMRFBiolinks)


#downloading clinical data
clin.mm<-MMRFqueryGDC_clinic(type = "clinical")

# grouping samples by therapy
bor.samples<-MMRFgetGDC_BarcodeTherapy("Bortezomib",clin.mm)
dexa.samples<- MMRFgetGDC_BarcodeTherapy("Dexamethasone",clin.mm)

#subsetting samples groups to make faster the run

bor.samples<-bor.samples[1:8]
dexa.samples<-dexa.samples[1:8]

# selecting for each group the samples that not are included in the other group

bor.samples<- bor.samples[!bor.samples %in% dexa.samples]
dexa.samples<- dexa.samples[!dexa.samples %in% bor.samples]


listSamples<-c(bor.samples,dexa.samples)

query.mm <- GDCquery(project = "MMRF-COMMPASS", 
                     data.category = "Transcriptome Profiling",
                     data.type = "Gene Expression Quantification",
                     experimental.strategy = "RNA-Seq",
                     workflow.type="HTSeq - FPKM",
                     barcode = listSamples)


#downloading samples
GDCdownload(query.mm)

# GDCdownload(query.mm, method = "api", files.per.chunk = 100)

# Prepare expression matrix with geneID in the rows and samples (barcode) in the columns

MMRFdata.prep <- MMRFGDC_prepare(query.mm,
                                 save = TRUE ,
                                 save.filename = "RNASeqSE.rda" ,
                                 directory = "GDCdata",
                                 summarizedExperiment = TRUE)


MMRFdataPrepro <- MMRFanalyzeGDC_Preprocessing(object = MMRFdata.prep,
                                               cor.cut = 0,
                                               datatype = "HTSeq - FPKM",
                                               filename ="MMRF_Preprocessing.png")



# extract the substring of the sample identifier in MMRFdataPrepro.log that map with sample identifier in clin.mm

colnames(MMRFdataPrepro) <- substr(colnames(MMRFdataPrepro),1,9)


G_list<-rownames(MMRFdataPrepro)

# subset of genes list to make faster the run
G_list<-G_list[1:100]





gr1<-bor.samples[bor.samples %in% colnames(MMRFdataPrepro)]
gr2<-dexa.samples[dexa.samples %in% colnames(MMRFdataPrepro)]


tabSurvKM <- MMRFanalyzeGDC_SurvivalKM(clin.mm,
                                      MMRFdataPrepro,
                                       Genelist = G_list,
                                       Survresult = TRUE,
                                       p.cut = 0.6,
                                       ThreshTop = 0.67,
                                       ThreshDown = 0.33,
                                       group1=gr1, 
                                       group2=gr2)


tabSurvKM <- tabSurvKM[order(tabSurvKM$pvalue, decreasing=F),]


col.names <- c("pvalue","Group1 Deaths","Group1 Deaths with Top","Group1 Deaths with Down", "Mean Group1 Top", "Mean Group1 Down","Mean Group2")
colnames(tabSurvKM) <- col.names

tabSurvKM<-tabSurvKM[1:10,]

