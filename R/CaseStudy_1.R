library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(DT)
library(MMRFBiolinks)
library(EnsDb.Hsapiens.v79)

#Get MMRF-COMMPASS Project Summary

summary<-MMRFgetProjectSummary()
summary$data_categories
summary$experimental_strategies


#downloading clinical data
clin.mm<-MMRFqueryGDC_clinic(type = "clinical")

# grouping samples by therapy
bor.samples<-MMRFGetGDC_IdentifierByTherapy("Bortezomib",clin.mm)
dexa.samples<- MMRFGetGDC_IdentifierByTherapy("Dexamethasone",clin.mm)



#subsetting samples groups to make faster the run

bor.samples<-bor.samples[1:15]
dexa.samples<-dexa.samples[1:15]

# selecting for each group the samples that not are included in the other group

bor.samples<- bor.samples[!bor.samples %in% dexa.samples]
dexa.samples<- dexa.samples[!dexa.samples %in% bor.samples]


listSamples<-c(bor.samples,dexa.samples)

query.fpkm <- GDCquery(project = "MMRF-COMMPASS", 
                     data.category = "Transcriptome Profiling",
                     data.type = "Gene Expression Quantification",
                     experimental.strategy = "RNA-Seq",
                     workflow.type="HTSeq - FPKM",
                     barcode = listSamples)




summary<-MMRFqueryGDC_Summary(query.fpkm)



#downloading samples
GDCdownload(query.fpkm)

# GDCdownload(query.fpkm, method = "api", files.per.chunk = 100)

# Prepare expression matrix with geneID in the rows and samples (barcode) in the columns

MMRFdata.prep <- MMRFGDC_prepare(query.fpkm,
                                 save = TRUE ,
                                 save.filename = "data/RNASeqSE.rda" ,
                                 directory = "GDCdata",
                                 summarizedExperiment = TRUE)






MMRFdataPrepro <- TCGAanalyze_Preprocessing(object = MMRFdata.prep,
                                            cor.cut = 0,
                                            datatype = "HTSeq - FPKM",
                                            filename ="img/MMRF_Preprocessing.png",width = 900, height = 900)


G_list<-rownames(MMRFdataPrepro)
symbol.gene <- ensembldb::select(EnsDb.Hsapiens.v79, keys= G_list, keytype = "GENEID", columns = c("SYMBOL","GENEID"))

row.names(MMRFdataPrepro)<-symbol.gene$SYMBOL
G_list<-symbol.gene$SYMBOL

# subset of genes list to make faster the run
G_list<-G_list[1:1000]




gr1<-bor.samples[bor.samples %in% colnames(MMRFdataPrepro)]
gr2<-dexa.samples[dexa.samples %in% colnames(MMRFdataPrepro)]



dataMMcomplete <- log2(MMRFdataPrepro[1:56457,] + 1)

tabSurvKM <- TCGAanalyze_SurvivalKM(clin.mm,
                                    dataMMcomplete,
                                    Genelist = G_list,
                                    Survresult = TRUE,
                                    p.cut = 0.05,
                                    ThreshTop = 0.67,
                                    ThreshDown = 0.33,
                                    group1=gr1, 
                                    group2=gr2)











tabSurvKM <- tabSurvKM[order(tabSurvKM$pvalue, decreasing=F),]


col.names <- c("pvalue","Group1 Deaths","Group1 Deaths with Top","Group1 Deaths with Down", "Mean Group1 Top", "Mean Group1 Down","Mean Group2")
colnames(tabSurvKM) <- col.names

#tabSurvKM<-tabSurvKM[1:10,]


datatable(tabSurvKM)

#For example we select the top two genes in tabSurvKM


top.gene1<-rownames(tabSurvKM)[1]
top.gene2<-rownames(tabSurvKM)[2]



tabSurvKM.gene1 <- TCGAanalyze_SurvivalKM(clin.mm,
                                          dataMMcomplete,
                                          Genelist = top.gene1,
                                          Survresult = TRUE,
                                          p.cut = 0.05,
                                          ThreshTop = 0.67,
                                          ThreshDown = 0.33,
                                          group1=gr1, 
                                          group2=gr2)



tabSurvKM.gene2 <- TCGAanalyze_SurvivalKM(clin.mm,
                                          dataMMcomplete,
                                          Genelist = top.gene2,
                                          Survresult = TRUE,
                                          p.cut = 0.05,
                                          ThreshTop = 0.67,
                                          ThreshDown = 0.33,
                                          group1=gr1, 
                                          group2=gr2)













