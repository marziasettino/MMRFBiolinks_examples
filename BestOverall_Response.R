library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(DT)
library(MMRFBiolinks)


# IMPORTANT 1! Data used in this example script are fictitious because of registration required by 
# MMRF Researcher Gateway (RG) for accessing data.
# Fictitious data can be loaded by using clinMMGateway <- get(load("clinMMGateway")). 


# IMPORTANT 1! If you want to reproduce the case study 2 by using complete MMRF Researcher Gateway (RG) clinical data, 
#you would need to download data from MMRF Researcher Gateway (RG) and import them into your R environmnet for feeding the following Workflow.
# Before feeding the workflow with downloaded data, you would need to comment <get(load("data/clinMMGateway.rda"))>, replace <listSamples> with
# your list of samples. Once this is done, you would can import downloaded data into R environmnet for using them instead of 
#<clinMMGateway>.


# Note1:at the time of this script, the response-treatment file name to download from MMRF Researcher Gateway (RG) 
# is "MMRF_CoMMpass_IA14_STAND_ALONE_TRTRESP.csv".

#Note2: pay attention to import MMRF-RG dataset into R environment including the heading of columns.




listSamples <- c("MMRF_0001","MMRF_0002",
                 "MMRF_0003","MMRF_0004",
                 "MMRF_0005","MMRF_0006",
                 "MMRF_0007","MMRF_0008",
                 "MMRF_0009","MMRF_0010",
                 "MMRF_0011","MMRF_0012",
                 "MMRF_0013")




clinMMGateway <- get(load("data/clinMMGateway.rda")) 

#filter clinMMGateway by samples identifiers
bestOveall<-MMRFgetGateway_BOresponse(listSamples, clinMMGateway)




#bar.dexa<-MMRFGetGDC_BarcodeTherapy("Dexamethasone",clin.mm)
#MMRFGetGDC_Treatments(clin.mm)
#bestOveallType<-MMRFGetGateway_BOresponseType(clinMMGateway,"PR" ) 
#MMRFget_InfoCohort(query.mm)
#MMRFget_NCasesCohort(query.mm)
#Convert_toGeneSymbol(ensembl.genes)
#Convert_toGeneEnsembl(symbol.gene)


# Draw plot of the Best Overall Response to the Treatment: only the subset of samples filtered by therapyname="Bortezomib" is considered.

bestOveallPlot1<-MMRFGetGateway_BOresponsePlot(clinMMGateway,"Bortezomib",height=5, width=8, filename = "img/BestOverallPlot")

#Draw plot of the Best Overall Response to the Treatment.
bestOveallPlot2<-MMRFGetGateway_BOresponsePlot(clinMMGateway,topN=40, height=15, width=15,filename = "img/BestOverallPlot_2")



#Draw plot of Time (days) Vs the Best Overall Response
MMRFGetGateway_TimeBOresponsePlot(clinMMGateway,"Dexamethasone","days", filename="img/TimeBestOverall_responsePlot")



#Draw plot of the Treatment duration (cycle)
MMRFgetGateway_TrtBOduration(clinMMGateway,"Bortezomib",ttime="cycles",bor="PR",filename="img/Trt_DurationPlot",height=10, width=10)


