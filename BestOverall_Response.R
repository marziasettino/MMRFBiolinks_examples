library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(DT)
library(MMRFBiolinks)


# IMPORTANT 1! Data used in this example script are fictitious because of registration required by MMRF Researcher Gateway (RG) for accessing data.
# Fictitious data can be loaded by using clinMMGateway <- get(load("mmrf_case2.rda")). 

# IMPORTANT 1! If you want to generalize the the case study 2, you would need to download data from MMRF Researcher Gateway (RG) and 
# import them into your R workspace for feeding the following Workflow.
# Before feeding the workflow with downloaded data, you would need to comment listSamples and get(load("mmrf_case2.rda")), 
# replace mmrf_case2.rda with with downloaded data and eventually replace listSamples with your chosen filtering list of samples.
# At the time of this script, the response-treatment file name to download from MMRF Researcher Gateway (RG) 
# is "MMRF_CoMMpass_IA14_STAND_ALONE_TRTRESP.csv".




listSamples <- c("MMRF_0001","MMRF_0002",
                 "MMRF_0003","MMRF_0004",
                 "MMRF_0005","MMRF_0006",
                 "MMRF_0007","MMRF_0008",
                 "MMRF_0009","MMRF_0010")



clinMMGateway <- get(load("mmrf_case2.rda")) 

#filter clinMMGateway by samples identifiers
bestOveall<-MMRFgetGateway_BestOverallResponse(listSamples, clinMMGateway)



# Draw plot of the Best Overall Response to the Treatment: only the subset of samples filtered by therapyname="Bortezomib" is considered.

bestOveallPlot1<-MMRFgetGateway_BestOverallResponsePlot(clinMMGateway,"Bortezomib",height=5, width=8, filename = "BestOverallPlot1")

#Draw plot of the Best Overall Response to the Treatment.
bestOveallPlot2<-MMRFgetGateway_BestOverallResponsePlot(clinMMGateway,topN=40, height=15, width=15,filename = "BestOverallPlot2")