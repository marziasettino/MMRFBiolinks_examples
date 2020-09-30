library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(DT)
library(MMRFBiolinks)


# IMPORTANT 1! If you want to reproduce exaclty the case study 2 by using complete MMRF Researcher Gateway (RG) clinical data, 
# you should gain access to MMRF Researcher Gateway (RG) for downloading clinical data (MMRF_CoMMpass_IA14_STAND_ALONE_TRTRESP.csv at the time of this script). 
# Once this is done, you should import them into your R environmnet for feeding the following Workflow.
# For feeding the workflow with downloaded data, you should leave commented <get(load("data/clinMMGateway.rda"))> and <listSamples>.

#Note1: pay attention to import MMRF-RG dataset into R environment including the heading of columns.


# IMPORTANT 2! If you want initially just test the proposed MMRFBiolinks functions without using MMRF Researcher Gateway (RG) clinical data,
# you can use fictitious data that we provide you in <clinMMGateway> and yiu should uncomment <listSamples>.



mmrfListSamples<- c("MMRF_1014","MMRF_1017",
                    "MMRF_1024","MMRF_1038",
                    "MMRF_1033","MMRF_1007",
                    "MMRF_1052","MMRF_1082",
                    "MMRF_1094","MMRF_1093",
                    "MMRF_1129","MMRF_1130")



#listSamples <- c("MMRF_0001","MMRF_0002",
#                 "MMRF_0003","MMRF_0004",
#                 "MMRF_0005","MMRF_0006",
#                 "MMRF_0007","MMRF_0008",
#                 "MMRF_0009","MMRF_0010",
#                 "MMRF_0011","MMRF_0012")



clinMMGateway<- read.csv("~/MMRF_CoMMpass_IA14_STAND_ALONE_TRTRESP.csv")
#clinMMGateway <- get(load("data/clinMMGateway.rda")) 

#filter clinMMGateway by samples identifiers
bestOverall<-MMRFGetGateway_BOresponse(mmrfListSamples, clinMMGateway)
#bestOverall<-MMRFGetGateway_BOresponse(listSamples, clinMMGateway)
bestOveall<-unique(bestOverall)



#filter clinical data by Best Overall Response (BOR) type
bestOveallType<-MMRFGetGateway_BOresponseType(clinMMGateway,"CR")  

# Draw plot of the Best Overall Response to the Treatment: only the subset of samples filtered by therapyname="Bortezomib" is considered.
bestOveallPlot1<-MMRFGetGateway_BOresponsePlot(clinMMGateway,"Bortezomib",height=4, width=8, filename = "img/BestOverallPlot")

#Draw plot of the Best Overall Response to the Treatment.
bestOveallPlot2<-MMRFGetGateway_BOresponsePlot(clinMMGateway,height=4, width=8,filename = "img/BestOverallPlot_2")



#Draw plot of Time (days) Vs the Best Overall Response
MMRFGetGateway_TimeBOresponsePlot(clinMMGateway,"Dexamethasone","days", filename="img/TimeBestOverall_responsePlot",height=10, width=15)



#Draw plot of the Treatment duration (cycle)
MMRFGetGateway_TrtBOdurationPlot(clinMMGateway,"Bortezomib",ttime="days",bor="PR",filename="img/Trt_DurationPlot",height=10, width=10)


