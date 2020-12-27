library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(DT)
library(stringr)
library(ggplot2)
library(MMRFBiolinks)


# IMPORTANT 1! If you want to reproduce exaclty the case study 3 by using complete MMRF Researcher Gateway (RG) clinical data, 
# you should gain access to MMRF Researcher Gateway (RG) for downloading the following files:
# i)MMRF_CoMMpass_IA14a_All_Canonical_Variants.csv
#ii)MMRF_CoMMpass_IA14_PER_PATIENT.csv
#iii)MMRF_CoMMpass_IA14_STAND_ALONE_TRTRESP.csv


# Once this is done, you should import them into your R environmnet for feeding the following Workflow.
# For feeding the workflow with downloaded data, you should leave commented <get(load("data/CaseStudy3.rda"))> .

#Note1: pay attention to import MMRF-RG dataset into R environment including the heading of columns.


# IMPORTANT 2! If you want initially just test the proposed MMRFBiolinks functions without using MMRF Researcher Gateway (RG) clinical data,
# you can use fictitious data that we provide you in <clinMMGateway> and you should uncomment <listSamples>.

 
get(load("data/CaseStudy3.rda")) 
 


summary.var<-MMRF_RG_VariantCountPlot(variant.ann,trt,topN=50,filenm=NULL)


variant <- c("rs755588843", "rs569344016","rs2066497")



patient.var<-MMRF_RG_GetIDSamplebyVariant(variant.ann,patient,variant)




MMRF_RG_SurvivalKM(patient.var,
                   trt,
                   FilterBy="treatment", 
                   filename=NULL,
                   xlim = c(100,3000),
                   conf.range = FALSE,
                   color = c("Dark2"))