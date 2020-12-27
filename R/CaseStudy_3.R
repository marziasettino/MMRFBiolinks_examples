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
# For feeding the workflow with downloaded data, you should leave commented <get(load("data/CaseStudy3.rda"))> and
# <variant>.

#Note1: pay attention to import MMRF-RG dataset into R environment including the heading of columns.


# IMPORTANT 2! If you want initially just test the proposed MMRFBiolinks functions without using MMRF Researcher Gateway (RG) data,
# you can use fictitious data that we provide you in <CaseStudy3.rda> and you should uncomment <get(load("data/CaseStudy3.rda"))> and
# <variant>.


 
#get(load("data/CaseStudy3.rda")) 

#variant <- c("rs755588843", "rs569344016","rs2066497")

 
variant.ann<- read.csv("~/MMRF_CoMMpass_IA14a_All_Canonical_Variants.csv")
trt<- read.csv("~/MMRF_CoMMpass_IA14_STAND_ALONE_TRTRESP.csv")
patient<- read.csv("~/MMRF_CoMMpass_IA14_PER_PATIENT.csv")


summary.var<-MMRF_RG_VariantCountPlot(variant.ann,trt,topN=50,filenm=NULL)

#We select in summary.var the top four dbSNP ID with the high count related to "Complete Response". 

variant <- c("rs111362472", "rs181828141","rs2157615","rs2157615")

#We select in "patient" dataframe only the samples having variants in "variant".

patient.var<-MMRF_RG_GetIDSamplebyVariant(variant.ann,patient,variant)

#We perform Survival analysis focused on patients having those variants.


MMRF_RG_SurvivalKM(patient.var,
                   trt,
                   FilterBy="treatment", 
                   filename=NULL,
                   conf.range = FALSE,
                   color = c("Dark2"))