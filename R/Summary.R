library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(DT)
library(MMRFBiolinks)
library(EnsDb.Hsapiens.v79)
library(jsonlite)
library(RCurl)
library(httr)
library(XML)


summary<-TCGAbiolinks:::getProjectSummary("MMRF-COMMPASS")
summary.df<-summary$data_categories

summary.dataCategory <- getDataCategorySummary("MMRF-COMMPASS")

#--------------------------------------------------
project<-"MMRF-COMMPASS"

getDataCategorySummary <- function(project){
  baseURL <- "https://api.gdc.cancer.gov/files/?"
  url <- paste0(baseURL,"&expand=cases&size=100000&fields=cases.submitter_id,data_category&filters=",
                URLencode('{"op":"and","content":[{"op":"in","content":{"field":"cases.project.project_id","value":["'),
                URLencode(project),
                URLencode('"]}}]}'))
  
  url <- paste0(url, URLencode(']}'))
  
  json<-getURL(url,fromJSON,timeout(600),simplifyDataFrame = TRUE)
  
  
  resp<-GET(url)
  resp2<-data <- xmlParse(resp)
  
  jsonRespText<-content(resp,as="text")
  jsonRespParsed<-content(resp,as="parsed")
  fromJSON(jsonRespText)
  json<-fromJSON(jsonRespText)
  
  
  
  json.df<-json$data$summary$file_count
  json.df<-json$data$summary$experimental_strategies
  json.df<-json$data$summary$data_categories
  json.df<-json$data$summary$case_count 
  
  
  
  
  
  
  
  
  
  json  <- tryCatch(
    getURL(url,fromJSON,timeout(600),simplifyDataFrame = TRUE),
    error = function(e) {
      fromJSON(content(getURL(url,GET,timeout(600)), as = "text", encoding = "UTF-8"), simplifyDataFrame = TRUE)
    }
  )
  json <- json$data$hits
  json$submitter_id <- unlist(lapply(json$cases, function(x) paste0(x$submitter_id,collapse = ",")))
  json$cases <- NULL
  json <- json[!duplicated(json),]
  json <- json[stringr::str_length(json$submitter_id) == 12,]
  ret <- as.data.frame.matrix(xtabs(~ submitter_id + data_category , json))
  return(ret)
}
#------------------------------------------------------

query.fpkm.platform.tp <- GDCquery(project = "MMRF-COMMPASS", 
                                data.category = "Transcriptome Profiling")


query.fpkm.platform.snv <- GDCquery(project = "MMRF-COMMPASS", 
                                   data.category = "Simple Nucleotide Variation")


query.fpkm.platform.sr <- GDCquery(project = "MMRF-COMMPASS", 
                                    data.category = "Sequencing Reads")



summary.tp<-MMRFqueryGDC_Summary(query.fpkm.platform.tp)
summary.snv<-MMRFqueryGDC_Summary(query.fpkm.platform.snv)
summary.sr<-MMRFqueryGDC_Summary(query.fpkm.platform.sr)





