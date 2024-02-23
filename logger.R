library(logr)
library(rjson)

logger=function(log_name){
  
  myData=fromJSON(file="config.json")
  output_folder=myData$Paths["output_folder"][[1]]
  
  log_open(file.path(output_folder,log_name))
}

log4r_info <- function(message) {
  log_print(paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " INFO: ",message),hide_notes = T)
}

log4r_warn <- function(message) {
  log_print(paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " WARNING: ",message),hide_notes = T)
}

log4r_error <- function(message) {
  log_print(paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ERROR: ",message),hide_notes = T)
}

