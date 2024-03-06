library(readxl)
library(tidyverse)
library(rjson)
library(phenoptr)
library(phenoptrReports)
library(ggplot2)
library(colorRamp2)

source("maps_plot_functions.R")
source("logger.R")

maps=function(){
  log_name=paste0("maps_plot_",format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))
  logger(log_name)
  cat("\n")
  log4r_info("Start mapping plot: This step will plot all phenotypes positive cells coordinates")
  
  # Load json file
  log4r_info("Loading configuration file...")
  
  tryCatch(
    {
      myData=fromJSON(file="config.json")
    }, error=function(cond){
      log4r_error(paste0("Something went wrong while reading config.json file. Try check the fields of your json file!"))
      log4r_error("Exiting from clean_data.R script")
      return(file.path(myData$Paths["output_folder"][[1]],"Log",paste0(log_name,".log")))
    })
  
  if (file.exists(file.path(myData$Paths["output_folder"][[1]],"Merged_clean"))){
    log4r_info("Merged_clean folder found! Data will be taken by this folder")
    output_f=file.path(myData$Paths["output_folder"][[1]],"Merged_clean")
  } else {
    log4r_info("Only Merged folder found! Data will be taken by this folder")
    output_f=file.path(myData$Paths["output_folder"][[1]],"Merged")
  }
  cat("\n")
  if (length(myData$Map_plot["pheno_list"][[1]])==0){
    log4r_info("The plot of all phenotypes was selected!")
    pheno_list=myData$Phenotypes["pheno_list"][[1]]
    log4r_info(paste0("The following phenotypes will be printed: ",paste(pheno_list, collapse = '; ')))
  }else{
    log4r_info("The plot of a sublist of phenotypes was selected!")
    pheno_list=myData$Map_plot["pheno_list"][[1]]
    log4r_info(paste0("The following phenotypes will be printed: ",paste(pheno_list, collapse = '; ')))
  }
  
  groups=list.files(output_f)
  log4r_info(paste0(length(groups)," groups found!"))
  
  #plot generation
  for (group in groups){
    out_f=file.path(myData$Paths["output_folder"],"Maps_plot")
    dir.create(out_f)
    
    out_f_gr=file.path(out_f,group)
    dir.create(out_f_gr)
    
    for (file in list.files(file.path(output_f,group), full.names = T)){
      print("################################")
      cat("\n")
      print(paste0("Analyzing ",file, "..."))
      id=gsub(".txt","",strsplit(file,"_")[[1]][length(strsplit(file,"_")[[1]])])
      data<-read.delim(file)
      
      print("starting multi images plot")
      
      filt_data=filt_data_Pheno(data, id, pheno_list)
      if (length(filt_data$Pheno)==0){
        log4r_warn("No phenotype(s) found for subject! Skip to the next one")
        next
      }
      multi_maps_plot(filt_data, id ,out_f_gr, myData$Map_plot["interactive"][[1]])
  
    }
    log4r_info("All done!")
    
    return(file.path(myData$Paths["output_folder"][[1]],"Log",paste0(log_name,".log")))
  }
}
