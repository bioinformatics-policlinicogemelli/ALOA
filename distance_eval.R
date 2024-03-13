library(rjson)
library(phenoptr)
library(phenoptrReports)
library(plotly)
library(ggplot2)
library(dplyr)
library(htmlwidgets)

source("distance_match_func.R")
source("maps_plot_functions.R")
source("logger.R")

distance=function(){
  
  cat("\n################## DISTANCE ###################\n")
  
  log_name=paste0("distance_",format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))
  logger(log_name)
  cat("\n")
  log4r_info("Start distance evaluation: This step will will evaluate distances between couples of phenotypes positive cells")
  cat("\n")
  
  #load info from json
  log4r_info("Loading configuration file...")
  myData=fromJSON(file="config.json")
  
  out_f=file.path(myData$Paths["output_folder"],"Distance")
  dir.create(out_f)
  log4r_info(paste0("Distances csv file will be stored in ", out_f))
  
  output_f=file.path(myData$Paths["output_folder"],"Merged_clean")

  if (length(myData$Distance["pheno_list"][[1]])==0){
    gene_list=myData$Phenotypes["pheno_list"][[1]]
  }else{
    gene_list=myData$Distance["pheno_list"][[1]]
  }
  
  groups=list.files(output_f)
  
  #plot generation
  for (group in groups){
    cat("\n###################\n")
    log4r_info(paste0("Checking Group: ",group))
    
    #creating output folders
    out_f_gr=file.path(out_f,group)
    dir.create(out_f_gr)
    
    for (file in list.files(file.path(output_f,group), full.names = T)){
      
      log4r_info(paste0("Checking ",file))
      data<-read.delim(file)
      
      id=gsub(".txt","",gsub("Merge_cell_seg_data_clean_","",basename(file)))
      
      sub_data_cl=filt_data_Pheno(data, id, gene_list)
      if (length(sub_data_cl$Pheno)==0){
        next
      }
      
      names(sub_data_cl)[names(sub_data_cl) == "Cell.X.Position"] <- "Cell X Position"
      names(sub_data_cl)[names(sub_data_cl) == "Cell.Y.Position"] <- "Cell Y Position"
      sub_data_cl$Pheno=gsub("+","+,",sub_data_cl$Pheno, fixed = T)
      sub_data_cl$Pheno = substr(sub_data_cl$Pheno, 1, nchar(sub_data_cl$Pheno)-1)
      
      log4r_info("Evaluating distances...")
      csd_with_distance = distance_eval(sub_data_cl)
      
      if (myData$Distance["save_csv"][[1]]){
        name_csv=paste0(id,"_Distance.txt")
        
        print(paste0("Saving ",name_csv," file"))
        write.table(csd_with_distance, file.path(out_f_gr,name_csv), append = FALSE, sep = "\t", dec = ".",
                    row.names = F, col.names = TRUE, quote = F)
      }
      cat("\n")
    }
  }
  log4r_info("End distance evaluation!")
  cat("\n")

  return(file.path(myData$Paths["output_folder"][[1]],"Log",paste0(log_name,".log")))
}

