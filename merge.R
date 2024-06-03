library(tidyverse)
library(rjson)
library(phenoptr)
library(phenoptrReports)
library(logr)
source("logger.R")

merge=function(){
  
  cat("\n################## MERGE ###################\n")
  
  log_name=paste0("Merge_",format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))
  logger(log_name)
  cat("\n")
  log4r_info("Start merging process: This step will merge together the cell_seg_data files for each subject")
  cat("\n")
  
  # Load json file
  log4r_info("Loading configuration file...")
  
  tryCatch(
    {
      myData=fromJSON(file="config.json")
    }, error=function(cond){
      log4r_error(paste0("Something went wrong while reading config.json file. Try check the fields of your json file!"))
      log4r_error("Exiting from clean_data.R script")
      stop()
    })
  
  #### Input parsing
  data_input_folder=file.path(myData$Paths["data_input_folder"][[1]],"raw_data")
  output_folder=myData$Paths["output_folder"][[1]]
  sample_sheet=file.path(myData$Paths["data_input_folder"][[1]],"sample_sheet.tsv")
  data_input=read.table(file = sample_sheet, sep = '\t', header = TRUE)
  
  #### Output organization folder
  dir.create(output_folder)
  dir.create(file.path(output_folder,"Merged"))
  log4r_info(paste0("Merged file will be stored in ",file.path(output_folder,"Merged")))
  cat("\n")
  
  #### Data merge
  log4r_info(paste0("Find ",nrow(data_input), " subjects!"))
  log4r_info("Start Merging process...")
  
  for (s in 1:nrow(data_input)){ #for loop on dataframe row (cells)
    
    cat("\n#############################\n")
    id=data_input$sbj_ID[s]
    group=data_input$Group[s]
  
    log4r_info(paste0("Analyzing Subject ", id, " - Group: ", group))
    
    group_out_path=file.path(output_folder,"Merged",group)
    dir.create(group_out_path)
    
    csv_file=list.files(data_input_folder, full.names = T)[grep(id, list.files(data_input_folder))]
    
    log4r_info("Check for already existing Merge_cell_seg_data.txt")
    
    if ("Merge_cell_seg_data.txt" %in% list.files(csv_file)){
      log4r_info("Removing old Merge_cell_seg_data.txt")
      file.remove(file.path(csv_file,"Merge_cell_seg_data.txt"))
    }
    
    cat("\n")
    
    files_merged=phenoptrReports::merge_cell_seg_files(csv_file)
    files_merged=list.files(csv_file,patter="Merge")
    
    if (is.null(files_merged)){
      log4r_warn(paste0("Something went wrong. Check the folder ",csv_file))
      next
      }
  
    seg_data=read.csv(file.path(csv_file,"Merge_cell_seg_data.txt"), sep="\t")
    pheno_list=colnames(seg_data[ , grepl( "Phenotype" , names( seg_data ) ) ])
    if (is.null(pheno_list)){
      pheno_list=seg_data[ , grepl( "Phenotype" , names( seg_data ) ) ]
    }
    seg_data$Pheno=seg_data[ ,  "Phenotype"]
    
    if (length(seg_data$Pheno)==0){
      log4r_warn("No Phenotypes found for this subject. Skip to the next one!")
      next
    }
    
    #write merged file
    cat("\n")
    write.table(seg_data, file.path(group_out_path,paste0("Merge_cell_seg_data_",id,".txt")), sep="\t", row.names = F, quote = F)
    log4r_info(paste0("File ",file.path(group_out_path,paste0("Merge_cell_seg_data_",id,".txt")), " saved correctly!"))
  
  } # end for loop
  
  cat("\n")
  log4r_info("End merge step!")
  cat("\n")
  
  return(file.path(myData$Paths["output_folder"][[1]],"Log",paste0(log_name,".log")))
}