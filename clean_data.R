library(readxl)
library(tidyverse)
library(rjson)
library(phenoptr)
library(phenoptrReports)

source("logger.R")
logger(paste0("clean_data_log_",format(Sys.time(), "%Y-%m-%d")))

log4r_info("Start cleaning process: This step will remove cell negative to all phenotype(s)")

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
output_folder=myData$Paths["output_folder"][[1]]
sample_sheet=myData$Paths["sample_sheet"][[1]]

data_input_folder=file.path(output_folder,"Merged")

#### Output organization folder
dir.create(output_folder)
dir.create(file.path(output_folder,"Merged_clean"))
log4r_info(paste0("Clean merged file will be stored in ",file.path(output_folder,"Merged_clean")))

group= list.files(data_input_folder)

if (length(group)==0){
  log4r_warn("no group found. Check the output folder!")
  log4r_info("Exit clean_data script!")
  stop()
}

#### Data merge
for (g in group){ #for loop on dataframe row (cells)
  log4r_info("#############################")
  log4r_info(paste0("Group: ",g))
  log4r_info("#############################")
  group_out_path=file.path(output_folder,"Merged_clean",g)
  dir.create(group_out_path)
  
  csv_file=list.files(file.path(data_input_folder,g), full.names = T)
  
  if (length(csv_file)==0){
    log4r_warn(paste0("No csv file found for the group ", g,". Check the output folder!"))
    log4r_info("Skipping to the next group!")
    next
  }
  #i=1
  for (csv in csv_file){
    
    id=strsplit(sapply(strsplit(csv,"_"),tail,1),"[.]")[[1]][1]
    #id=paste0("roi_",i)
    
    log4r_info(paste0(" -> Subject: ",id))
    seg_data=read.csv(csv, sep="\t")
    
    #if(!"Pheno" %in% colnames(seg_data)){
    pheno_list=colnames(seg_data[ , grepl( "Pheno" , names( seg_data ) ) ])
    pheno_list=pheno_list[pheno_list != "Pheno"]
    if (length(pheno_list)>1){
      seg_data$Pheno=apply( seg_data[ ,  pheno_list] , 1 , paste , collapse = "," )
    }else{
      seg_data$Pheno=seg_data[ ,  pheno_list]
    }
    #}
    
    #pheno_list=colnames(seg_data[ , grepl( "Pheno" , names( seg_data ) ) ])
    
    #pheno_list=colnames(seg_data[ , grepl( "Phenotype" , names( seg_data ) ) ])
    seg_clean=subset(seg_data,seg_data$Pheno!=str_c(rep("other",each=length(pheno_list)),collapse=","))
    seg_clean=subset(seg_data,seg_data$Pheno!=str_c(rep("OTHER",each=length(pheno_list)),collapse=","))
    seg_clean=subset(seg_clean,seg_clean$Pheno!=str_c(rep(",",each=length(pheno_list)-1),collapse=""))
    
    log4r_info(paste0("     All empty/all other Row removed: ", nrow(seg_data)-nrow(seg_clean), " of ", nrow(seg_data)))
    
    log4r_info(paste0("Writing ", "Merge_cell_seg_data_clean_",id,".txt file..."))
    write.table(seg_clean, file.path(group_out_path,paste0("Merge_cell_seg_data_clean_",id,".txt")), sep="\t", row.names = F, quote = F)
    #i=i+1
    
  }
}

print("Clean Step completed!")
