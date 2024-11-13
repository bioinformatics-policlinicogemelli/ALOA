#Copyright 2024 bioinformatics-policlinicogemelli

#Licensed under the Apache License, Version 2.0 (the "License");
#you may not use this file except in compliance with the License.
#You may obtain a copy of the License at

#    http://www.apache.org/licenses/LICENSE-2.0

#Unless required by applicable law or agreed to in writing, software
#distributed under the License is distributed on an "AS IS" BASIS,
#WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#See the License for the specific language governing permissions and
#limitations under the License.

library(rjson)
library(phenoptr)
library(phenoptrReports)
library(plotly)
library(ggplot2)
library(dplyr)
library(htmlwidgets)
library(foreach)
vignette("nested")

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
  log4r_info(paste0("Distances output file will be stored in ", out_f))
  
  output_f=file.path(myData$Paths["output_folder"],"Merged_clean")

  if (length(myData$Distance["pheno_list"][[1]])==0){
    gene_list=myData$Phenotypes["pheno_list"][[1]]
  }else{
    gene_list=myData$Distance["pheno_list"][[1]]
  }
  
  if (length(gene_list)){
    log4r_info
  }
  
  groups=list.files(output_f)
  
  #plot generation
  for (group in groups){
    gc()
    cat("\n###################\n")
    log4r_info(paste0("Checking Group: ",group))
    
    #creating output folders
    out_f_gr=file.path(out_f,group)
    dir.create(out_f_gr)
    
    list_of_files=list.files(file.path(output_f,group), full.names = T)
    foreach(i=1:NROW(list_of_files), .combine=c) %do%{
      file=list_of_files[i]
    #for(file in list_of_files){  
      log4r_info(paste0("Checking ",file))
      #data<-read.delim(file)
      data <- data.table(fread(file, showProgress = FALSE))
      data=data[,c("Cell.ID","Cell.X.Position", "Cell.Y.Position", "Pheno")]

      id=gsub(".txt","",gsub("Merge_cell_seg_data_clean_","",basename(file)))
      
      sub_data_cl=filt_data_Pheno(data, id, gene_list)
      if (length(sub_data_cl$Pheno)==0){
        next
      }
      #sub_data_cl=sub_data_cl[,c("Cell.ID","Cell.X.Position", "Cell.Y.Position", "Pheno")]
      #sub_data_cl=sub_data_cl[c("Cell.ID","Cell.X.Position", "Cell.Y.Position", "Pheno")]
      rm(data)
      
      names(sub_data_cl)[names(sub_data_cl) == "Cell.X.Position"] <- "Cell X Position"
      names(sub_data_cl)[names(sub_data_cl) == "Cell.Y.Position"] <- "Cell Y Position"
      sub_data_cl$Pheno=gsub("+","+,",sub_data_cl$Pheno, fixed = T)
      sub_data_cl$Pheno = substr(sub_data_cl$Pheno, 1, nchar(sub_data_cl$Pheno)-1)
      gc()
      log4r_info("Evaluating distances...")
      csd_with_distance = distance_eval(sub_data_cl)
      
      if (myData$Distance["save_csv"][[1]]){
        name_csv=paste0(id,"_Distance.txt")
        
        print(paste0("Saving ",name_csv," file"))
        write.table(csd_with_distance, file.path(out_f_gr,name_csv), append = FALSE, sep = "\t", dec = ".",
                    row.names = F, col.names = TRUE, quote = F)
      }
      rm(csd_with_distance)
      gc()
      cat("\n")
    }
  }
  log4r_info("End distance evaluation!")
  cat("\n")

  return(file.path(myData$Paths["output_folder"][[1]],"Log",paste0(log_name,".log")))
}

