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

library(readxl)
library(tidyverse)
library(rjson)
library(phenoptr)
library(phenoptrReports)
library(fs)

source("logger.R")

clean=function(){
  
  cat("\n################## CLEAN DATA ###################\n")
  
  log_name=paste0("clean_data_",format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))
  logger(log_name)
  cat("\n")
  log4r_info("Start cleaning process: This step will remove cell negative to all phenotype(s)")
  cat("\n")
  
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
    return(file.path(myData$Paths["output_folder"][[1]],"Log",paste0(log_name,".log")))
  }

  #### Data merge
  for (g in group){ #for loop on dataframe row (cells)
    
    cat("\n#############################\n")
    log4r_info(paste0("Group: ",g))
    group_out_path=file.path(output_folder,"Merged_clean",g)
    dir.create(group_out_path)
    
    csv_file=list.files(file.path(data_input_folder,g), full.names = T)
    
    if (length(csv_file)==0){
      log4r_warn(paste0("No csv file found for the group ", g,". Check the output folder!"))
      log4r_info("Skipping to the next group!")
      next
    }

    for (csv in csv_file){
      cat("\n")
      id=gsub(".txt","",gsub("Merge_cell_seg_data_","",path_file(csv)))

      log4r_info(paste0("Subject: ",id))
      seg_data=read.csv(csv, sep="\t")

      pheno_list=colnames(seg_data[ , grepl( "Pheno" , names( seg_data ) ) ])
      pheno_list=pheno_list[pheno_list != "Pheno"]
      if (length(pheno_list)>1){
        
        patterns_to_remove <- c("other", "OTHER", "others", "OTHERS", "Other")
        seg_clean <- subset(seg_data, !seg_data$Pheno %in% sapply(patterns_to_remove, function(x) str_c(rep(x, each = length(pheno_list)), collapse = ",")))

        for (pattern in patterns_to_remove) {
          seg_clean$Pheno <- gsub(paste0("[, ]*", pattern, "[, ]*"), "", seg_clean$Pheno, ignore.case = TRUE)
        }

        seg_clean$Pheno <- gsub(",", "", seg_clean$Pheno, fixed = TRUE)
        seg_clean$Pheno <- trimws(seg_clean$Pheno)
        to_exclude <- c("", "other")
        seg_clean <- seg_clean[!(is.na(seg_clean$Pheno) | seg_clean$Pheno %in% to_exclude), ]

        # seg_clean=subset(seg_data,seg_data$Pheno!=str_c(rep("other",each=length(pheno_list)),collapse=","))
        # seg_clean=subset(seg_clean,seg_clean$Pheno!=str_c(rep("OTHER",each=length(pheno_list)),collapse=","))
        # seg_clean=subset(seg_clean,seg_clean$Pheno!=str_c(rep("others",each=length(pheno_list)),collapse=","))
        # seg_clean=subset(seg_clean,seg_clean$Pheno!=str_c(rep("OTHERS",each=length(pheno_list)),collapse=","))
        # seg_clean=subset(seg_clean,seg_clean$Pheno!=str_c(rep("Other",each=length(pheno_list)),collapse=","))
        # seg_clean=subset(seg_clean,seg_clean$Pheno!=str_c(rep(",",each=length(pheno_list)-1),collapse=""))
        
        # seg_clean$Pheno=gsub("OTHERS","",seg_clean$Pheno)
        # seg_clean$Pheno=gsub("others,","",seg_clean$Pheno)
        # seg_clean$Pheno=gsub("OTHER","",seg_clean$Pheno)
        # seg_clean$Pheno=gsub("other,","",seg_clean$Pheno)
        # seg_clean$Pheno=gsub("Other","",seg_clean$Pheno)
        # seg_clean$Pheno=gsub(",OTHERS","",seg_clean$Pheno)
        # seg_clean$Pheno=gsub(",others","",seg_clean$Pheno)
        # seg_clean$Pheno=gsub(",OTHER","",seg_clean$Pheno)
        # seg_clean$Pheno=gsub(",other","",seg_clean$Pheno)
        # seg_clean$Pheno=gsub("Other,","",seg_clean$Pheno)
        # seg_clean$Pheno=gsub(",Other","",seg_clean$Pheno)
        
        # seg_clean$Pheno=gsub(",","",seg_clean$Pheno, fixed = T)
        # seg_clean=seg_clean[!(is.na(seg_clean$Pheno) | seg_clean$Pheno=="" | seg_clean$Pheno=="other"), ]
        
      }else if ((length(pheno_list)==1) & "Phenotype" %in% pheno_list){
        seg_data$Pheno=seg_data$Phenotype
        
        seg_clean=subset(seg_data, seg_data$Pheno!="other" & seg_data$Pheno!="OTHER" & seg_data$Pheno!="others" & seg_data$Pheno!="OTHERS")
      }
      
      

      log4r_info(paste0("All empty/only other Row removed: ", nrow(seg_data)-nrow(seg_clean), " of ", nrow(seg_data)))
      
      log4r_info(paste0("Writing ", "Merge_cell_seg_data_clean_",id,".txt file..."))
      write.table(seg_clean, file.path(group_out_path,paste0("Merge_cell_seg_data_clean_",id,".txt")), sep="\t", row.names = F, quote = F)
      
    }
  }
  log4r_info("End cleaning merge step!")
  cat("\n")
  
  return(file.path(myData$Paths["output_folder"][[1]],"Log",paste0(log_name,".log")))
}
