library(rjson)
library(phenoptr)
library(phenoptrReports)
library(plotly)
library(ggplot2)
library(dplyr)
library(htmlwidgets)

source("distance_match_func.R")
source("maps_plot_functions.R")
#load info from json
myData=fromJSON(file="config.json")

if (myData$Clean_data["other_rm"][[1]]){
  output_f=file.path(myData$Paths["output_folder"],"Merged_clean")
} else {
  output_f=file.path(myData$Paths["output_folder"],"Merged")
}

if (length(myData$Distance["pheno_list"][[1]])==0){
  gene_list=myData$Phenotypes["pheno_list"][[1]]
}else{
  gene_list=myData$Distance["pheno_list"][[1]]
}

groups=list.files(output_f)

#plot generation
for (group in groups){
  
  #group=groups[1] #da togliere
  
  print(paste0("Checking ",group))
  
  #creating output folders
  out_f=file.path(myData$Paths["output_folder"],"Distance")
  dir.create(out_f)
  out_f_gr=file.path(out_f,group)
  dir.create(out_f_gr)
  out_f_gr_csv=file.path(out_f_gr,"csv")
  dir.create(out_f_gr_csv)
  
  for (file in list.files(file.path(output_f,group), full.names = T)){
    
    print(paste0("Checking ",file))
    data<-read.delim(file)
    
    id=gsub(".txt","",strsplit(file,"_")[[1]][length(strsplit(file,"_")[[1]])])
    
    sub_data_cl=filt_data_Pheno(data, id, gene_list)
    if (length(sub_data_cl$Pheno)==0){
      next
    }
    
    names(sub_data_cl)[names(sub_data_cl) == "Cell.X.Position"] <- "Cell X Position"
    names(sub_data_cl)[names(sub_data_cl) == "Cell.Y.Position"] <- "Cell Y Position"
    sub_data_cl$Pheno=gsub("+","+,",sub_data_cl$Pheno, fixed = T)
    sub_data_cl$Pheno = substr(sub_data_cl$Pheno, 1, nchar(sub_data_cl$Pheno)-1)
    
    csd_with_distance = distance_eval(sub_data_cl)
    
    csd_with_distance %>% group_by(Phenotype) %>% 
      select(Phenotype, starts_with('Distance to')) %>% 
      summarize_all(~round(mean(.), 1))
    
    if (myData$Distance["save_csv"][[1]]){
      name_csv=paste0(id,"_Distance.txt")
      
      print(paste0("Saving ",name_csv))
      write.table(csd_with_distance, file.path(out_f_gr_csv,name_csv), append = FALSE, sep = "\t", dec = ".",
                  row.names = F, col.names = TRUE, quote = F)
    }
  }
}



