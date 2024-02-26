library(readxl)
library(tidyverse)
library(phenoptr)
library(phenoptrReports)
library(ggplot2)
library(colorRamp2)
library(dplyr)
library(plotly)
library(htmlwidgets)
library(data.table)

filt_data_Pheno=function(data, id, gene_list){
  
  #clean "other" tag from data$Pheno column
  data$Pheno=gsub("OTHER","",data$Pheno)
  data$Pheno=gsub("other,","",data$Pheno)
  data$Pheno=gsub(",other","",data$Pheno)
  data$Pheno=gsub(",","",data$Pheno, fixed = T)
  data=data[!(is.na(data$Pheno) | data$Pheno=="" | data$Pheno=="other"), ]
  
  pheno_list_mod <- list()
  for( p in 1:length(gene_list)){
    
    if (length(Reduce(intersect,(list(gene_list[p], unique(data$Pheno)))))>0){
      pheno_list_mod[p]=gene_list[p]
    }else if (lengths(str_split(gene_list[p],","))>1){
      l=str_split(gene_list[p],",")
      cmd=paste0("crossing(",paste(rep("l[[1]]",lengths(l)), collapse=","),")")
      comb=eval(parse(text=cmd))
      
      list_c=apply( comb , 1 , paste , collapse = "" )

      for (i in 1:length(list_c)){
        cc=list_c[i]
        if (cc %in% unique(data$Pheno)){
          pheno_list_mod[p]=cc
          }
      }
    }
  }

  idx=list()
  i=1
  for (cc in pheno_list_mod){
    print(cc)

    idx[[i]]=which(data$Pheno == cc)
    i=i+1
  }
  
  idx=unlist(idx)
  
  return(data[idx,])
}

multi_maps_plot=function(data, id, out_folder, interact=F){
  
  gene_list=unique(data$Pheno)
  
  pdf(file = file.path(out_folder,paste0(id,"_All_Pheno_",paste(gene_list, collapse=''),".pdf")))
  #jpeg(file = file.path(out_folder,paste0(id,"_All_Pheno_",paste(gene_list, collapse=''),".jpeg")))
  print(ggplot(data)+
          geom_point(aes(`Cell.X.Position`,`Cell.Y.Position`, color=Pheno),size=2, alpha = 0.6)+
          xlab("X Position") +
          ylab("Y Position") +
          guides(color = guide_legend("Phenotypes",override.aes = list(size = 5)))+
          theme(panel.background = element_rect(fill = 'white'),legend.text = element_text(size=10)))
  
  dev.off()
}
