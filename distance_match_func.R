library(phenoptr)
library(tiff)
library(ggplot2)
library(rjson)
library(tidyverse)

distance_eval=function(csd){
  
  phenocol=colnames(csd[ , grepl( "Pheno" , names( csd ) ) ])
  if (!("Phenotype" %in% phenocol)){
    names(csd)[names(csd) == "Pheno"] <- "Phenotype"
  }
  
  names(csd)[names(csd) == "Cell.ID"] <- "Cell ID"
  names(csd)[names(csd) == "Cell.X.Position"] <- "Cell X Position"
  names(csd)[names(csd) == "Cell.Y.Position"] <- "Cell Y Position"
  
  distances <- find_nearest_distance(csd)
  csd_with_distance <- bind_cols(csd, distances)
  
  return (csd_with_distance)
  
}




