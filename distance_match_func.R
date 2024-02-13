library(phenoptr)
library(tiff)
library(ggplot2)
library(rjson)
library(tidyverse)

distance_eval=function(csd){
  
  names(csd)[names(csd) == "Pheno"] <- "Phenotype"
  names(csd)[names(csd) == "Cell.ID"] <- "Cell ID"
  
  distances <- find_nearest_distance(csd)
  csd_with_distance <- bind_cols(csd, distances)
  
  return (csd_with_distance)
  
}




