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




