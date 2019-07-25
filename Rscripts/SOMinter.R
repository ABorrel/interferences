#!/usr/bin/env Rscript
require(kohonen)
library(ggplot2)

#################
# Color reverse #
#################

colors <- function(n, alpha = 1) {
  rev(heat.colors(n, alpha))
}

coolBlueHotRed <- function(n, alpha = 1) {rainbow(n, end=4/6, alpha=alpha)[n:1]}


#######
# SOM #
#######

generateSOMinter = function(ddesc, dinter, som_model, prout){
  
  SOM_temp = som_model
  rownames(dinter) = dinter[,1]
  
  dclust = som_model$unit.classif
  names(dclust) = rownames(ddesc)
  
  
  lcas = NULL
  for(cas in rownames(dinter)){
    lcas = append(lcas, dclust[which(names(dclust) == cas)])
  }
  
  write.csv(lcas, paste(prout, "interClust.csv", sep = ""))
  SOM_temp$unit.classif = as.vector(lcas)
  
  png(paste(prout, ".png", sep = ""))
  plot(SOM_temp, type = "counts", palette.name = colors, heatkey = TRUE)
  dev.off()
  
}


################
#     MAIN     #
################

args <- commandArgs(TRUE)
pdesc = args[1]
pinter = args[2]
pmodel = args[3]
prout = args[4]

#pdesc = "/home/borrela2/interference/ToxCast_analysis/IC50Cluster/detection_technology_type/Interfers_-2_4/SOM/descClean.csv"
#pinter = "/home/borrela2/interference/ToxCast_analysis/IC50Cluster/detection_technology_type/Interfers_-2_4/Zscore-2_4"
#prout = "/home/borrela2/interference/ToxCast_analysis/IC50Cluster/detection_technology_type/Interfers_-2_4/SOM/"
  
load(pmodel)

##############################
# Process descriptors matrix #
##############################
ddesc = read.csv(pdesc, header = TRUE)
rownames(ddesc) = ddesc[,1]
ddesc = ddesc[,-1]

dinter = read.csv(pinter, sep = "\t", header = TRUE)
lintertype = unique(dinter[,c("Interfer")])

for(intertype in lintertype){
  generateSOMinter(ddesc, dinter[which(dinter[,c("Interfer")] == intertype),], model, paste(prout, intertype, sep = ""))
}

