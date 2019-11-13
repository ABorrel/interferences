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

generateSOM = function(ddesc, xdim, ydim, prout, modelin){
  ddesc = as.matrix(scale(ddesc))
  if(is.double(modelin)){
    som_grid <- somgrid(xdim=xdim, ydim=ydim, topo="hexagonal")
    som_model <- som(ddesc, 
                     grid=som_grid, 
                     rlen=100, 
                     alpha=c(0.05,0.01), 
                     keep.data = TRUE)  
  }else{
    som_model = modelin
  }
  
  svg(paste(prout, "SOM_count.svg", sep = ""))
  plot(som_model, type = "count", palette.name=coolBlueHotRed, main = "", dpi=300, height = 20, width = 20, bg = "transparent")
  dev.off()
  
  
  #write cluster
  dclust = som_model$unit.classif
  names(dclust) = rownames(ddesc)
  write.csv(dclust, paste(prout, "SOMClust", sep = ""))
  
  return (som_model)
}



################
#     MAIN     #
################

args <- commandArgs(TRUE)
pdesc = args[1]
prout = args[2]
size = args[3]
pmodel = args[4]

#pdesc = "/home/borrela2/interference/spDataAnalysis/tox21-spec-hepg2-p1/SOM/descClean.csv"
#prout = "/home/borrela2/interference/spDataAnalysis/tox21-spec-hepg2-p1/SOM/"
#size = 15
if(pmodel != "0"){
  load(pmodel)  
}else{
  model = 0.0
}


##############################
# Process descriptors matrix #
##############################
ddesc = read.csv(pdesc, header = TRUE)
rownames(ddesc) = ddesc[,1]
ddesc = ddesc[,-1]

model = generateSOM(ddesc, size, size, prout, model)

save(model, file = paste(prout, "SOMmodel.Rdata", sep = ""))
