#!/usr/bin/env Rscript
require(kohonen)


generateModelSOM = function(ddesc, ){
  
  
  
}





################
#     MAIN     #
################

args <- commandArgs(TRUE)
pdesc = args[1]
prout = args[2]

pdesc = "/home/borrela2/interference/spDataAnalysis/SOM/descClean.csv"
prout = "/home/borrela2/interference/spDataAnalysis/SOM/"

##############################
# Process descriptors matrix #
##############################
ddesc = read.csv(pdesc, header = TRUE)
rownames(ddesc) = ddesc[,1]
ddesc = ddesc[,-1]


desc = as.matrix(scale(ddesc))

som_grid <- somgrid(xdim = 40, ydim=40, topo="hexagonal")

som_model <- som(desc, 
                 grid=som_grid, 
                 rlen=100, 
                 alpha=c(0.05,0.01), 
                 keep.data = TRUE)#,
   #              n.hood="circular")



plot(som_model, type="counts")

plot(som_model, type="codes")
