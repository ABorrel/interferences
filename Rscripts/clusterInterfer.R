#!/usr/bin/env Rscript
source("dendocircular.R")


################
#     MAIN     #
################

args <- commandArgs(TRUE)
pdesc = args[1]
pinter = args[2]#classes of interferences 
prout = args[3]


#pdesc = "/home/borrela2/interference/ToxCast_analysis/IC50Cluster/detection_technology_type/Interfers_-2_4/descClean.csv"
#pinter = "/home/borrela2/interference/ToxCast_analysis/IC50Cluster/detection_technology_type/Interfers_-2_4/Zscore-2_4"
#prout = "/home/borrela2/interference/ToxCast_analysis/IC50Cluster/detection_technology_type/Interfers_-2_4/"


##############################
# Process descriptors matrix #
##############################
ddesc = read.csv(pdesc, header = TRUE)
rownames(ddesc) = ddesc[,1]
ddesc = ddesc[,-1]
dinter = read.csv(pinter, sep = "\t", header = TRUE)

ltype = unique(dinter[,"Interfer"])
for(t in ltype){
  print(t)
  dinterSp = dinter[which(dinter[,"Interfer"] == t),]
  rownames(dinterSp) = dinterSp[,1]
  lchem = intersect(dinterSp[,1], rownames(ddesc))
  ddesctemp = ddesc[lchem,]
  dinterSp = dinterSp[lchem,]
  if(dim(ddesctemp)[1] >2){
    dendogramInterfer(ddesctemp, dinterSp, t, prout)

    }
}

