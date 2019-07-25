#!/usr/bin/env Rscript
source("dendocircular.R")


################
#     MAIN     #
################

args <- commandArgs(TRUE)
pdesc = args[1]
pAC50Luc = args[2]
pAC50Hek = args[3]
pAC50HepG = args[4]
prresult = args[5]

#pdesc = "/home/borrela2/interference/spDataAnalysis/clusters/interfer/13_15/desc.csv"
#pAC50Luc = "/home/borrela2/interference/spDataAnalysis/tox21-luc-biochem-p1/AC50_combine"
#pAC50Hek = "/home/borrela2/interference/spDataAnalysis/tox21-spec-hek293-p1/AC50_sample"
#pAC50HepG = "/home/borrela2/interference/spDataAnalysis/tox21-spec-hepg2-p1/AC50_sample"
#prresult = "/home/borrela2/interference/spDataAnalysis/clusters/interfer/13_15/"
  
ddesc = read.csv(pdesc, header = TRUE, sep = "\t")
rownames(ddesc) = ddesc[,1]

dAC50luc = read.csv(pAC50Luc, header = TRUE, sep = "\t")
rownames(dAC50luc) = dAC50luc[,1]

dAC50Hek = read.csv(pAC50Hek, header = TRUE, sep = "\t")
rownames(dAC50Hek) = dAC50Hek[,1]

dAC50HepG = read.csv(pAC50HepG, header = TRUE, sep = "\t")
rownames(dAC50HepG) = dAC50HepG[,1]

ltypecolor = colnames(dAC50Hek)[-1]

for (typecol in ltypecolor){
  dACplot = dAC50luc[rownames(ddesc),2]
  dACplot = cbind(dACplot, dAC50Hek[rownames(ddesc), typecol])
  dACplot = cbind(dACplot, dAC50HepG[rownames(ddesc), typecol])
  rownames(dACplot) = rownames(ddesc)
  
  colnames(dACplot) = c("Luc", "hek293", "hepg2")
  dACplot = as.data.frame(dACplot)
  dACplot = -log10(dACplot)
  
  dcluster = rep(1,dim(dACplot)[1])
  dcluster = cbind(rownames(dACplot), dcluster)
  rownames(dcluster) = rownames(dACplot)
  colnames(dcluster) = c("ID", "cluster")
  dcluster = as.data.frame(dcluster)
  ddesc = ddesc[,-1]
  
  dendogramCluster(ddesc, dACplot, dcluster, paste(prresult, typecol, sep = ""))
  
}

#print(dim(dAC50HepG))
#print(dim(dAC50Hek))
#print(dim(dAC50luc))

#print(dim(ddesc))



  