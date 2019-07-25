#!/usr/bin/env Rscript
library(ggplot2)




MDSplot = function(din, vcol, type_dist, prout){
  
  
  din = na.omit (din)
  if (type_dist == "corr"){
    MC = cor(din)
    dist1 = abs(1-MC)  
  }else{
    dist1 = dist(din, method = typeMDS)
  }
  name_descriptor = colnames(din)
  
  color_desc = "black" #colorDesc(name_descriptor)
  #print (color_desc)
  
  png (paste (prout, "MDS_desc_", type_dist, ".png", sep = ""), 3500, 3500, res = 180, bg = "white")
  par( mar=c(6,6,6,6))
  fit <- cmdscale(as.dist(dist1),eig=TRUE, k=2)
  #c = as.vector (col.desc[,colnames (data)])
  plot (fit$points[,1], fit$points[,2], main="MDS Descriptor", xlab = "DIM 1", ylab = "DIM 2", cex.lab = 2.75, cex.axis = 2, cex.main = 2.75, type = "n")
  #text (fit$points[,1], fit$points[,2]+0.05, labels = num_desc,  cex = 1.8, col = color_desc, font = 12)
  #text (fit$points[descriptor_selected,1], fit$points[descriptor_selected,2]+0.02, labels = "*",  cex = 4, col = color_desc)
  text (fit$points[,1], fit$points[,2], labels = name_descriptor,  cex = 2.6, col = color_desc)
  dev.off()
  
  
  din = t(din)
  if (type_dist == "corr"){
    MC = cor(din)
    dist1 = abs(1-MC)  
  }else{
    dist1 = dist(din, method = typeMDS)
  }
  
  png (paste (prout, "MDS_chem_", type_dist, ".png", sep = ""), 3500, 3500, res = 180, bg = "white")
  par( mar=c(6,6,6,6))
  fit <- cmdscale(as.dist(dist1),eig=TRUE, k=2)
  #c = as.vector (col.desc[,colnames (data)])
  plot (fit$points[,1], fit$points[,2], main="MDS chemicals", xlab = "DIM 1", ylab = "DIM 2", cex.lab = 2.75, cex.axis = 2, cex.main = 2.75, col = vcol)
  dev.off()
  
}

################
#     MAIN     #
################

args <- commandArgs(TRUE)
pdesc = args[1]
pAC50 = args[2]
prout = args[3]

#pdesc = "/home/borrela2/interference/spDataAnalysis/tox21-spec-hepg2-p1/PCA/descClean.csv"
#pAC50 = "/home/borrela2/interference/spDataAnalysis/tox21-spec-hepg2-p1/AC50_sample_curve"
#prout = "/home/borrela2/interference/spDataAnalysis/tox21-spec-hepg2-p1/PCA/"

ddesc = read.csv(pdesc, sep = ",", header = TRUE)
rownames(ddesc) = ddesc[,1]
ddesc = ddesc[,-1]
dAC50 = read.csv(pAC50, sep = "\t", header = TRUE)

rownames(dAC50) = dAC50[,1]
dAC50 = dAC50[,-1]


dAC50 = dAC50[rownames(ddesc),]

vcolor = rep("gray90", dim(ddesc)[1])
vcolor[which(dAC50[,"cell_blue"]>= 0)] = "cyan"
vcolor[which(dAC50[,"med_blue"]>= 0)] = "blue"
vcolor[which(dAC50[,"med_green"]>= 0)] = "darkgreen"
vcolor[which(dAC50[,"cell_green"]>= 0)] = "darkolivegreen3"
vcolor[which(dAC50[,"med_red"]>= 0)] = "red"
vcolor[which(dAC50[,"cell_red"]>= 0)] = "firebrick"



MDSplot(ddesc, vcolor, "corr", prout)
