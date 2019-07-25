#!/usr/bin/env Rscript
library(ggplot2)
source("~/development/Rglobal/source/PCAdrawer.R")

PCAplot = function (din, vcolor, prout){
  coord = generatePCAcoords(din)
  data_plot = coord[[1]]
  var_cap = coord[[2]]
  cp = coord[[3]]
  
  col.desc = "black"
  
  
  png (paste (prout, "PCA_text.png", sep = ""), 1700, 1500)
  factor = factorACP (data_plot, cp)
  
  color_arrow = col.desc[rownames(cp)]
  par(mar=c(8,8,8,8))
  plot(data_plot[,1],data_plot[,2], pch=20, main = paste (var_cap[1],var_cap[2], sep = "_" ), xlab = paste("PC1: ", signif (var_cap[1], 4), "%", sep = ""), ylab = paste("PC2: ", signif (var_cap[2], 4), "%", sep = ""), cex.lab = 4, cex.main = 4, cex.axis = 1.75, cex = 4, type = "n")
  text (data_plot[,1],data_plot[,2], label = rownames (din), cex = 1.2)
  abline(h=0,v=0)
  warnings ()
  dev.off()
  
  print (length(vcolor))
  print (dim(data_plot))
  
  png (paste (prout, "PCA_color.png", sep = ""), 1700, 1500)
  factor = factorACP (data_plot, cp)
  color_arrow =col.desc
  par(mar=c(8,8,8,8))
  plot(data_plot[,1],data_plot[,2], pch=20, col = vcolor, xlab = paste("PC1: ", round (var_cap[1], 2), "%", sep = ""), ylab = paste("PC2: ", round(var_cap[2], 2), "%", sep = ""), cex.lab = 4, cex.main = 4, cex.axis = 1.75, cex = 4)
  points(data_plot[which(vcolor != "gray90"),1], data_plot[which(vcolor != "gray90"),2], pch=20, col = vcolor[which(vcolor != "gray90")], cex = 4)
  abline(h=0,v=0)
  warnings ()
  dev.off()
  
  
  png (paste (prout, "PCA_descriptor.png", sep = ""), 1700, 1500)
  par(mar=c(8,8,8,8))
  plot(data_plot[,1],data_plot[,2], xlab = paste("CP1: ", signif (var_cap[1], 4), "%", sep = ""), ylab = paste("CP2: ", signif (var_cap[2], 4), "%", sep = ""), pch=20, cex.lab = 4, cex.main = 4, cex.axis = 1.75, cex = 4, type = "n")
  #points (data_plot[,1][length (color_point):dim(data_plot)[1]],data_plot[,2][length (color_point):dim(data_plot)[1]], pch=17, cex = 4, col = color_point2)
  abline(h=0,v=0)
  arrows (0,0,cp[,1]*factor,cp[,2]*factor, col = color_arrow, lwd = 4 )
  text (cp[,1]*factor, cp[,2]*factor,rownames(cp), col = color_arrow, cex = 2.5)
  dev.off()
  
  
  svg (file = paste (prout, "PCA_descriptor.svg", sep = ""), 25, 25)
  par(mar=c(8,8,8,8))
  plot(data_plot[,1],data_plot[,2], main = "", xlab = paste("CP1: ", round (var_cap[1], 2), "%", sep = ""), ylab = paste("CP2: ", round (var_cap[2], 2), "%", sep = ""), cex.lab = 6, cex.main = 4, cex.axis = 1.75, cex = 6, type = "n")
  #points (data_plot[,1][length (color_point):dim(data_plot)[1]],data_plot[,2][length (color_point):dim(data_plot)[1]], pch=17, cex = 4, col = color_point2)
  abline(h=0,v=0)
  arrows (0,0,cp[,1]*factor,cp[,2]*factor, col = color_arrow, lwd = 3 )
  text (cp[,1]*factor, cp[,2]*factor,rownames(cp), col = color_arrow, cex = 3.5)
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
#pAC50 = "/home/borrela2/interference/spDataAnalysis/tox21-spec-hepg2-p1/AC50_sample_curve_eff_burst"
#prout = "/home/borrela2/interference/spDataAnalysis/tox21-spec-hepg2-p1/PCA/"
  
#pdesc = "/home/borrela2/interference/spDataAnalysis/tox21-luc-biochem-p1/PCA/descClean.csv"
#pAC50 = "/home/borrela2/interference/spDataAnalysis/tox21-luc-biochem-p1/AC50_combine"
#prout = "/home/borrela2/interference/spDataAnalysis/tox21-luc-biochem-p1/PCA/"

#pdesc = "/home/borrela2/interference/spDataAnalysis/tox21-spec-hepg2-p1/PCA/descClean.csv"
#pAC50 = "/home/borrela2/interference/spDataAnalysis/tox21-spec-hepg2-p1/AC50_sample_curve"
#prout = "/home/borrela2/interference/spDataAnalysis/tox21-spec-hepg2-p1/PCA/"

ddesc = read.csv(pdesc, sep = ",", header = TRUE)
rownames(ddesc) = ddesc[,1]
ddesc = ddesc[,-1]
dAC50 = read.csv(pAC50, sep = "\t", header = TRUE)
print(head(dAC50))

rownames(dAC50) = dAC50[,1]
#dAC50 = dAC50[,-1]


lID = intersect(rownames(dAC50), rownames(ddesc))
print(lID)

ddesc = ddesc[lID,]
dAC50 = dAC50[lID,]

print (dim(dAC50))
print (dim(ddesc))

if (dim(dAC50)[2] == 2){
  dAC50 = as.matrix(dAC50)
  vcolor = rep("gray90", dim(ddesc)[1])
  print(which(dAC50[,2] >= 0))
  vcolor[which(dAC50[,2] >= 0)] = "gray25"
  print (vcolor)
  
}else{
  vcolor = rep("gray90", dim(ddesc)[1])
  vcolor[which(dAC50[,"cell_blue_n"]>= 0)] = "cyan"
  vcolor[which(dAC50[,"med_blue_n"]>= 0)] = "blue"
  vcolor[which(dAC50[,"med_green_n"]>= 0)] = "darkgreen"
  vcolor[which(dAC50[,"cell_green_n"]>= 0)] = "darkolivegreen3"
  vcolor[which(dAC50[,"med_red_n"]>= 0)] = "red"
  vcolor[which(dAC50[,"cell_red_n"]>= 0)] = "firebrick"
}

PCAplot(ddesc, vcolor, prout)
