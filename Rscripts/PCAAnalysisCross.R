#!/usr/bin/env Rscript
library(ggplot2)
source("~/development/Rglobal/source/PCAdrawer.R")


PCAplot = function (din, vcolor, vpoint, prout){
  
  coord = generatePCAcoords(din)
  data_plot = coord[[1]]
  var_cap = coord[[2]]
  cp = coord[[3]]
  
  col.desc = "black"
  
  
  png (paste (prout, "PCA_text.png", sep = ""), 1700, 1500)
  factor = factorACP (data_plot, cp)
  
  color_arrow = col.desc[rownames(cp)]
  par(mar=c(8,8,8,8))
  plot(data_plot[,1],data_plot[,2], pch=20, main = paste (var_cap[1],var_cap[2], sep = "_" ), xlab = paste("CP1: ", signif (var_cap[1], 4), "%", sep = ""), ylab = paste("CP2: ", signif (var_cap[2], 4), "%", sep = ""), cex.lab = 4, cex.main = 4, cex.axis = 1.75, cex = 4, type = "n")
  text (data_plot[,1],data_plot[,2], label = rownames (din), cex = 1.2)
  abline(h=0,v=0)
  warnings ()
  dev.off()
  
  
  png (paste (prout, "PCA_color.png", sep = ""), 1700, 1500)
  factor = factorACP (data_plot, cp)
  color_arrow =col.desc
  par(mar=c(8,8,8,8))
  plot(data_plot[,1],data_plot[,2], pch=vpoint, col = vcolor, xlab = paste("CP1: ", round (var_cap[1], 2), "%", sep = ""), ylab = paste("CP2: ", round(var_cap[2], 2), "%", sep = ""), cex.lab = 4, cex.main = 4, cex.axis = 1.75, cex = 2)
  points(data_plot[which(vcolor != "gray90"),1], data_plot[which(vcolor != "gray90"),2], pch=vpoint[which(vcolor != "gray90")], col = vcolor[which(vcolor != "gray90")], cex = 2.5)
  abline(h=0,v=0)
  legend(-40, -2, legend=c("Hepg2", "Hek293",  "Hepg2+Hek293", "cell based", "cell free",  "cell based", "cell free",  "cell based", "cell free"),
         col=c("black", "black", "black", "blue", "cyan", "darkgreen", "darkolivegreen3", "red", "firebrick"), cex=2.5,
         box.lty=2, box.lwd=2, box.col="grey", pch = c(19, 18, 8, 15, 15, 15, 15, 15, 15))
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
pAC501 = args[2]
pAC502 = args[3]
prout = args[4]

#pdesc = "/home/borrela2/interference/spDataAnalysis/CrossPCA/descClean.csv"
#pAC501 = "/home/borrela2/interference/spDataAnalysis/tox21-spec-hepg2-p1/AC50_sample_curve"
#pAC502 = "/home/borrela2/interference/spDataAnalysis/tox21-spec-hek293-p1/AC50_sample_curve"
#prout = "/home/borrela2/interference/spDataAnalysis/CrossPCA/"

ddesc = read.csv(pdesc, sep = ",", header = TRUE)
rownames(ddesc) = ddesc[,1]
ddesc = ddesc[,-1]

dAC501 = read.csv(pAC501, sep = "\t", header = TRUE)
rownames(dAC501) = dAC501[,1]
dAC501 = dAC501[,-1]
dAC501 = dAC501[rownames(ddesc),]


dAC502 = read.csv(pAC502, sep = "\t", header = TRUE)
rownames(dAC502) = dAC502[,1]
dAC502 = dAC502[,-1]
dAC502 = dAC502[rownames(ddesc),]

vcolor = rep("gray90", dim(ddesc)[1])
vcolor[which(dAC501[,"cell_blue_n"]>= 0)] = "cyan"
vcolor[which(dAC502[,"cell_blue_n"]>= 0)] = "cyan"
vcolor[which(dAC501[,"med_blue_n"]>= 0)] = "blue"
vcolor[which(dAC502[,"med_blue_n"]>= 0)] = "blue"
vcolor[which(dAC501[,"med_green_n"]>= 0)] = "darkgreen"
vcolor[which(dAC502[,"med_green_n"]>= 0)] = "darkgreen"
vcolor[which(dAC501[,"cell_green_n"]>= 0)] = "darkolivegreen3"
vcolor[which(dAC502[,"cell_green_n"]>= 0)] = "darkolivegreen3"
vcolor[which(dAC501[,"med_red_n"]>= 0)] = "red"
vcolor[which(dAC502[,"med_red_n"]>= 0)] = "red"
vcolor[which(dAC501[,"cell_red_n"]>= 0)] = "firebrick"
vcolor[which(dAC502[,"cell_red_n"]>= 0)] = "firebrick"

vpoint = rep(19, dim(ddesc)[1])
vpoint[which(dAC502[,"cell_blue_n"]>= 0)] = 18
vpoint[which(dAC502[,"cell_blue_n"]>= 0 & dAC501[,"cell_blue_n"]>= 0 )] = 8
vpoint[which(dAC502["med_blue_n"]>= 0)] = 18
vpoint[which(dAC502[,"med_blue_n"]>= 0 & dAC501[,"med_blue_n"]>= 0)] = 8
vpoint[which(dAC502[,"med_green_n"]>= 0)] = 18
vpoint[which(dAC502[,"med_green_n"]>= 0 & dAC501[,"med_green_n"]>= 0)] = 8
vpoint[which(dAC502[,"cell_green_n"]>= 0)] = 18
vpoint[which(dAC502[,"cell_green_n"]>= 0 & dAC501[,"cell_green_n"]>= 0)] = 8
vpoint[which(dAC502[,"med_red_n"]>= 0)] = 18
vpoint[which(dAC502[,"med_red_n"]>= 0 & dAC501[,"med_red_n"]>= 0)] = 8
vpoint[which(dAC502[,"cell_red_n"]>= 0)] = 18
vpoint[which(dAC502[,"cell_red_n"]>= 0 & dAC501[,"cell_red_n"]>= 0)] = 8

PCAplot(ddesc, vcolor, vpoint, prout)
