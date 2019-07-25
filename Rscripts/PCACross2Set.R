#!/usr/bin/env Rscript
library(ggplot2)
#source("~/development/Rglobal/source/PCAdrawer.R")

library(factoextra)



generatePCAcoords = function(din){
  
  dinScale = scale(din)
  data.cor=cor(dinScale)
  data.eigen=eigen(data.cor)
  lambda = data.eigen$values
  var_cap = lambda/sum(lambda)*100
  cp = data.eigen$vectors
  rownames (cp) = colnames (dinScale)
  colnames (cp) = colnames (dinScale)
  data_plot = as.matrix(dinScale)%*%cp
  
  return(list(data_plot, var_cap, cp))
}



################
#     MAIN     #
################

args <- commandArgs(TRUE)
pdesc1 = args[1]
pAC50 = args[2]
pdesc2 = args[3]
prout = args[4]

#pdesc1 = "/home/borrela2/interference/spDataAnalysis/CrossPCA/descClean.csv"
#pAC50 = "/home/borrela2/interference/spDataAnalysis/AC50_all"
#pdesc2 = "/home/borrela2/interference/testing/588_resorufin/descMat"
#prout = "/home/borrela2/interference/testing/588_resorufin/PCA/"

# for window dev
#pdesc1 = "./../../trash/descClean.csv"
#pAC50 = "./../../trash/AC50_all"
#pdesc2 = "./../../trash/descMat"
#prout = "./../../trash/"


ddesc1 = read.csv(pdesc1, sep = ",", row.names = 1)
ddesc2 = read.csv(pdesc2, sep = ",", row.names = 1)

ddesc2_aff = ddesc2$Aff
ddesc2 = ddesc2[-which(colnames(ddesc2) == "Aff")]

res.pca <- prcomp(ddesc1, scale = TRUE)
var_cap = generatePCAcoords(ddesc1)[[2]]

add_coord = predict(res.pca, newdata = ddesc2)

# color based on AC50
dAC50 = read.csv(pAC50, sep = "\t", header = TRUE)
rownames(dAC50) = dAC50[,1]
dAC50 = dAC50[,-1]

dAC50 = dAC50[rownames(ddesc1),]



vcolor = rep("darkgray", dim(ddesc1)[1])
vcolor[which(dAC50[,"hepg2_cell_blue_n"]>= 0)] = "cyan"
vcolor[which(dAC50[,"hek293_cell_blue_n"]>= 0)] = "cyan"
vcolor[which(dAC50[,"hepg2_med_blue_n"]>= 0)] = "blue"
vcolor[which(dAC50[,"hek293_med_blue_n"]>= 0)] = "blue"
vcolor[which(dAC50[,"hepg2_med_green_n"]>= 0)] = "darkgreen"
vcolor[which(dAC50[,"hek293_med_green_n"]>= 0)] = "darkgreen"
vcolor[which(dAC50[,"hepg2_cell_green_n"]>= 0)] = "darkolivegreen3"
vcolor[which(dAC50[,"hek293_cell_green_n"]>= 0)] = "darkolivegreen3"
vcolor[which(dAC50[,"hepg2_med_red_n"]>= 0)] = "red"
vcolor[which(dAC50[,"hek293_med_red_n"]>= 0)] = "red"
vcolor[which(dAC50[,"hepg2_cell_red_n"]>= 0)] = "firebrick"
vcolor[which(dAC50[,"hek293_cell_red_n"]>= 0)] = "firebrick"


# vcolor
vcoladd = rep("bisque", dim(add_coord)[1])
vcoladd[which(ddesc2_aff ==1)] = "black"


png (paste (prout, "PCA_color.png", sep = ""), 1700, 1500)
par(mar=c(8,8,8,8))
plot(add_coord[,1], add_coord[,2], pch=19, col = "bisque", xlab = paste("CP1: ", round (var_cap[1], 2), "%", sep = ""), ylab = paste("CP2: ", round(var_cap[2], 2), "%", sep = ""), cex.lab = 4, cex.main = 4, cex.axis = 1.75, cex = 4)
points(res.pca$x[,1], res.pca$x[,2], pch=19, col = vcolor, cex = 2.5)
points(res.pca$x[which(vcolor != "darkgray"),1], res.pca$x[which(vcolor != "darkgray"),2], pch=19, col = vcolor[which(vcolor != "darkgray")], cex = 2.5)
points(add_coord[which(vcoladd != "bisque"),1], add_coord[which(vcoladd != "bisque"),2], col = vcoladd[which(vcoladd != "bisque")], pch = 19, cex = 2.5)
abline(h=0,v=0)
warnings ()
dev.off()

#points(data_plot[which(vcolor != "gray90"),1], data_plot[which(vcolor != "gray90"),2], pch=20, col = vcolor[which(vcolor != "gray90")], cex = 4)
#plot(res.pca$x[,1], res.pca$x[,2], xlim = c(-40, 20), ylim = c(-20, 20))

#points(add_coord[,1], add_coord[,2], col = "blue")
#points(res.pca$x[,1], res.pca$x[,2], col = "red")
