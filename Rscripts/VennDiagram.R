#!/usr/bin/env Rscript
library(VennDiagram)


venPlot = function(xin, prout){
  
  colfill = rainbow(length(xin))
  
  venn.diagram(xin, 
               filename = paste(prout, ".tiff", sep = ""),
               col = "black",
               lty = "dotted",
               lwd = 4,
               fill = colfill,
               alpha = 0.50,
               cex = 1.5,
               fontfamily = "serif",
               fontface = "bold",
               cat.col = colfill,
               cat.cex = 1,
               cat.fontfamily = "serif")
}



################
#     MAIN     #
################

args <- commandArgs(TRUE)
pAC50 = args[1]
prout = args[2]

#pAC50 = "/home/borrela2/interference/spDataAnalysis/tox21-spec-hek293-p1/AC50_sample_curve"
#prout = "/home/borrela2/interference/spDataAnalysis/tox21-spec-hek293-p1/venn/"

dAC50 = read.csv(pAC50, sep = "\t", header = TRUE)
rownames(dAC50) = dAC50[,1]
dAC50 = dAC50[,-1]


#for green
x = list()
x$cell_green = which(dAC50[,c("cell_green")] >= 0)
x$med_green = which(dAC50[,c("med_green")] >= 0)
venPlot(x, paste(prout, "green", sep = ""))

#for green n
x = list()
x$cell_green_n = which(dAC50[,c("cell_green_n")] >= 0)
x$med_green_n = which(dAC50[,c("med_green_n")] >= 0)
venPlot(x, paste(prout, "green_n", sep = ""))

#for blue
x = list()
x$med_blue = which(dAC50[,c("med_blue")] >= 0)
x$cell_blue = which(dAC50[,c("cell_blue")] >= 0)
venPlot(x, paste(prout, "blue", sep = ""))

#for blue n
x = list()
x$cell_blue_n = which(dAC50[,c("cell_blue_n")] >= 0)
x$med_blue_n = which(dAC50[,c("med_blue_n")] >= 0)
venPlot(x, paste(prout, "bluen", sep = ""))


#for red
x = list()
x$med_red = which(dAC50[,c("med_red")] >= 0)
x$cell_red = which(dAC50[,c("cell_red")] >= 0)
venPlot(x, paste(prout, "red", sep = ""))

#for red n
x = list()
x$med_red = which(dAC50[,c("med_red_n")] >= 0)
x$cell_red = which(dAC50[,c("cell_red_n")] >= 0)
venPlot(x, paste(prout, "red_n", sep = ""))



# all blue
x = list()
x$med_blue = which(dAC50[,c("med_blue")] >= 0)
x$cell_blue = which(dAC50[,c("cell_blue")] >= 0)
x$cell_blue_n = which(dAC50[,c("cell_blue_n")] >= 0)
x$med_blue_n = which(dAC50[,c("med_blue_n")] >= 0)
venPlot(x, paste(prout, "allblue", sep = ""))

# all red
x = list()
x$med_red = which(dAC50[,c("med_red")] >= 0)
x$cell_red = which(dAC50[,c("cell_red")] >= 0)
x$cell_red_n = which(dAC50[,c("cell_red_n")] >= 0)
x$med_red_n = which(dAC50[,c("med_red_n")] >= 0)
venPlot(x, paste(prout, "allred", sep = ""))

#all green
x = list()
x$med_green = which(dAC50[,c("med_green")] >= 0)
x$cell_green = which(dAC50[,c("cell_green")] >= 0)
x$cell_green_n = which(dAC50[,c("cell_green_n")] >= 0)
x$med_green_n = which(dAC50[,c("med_green_n")] >= 0)
venPlot(x, paste(prout, "allgreen", sep = ""))


#for cell free
x = list()
x$med_red = which(dAC50[,c("med_red")] >= 0)
x$med_green = which(dAC50[,c("med_green")] >= 0)
x$med_blue = which(dAC50[,c("med_blue")] >= 0)
venPlot(x, paste(prout, "cellFree", sep = ""))

x = list()
x$med_red_n = which(dAC50[,c("med_red_n")] >= 0)
x$med_green_n = which(dAC50[,c("med_green_n")] >= 0)
x$med_blue_n = which(dAC50[,c("med_blue_n")] >= 0)
venPlot(x, paste(prout, "cellFree_n", sep = ""))


#for based cell
x = list()
x$cell_red = which(dAC50[,c("cell_red")] >= 0)
x$cell_green = which(dAC50[,c("cell_green")] >= 0)
x$cell_blue = which(dAC50[,c("cell_blue")] >= 0)
venPlot(x, paste(prout, "cellbased", sep = ""))

x$cell_red_n = which(dAC50[,c("cell_red_n")] >= 0)
x$cell_green_n = which(dAC50[,c("cell_green_n")] >= 0)
x$cell_blue_n = which(dAC50[,c("cell_blue_n")] >= 0)
venPlot(x, paste(prout, "cellbased_n", sep = ""))

