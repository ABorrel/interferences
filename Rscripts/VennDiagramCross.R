#!/usr/bin/env Rscript
library(VennDiagram)


venPlot = function(xin, colfill, prout){
  
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
pluc = args[1]
phepg2 = args[2]
phek293 = args[3]
prout = args[4]

#pluc = "/home/borrela2/interference/spDataAnalysis/tox21-luc-biochem-p1/AC50_combine"
#phepg2 = "/home/borrela2/interference/spDataAnalysis/tox21-spec-hepg2-p1/AC50_sample_curve"
#phek293 = "/home/borrela2/interference/spDataAnalysis/tox21-spec-hek293-p1/AC50_sample_curve"
#prout = "/home/borrela2/interference/spDataAnalysis/CrossVenn/"

dluc = read.csv(pluc, sep = "\t", header = TRUE)
dhepg2 = read.csv(phepg2, sep = "\t", header = TRUE)
dhek293 = read.csv(phek293, sep = "\t", header = TRUE)

rownames(dhepg2) = dhepg2[,1]
rownames(dhek293) = dhek293[,1]
rownames(dluc) = dluc[,1]
dhek293 = dhek293[rownames(dhepg2),]
dluc = dluc[rownames(dhepg2),]


#cell free hepg2
#x = list()
#x$hepg2_blue = which(dhepg2[,"med_blue_n"]>= 0)
#x$hepg2_green = which(dhepg2[,"med_green_n"]>= 0)
#x$hepg2_red = which(dhepg2[,"med_red_n"]>= 0)
#venPlot(x, c("blue", "green", "red"), paste(prout, "cell_free_hepg2", sep = ""))


#cell free hek293
#x = list()
#x$hek293_blue = which(dhek293[,"med_blue_n"]>= 0)
#x$hek293_green = which(dhek293[,"med_green_n"]>= 0)
#x$hek293_red = which(dhek293[,"med_red_n"]>= 0)
#venPlot(x, c("blue", "green", "red"), paste(prout, "cell_free_hek293", sep = ""))


# cell based hepg2
#x = list()
#x$hepg2_blue = which(dhepg2[,"cell_blue_n"]>= 0)
#x$hepg2_green = which(dhepg2[,"cell_green_n"]>= 0)
#x$hepg2_red = which(dhepg2[,"cell_red_n"]>= 0)
#venPlot(x, c("blue", "green", "red"), paste(prout, "cell_based_hepg2", sep = ""))


#cell free hek293
#x = list()
#x$hek293_blue = which(dhek293[,"cell_blue_n"]>= 0)
#x$hek293_green = which(dhek293[,"cell_green_n"]>= 0)
#x$hek293_red = which(dhek293[,"cell_red_n"]>= 0)
#venPlot(x, c("blue", "green", "red"), paste(prout, "cell_based_hek293", sep = ""))



#med_red
#x = list()
#x$luc = which(dluc[, "IC50"]>= 0)
#x$hepg2 = which(dhepg2[,"med_red"]>= 0)
#x$hek293 = which(dhek293[,"med_red"]>= 0)
#venPlot(x, c("grey", "red", "red"), paste(prout, "med_red", sep = ""))

#"med_blue"
#x = list()
#x$luc = which(dluc[, "IC50"]>= 0)
#x$hepg2 = which(dhepg2[,"med_blue"]>= 0)
#x$hek293 = which(dhek293[,"med_blue"]>= 0)
#venPlot(x, c("grey", "blue", "blue"), paste(prout, "med_blue", sep = ""))

#"med_green"
#x = list()
#x$luc = which(dluc[, "IC50"]>= 0)
#x$hepg2 = which(dhepg2[,"med_green"]>= 0)
#x$hek293 = which(dhek293[,"med_green"]>= 0)
#venPlot(x, c("grey", "green", "green"), paste(prout, "med_green", sep = ""))

#"cell_green"
#x = list()
#x$luc = which(dluc[, "IC50"]>= 0)
#x$hepg2 = which(dhepg2[,"cell_green"]>= 0)
#x$hek293 = which(dhek293[,"cell_green"]>= 0)
#venPlot(x, c("grey", "green", "green"), paste(prout, "cell_green", sep = ""))

#"cell_blue"
#x = list()
#x$luc = which(dluc[, "IC50"]>= 0)
#x$hepg2 = which(dhepg2[,"cell_blue"]>= 0)
#x$hek293 = which(dhek293[,"cell_blue"]>= 0)
#venPlot(x, c("grey", "blue", "blue"), paste(prout, "cell_blue", sep = ""))

#"cell_red"
#x = list()
#x$luc = which(dluc[, "IC50"]>= 0)
#x$hepg2 = which(dhepg2[,"cell_red"]>= 0)
#x$hek293 = which(dhek293[,"cell_red"]>= 0)
#venPlot(x, c("grey", "red", "red"), paste(prout, "cell_red", sep = ""))

#"cell_blue_n"
#x = list()
#x$luc = which(dluc[, "IC50"]>= 0)
#x$hepg2 = which(dhepg2[,"cell_blue_n"]>= 0)
#x$hek293 = which(dhek293[,"cell_blue_n"]>= 0)
#venPlot(x, c("grey", "blue", "blue"), paste(prout, "cell_blue_n", sep = ""))

#"med_blue_n
#x = list()
#x$luc = which(dluc[, "IC50"]>= 0)
#x$hepg2 = which(dhepg2[,"med_blue_n"]>= 0)
#x$hek293 = which(dhek293[,"med_blue_n"]>= 0)
#venPlot(x, c("grey", "blue", "blue"), paste(prout, "med_blue_n", sep = ""))


#"cell_green_n"
#x = list()
#x$luc = which(dluc[, "IC50"]>= 0)
#x$hepg2 = which(dhepg2[,"cell_green_n"]>= 0)
#x$hek293 = which(dhek293[,"cell_green_n"]>= 0)
#venPlot(x, c("grey", "green", "green"), paste(prout, "cell_green_n", sep = ""))

#"med_green_n
#x = list()
#x$luc = which(dluc[, "IC50"]>= 0)
#x$hepg2 = which(dhepg2[,"med_green_n"]>= 0)
#x$hek293 = which(dhek293[,"med_green_n"]>= 0)
#venPlot(x, c("grey", "green", "green"), paste(prout, "med_green_n", sep = ""))


#"cell_red_n"
#x = list()
#x$luc = which(dluc[, "IC50"]>= 0)
#x$hepg2 = which(dhepg2[,"cell_red_n"]>= 0)
#x$hek293 = which(dhek293[,"cell_red_n"]>= 0)
#venPlot(x, c("grey", "red", "red"), paste(prout, "cell_red_n", sep = ""))

#"med_red_n
#x = list()
#x$luc = which(dluc[, "IC50"]>= 0)
#x$hepg2 = which(dhepg2[,"med_red_n"]>= 0)
#x$hek293 = which(dhek293[,"med_red_n"]>= 0)
#venPlot(x, c("grey", "red", "red"), paste(prout, "med_red_n", sep = ""))

# cross four ways #
###################
#blue normalised
x = list()
x$hepg2_free = which(!is.na(dhepg2[,"med_blue_n"]))
x$hek293_free = which(!is.na(dhek293[,"med_blue_n"]))
x$hepg2_cell = which(!is.na(dhepg2[,"cell_blue_n"]))
x$hek293_cell = which(!is.na(dhek293[,"cell_blue_n"]))
venPlot(x, c("blue", "blue", "cyan", "cyan"), paste(prout, "blue_n", sep = ""))


#green normalised
x = list()
x$hepg2_free = which(dhepg2[,"med_green_n"]>= 0)
x$hek293_free = which(dhek293[,"med_green_n"]>= 0)
x$hepg2_cell = which(dhepg2[,"cell_green_n"]>= 0)
x$hek293_cell = which(dhek293[,"cell_green_n"]>= 0)
venPlot(x, c("darkgreen", "darkgreen", "darkolivegreen3", "darkolivegreen3"), paste(prout, "green_n", sep = ""))


# red normalised
x = list()
x$hepg2_free = which(dhepg2[,"med_red_n"]>= 0)
x$hek293_free = which(dhek293[,"med_red_n"]>= 0)
x$hepg2_cell = which(dhepg2[,"cell_red_n"]>= 0)
x$hek293_cell = which(dhek293[,"cell_red_n"]>= 0)
venPlot(x, c("red", "red", "firebrick", "firebrick"), paste(prout, "red_n", sep = ""))

# color
lID = intersect(rownames(dhek293), rownames(dhepg2))
dhek293 = dhek293[lID,]
dhepg2 = dhepg2[lID,]


#luc
lIDluc = intersect(rownames(dluc), lID)
dluc = dluc[lIDluc,]
liluc = which(dluc[,"IC50"]>= 0)


liblue = which(dhek293[,"cell_blue_n"]>= 0)
liblue = append(liblue, which(dhepg2[,"cell_blue_n"]>= 0))
liblue = append(liblue, which(dhek293[,"med_blue_n"]>= 0))
liblue = append(liblue,  which(dhepg2[,"med_blue_n"]>= 0))
liblue = unique(liblue)


ligreen = which(dhek293[,"cell_green_n"]>= 0)
ligreen = append(ligreen, which(dhepg2[,"cell_green_n"]>= 0))
ligreen = append(ligreen, which(dhek293[,"med_green_n"]>= 0))
ligreen = append(ligreen,  which(dhepg2[,"med_green_n"]>= 0))
ligreen = unique(ligreen)


lired = which(dhek293[,"cell_red_n"]>= 0)
lired = append(lired, which(dhepg2[,"cell_red_n"]>= 0))
lired = append(lired, which(dhek293[,"med_red_n"]>= 0))
lired = append(lired,  which(dhepg2[,"med_red_n"]>= 0))
lired = unique(lired)


x = list()
x$blue = liblue
x$green = ligreen
x$red = lired
venPlot(x, c("blue", "green", "red"), paste(prout, "color", sep = ""))


x = list()
x$luciferase = liluc
x$blue = liblue
x$green = ligreen
x$red = lired
venPlot(x, c("grey", "blue", "green", "red"), paste(prout, "color_luc", sep = ""))

