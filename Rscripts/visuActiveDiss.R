#!/usr/bin/env Rscript
source("dendocircular.R")




mergeCol = function(ddesc, lmerge, namecol){

  dtemp = ddesc[,lmerge]
  dout = rowMeans(dtemp, na.rm = TRUE)
  dout[which(dout == "NaN")] = NA
  
  return (dout)
}


################
#     MAIN     #
################

args <- commandArgs(TRUE)
pdiss = args[1]
pAC50 = args[2]
methAgg = args[3]
prout = args[4]


#pdiss = "~/interference/FP/MACCS-Tanimoto"
#pAC50 = "~/interference/spDataAnalysis/DescActive/AC50Clean.csv"
#methAgg = "ward.D2"
#prout = "~/interference/spDataAnalysis/DescActive/"


prout = paste(prout, basename(pdiss), sep = "")


#open ddesc #
#############
ddiss = read.csv(pdiss, sep= "\t", header = TRUE)
colnames(ddiss) = rownames(ddiss)

#open dIC50 #
#############
dAC50 = read.csv(pAC50, header = TRUE, sep = ",")
rownames(dAC50) = dAC50[,1]

# AC50 considered
ltypeCellChannel = c("hepg2_cell_blue_n", "hepg2_med_blue_n", "hek293_cell_blue_n",  "hek293_med_blue_n", 
                    "hepg2_cell_green_n", "hepg2_med_green_n", "hek293_cell_green_n", "hek293_med_green_n",
                    "hepg2_cell_red_n", "hepg2_med_red_n", "hek293_cell_red_n",  "hek293_med_red_n")


ltypeCellBlue = c("hepg2_cell_blue_n", "hepg2_med_blue_n", "hek293_cell_blue_n",  "hek293_med_blue_n") 
ltypeCellGreen = c("hepg2_cell_green_n", "hepg2_med_green_n", "hek293_cell_green_n",  "hek293_med_green_n") 
ltypeCellRed = c("hepg2_cell_red_n", "hepg2_med_red_n", "hek293_cell_red_n",  "hek293_med_red_n") 


# all
dAC50 = dAC50[,ltypeCellChannel]
ddiss = ddiss[rownames(dAC50), rownames(dAC50)]
dendoAtiveDiss(ddiss, dAC50, methAgg, prout)

# => not used
#by color
#dblue = dAC50[,ltypeCellBlue]
#ddescblue = ddesc[rownames(dblue),]
#dendoAtive(ddescblue, dblue, methdist, methAgg, paste(prout, "blue"))

#dgreen = dAC50[,ltypeCellGreen]
#ddescgreen = ddesc[rownames(dgreen),]
#dendoAtive(ddescgreen, dgreen, methdist, methAgg, paste(prout, "green"))

#dred = dAC50[,ltypeCellRed]
#ddescred = ddesc[rownames(dred),]
#dendoAtive(ddescred, dred, methdist, methAgg, paste(prout, "red"))



## merge color
#blue = mergeCol(dAC50, ltypeCellBlue, "blue")
#red = mergeCol(dAC50, ltypeCellRed, "red")
#green = mergeCol(dAC50, ltypeCellGreen, "green")
#dbycol = cbind(blue, green)
#dbycol = cbind(dbycol, red)
#dendoAtive(ddesc, dbycol, methdist, methAgg, paste(prout, "color"))

