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
pdesc = args[1]
pAC50 = args[2]
methdist = args[3]
methAgg = args[4]
prout = args[5]

#pdesc = "/home/borrela2/interference/spDataAnalysis/DescActive/descClean.csv"
#pAC50 = "/home/borrela2/interference/spDataAnalysis/DescActive/AC50Clean.csv"
#methdist = "euclidean"
#methAgg = "ward.D2"
#prout = "~/interference/spDataAnalysis/DescActive/"


#open ddesc #
#############
ddesc = read.csv(pdesc, sep= ",", header = TRUE)
rownames(ddesc) = ddesc[,1]
ddesc = ddesc[,-1]

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


# for luciferase exclusive
#print(colnames(dAC50))
dIC50 = as.data.frame(dAC50[,c("Luc_IC50")])
rownames(dIC50) = rownames(dAC50)
dendoAtive(ddesc, dIC50, methdist, methAgg, paste(prout, "luc"))
dendoAtive(ddesc, dIC50, methdist, methAgg, paste(prout, "luc"), 1)


# all
dAC50 = dAC50[,ltypeCellChannel]
ddesc = ddesc[rownames(dAC50),]
dendoAtive(ddesc, dAC50, methdist, methAgg, paste(prout, "all"))
dendoAtive(ddesc, dAC50, methdist, methAgg, paste(prout, "all"), 1)
MAXAC50 = max(dAC50, na.rm = TRUE)
MINAC50 = min(dAC50, na.rm = TRUE)


#by color
dblue = dAC50[,ltypeCellBlue]
dblue[dblue == max(dblue, na.rm = TRUE)] = MAXAC50
dblue[dblue == min(dblue, na.rm = TRUE)] = MINAC50
ddescblue = ddesc[rownames(dblue),]
dendoAtive(ddescblue, dblue, methdist, methAgg, paste(prout, "blue"))
dendoAtive(ddescblue, dblue, methdist, methAgg, paste(prout, "blue"), 1)

dgreen = dAC50[,ltypeCellGreen]
dgreen[dgreen == max(dgreen, na.rm = TRUE)] = MAXAC50
dgreen[dgreen == min(dgreen, na.rm = TRUE)] = MINAC50
ddescgreen = ddesc[rownames(dgreen),]
dendoAtive(ddescgreen, dgreen, methdist, methAgg, paste(prout, "green"))
dendoAtive(ddescgreen, dgreen, methdist, methAgg, paste(prout, "green"), 1)

dred = dAC50[,ltypeCellRed]
print(dred)
dred[dred == max(dred, na.rm = TRUE)] = MAXAC50
dred[dred == min(dred, na.rm = TRUE)] = MINAC50
ddescred = ddesc[rownames(dred),]
dendoAtive(ddescred, dred, methdist, methAgg, paste(prout, "red"))
dendoAtive(ddescred, dred, methdist, methAgg, paste(prout, "red"), 1)

## merge color
blue = mergeCol(dAC50, ltypeCellBlue, "blue")
red = mergeCol(dAC50, ltypeCellRed, "red")
green = mergeCol(dAC50, ltypeCellGreen, "green")
dbycol = cbind(blue, green)
dbycol = cbind(dbycol, red)
dendoAtive(ddesc, dbycol, methdist, methAgg, paste(prout, "color"))
dendoAtive(ddesc, dbycol, methdist, methAgg, paste(prout, "color"), 1)
