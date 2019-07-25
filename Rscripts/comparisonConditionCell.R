#!/usr/bin/env Rscript
source("~/development/Stat/source/testStats.R")




mergeCol = function(ddesc, lmerge, delNA = 1){
  
  dtemp = ddesc[,lmerge]
  if(delNA == 1){
    dout = rowMeans(dtemp, na.rm = TRUE)
  }else{
    dout = rowMeans(dtemp, na.rm = FALSE)  
  }
  
  dout[which(dout == "NaN")] = NA
  
  return (dout)
}




################
#     MAIN     #
################

args <- commandArgs(TRUE)
pAC50 = args[1]
presult = args[2]


pAC50 = "/home/borrela2/interference/spDataAnalysis/AC50_all"
presult = "/home/borrela2/interference/spDataAnalysis/TtestAC50/"

# AC50
dAC50 = read.csv(pAC50, sep="\t", header = TRUE)
rownames(dAC50) = dAC50[,1]
dAC50 = dAC50[,-1]


# compute list of chem active for any condition and cell
lcrossBlue = mergeCol(dAC50, c("hepg2_cell_blue_n", "hepg2_med_blue_n", "hek293_cell_blue_n",  "hek293_med_blue_n"), delNA = 0)
lcrossGreen = mergeCol(dAC50, c("hepg2_cell_green_n", "hepg2_med_green_n", "hek293_cell_green_n",  "hek293_med_green_n"), delNA = 0)
lcrossRed = mergeCol(dAC50, c("hepg2_cell_red_n", "hepg2_med_red_n", "hek293_cell_red_n",  "hek293_med_red_n"), delNA = 0)

dcross = cbind(lcrossBlue, lcrossGreen)
dcross = cbind(dcross, lcrossRed)
colnames(dcross) = c("blue", "green", "red")

lBlue = mergeCol(dAC50, c("hepg2_cell_blue_n", "hepg2_med_blue_n", "hek293_cell_blue_n",  "hek293_med_blue_n"), delNA = 1)
lGreen = mergeCol(dAC50, c("hepg2_cell_green_n", "hepg2_med_green_n", "hek293_cell_green_n",  "hek293_med_green_n"), delNA = 1)
lRed = mergeCol(dAC50, c("hepg2_cell_red_n", "hepg2_med_red_n", "hek293_cell_red_n",  "hek293_med_red_n"), delNA = 1)

dnotcross = cbind(lBlue, lGreen)
dnotcross = cbind(dnotcross, lRed)
colnames(dnotcross) = c("blue", "green", "red")

lcross = mergeCol(dcross, c("blue", "green", "red"), delNA = 1)
lnocross = mergeCol(dnotcross, c("blue", "green", "red"), delNA = 1)

#print(lcross)
ldel = names(lcross)[which(!is.na(lcross))]
#print (ldel)

lnocross[ldel] = NA

dout = computeComparison(lnocross, lcross)
print (dout)
pfilout = paste(presult, "TtestCrossactive.csv", sep = "")
write.csv(dout, pfilout)


# by color #
############
# BLUE
ldel = names(lcrossBlue)[which(!is.na(lcrossBlue))]
lBlue[ldel] = NA
dout = computeComparison(lBlue, lcrossBlue)
print (dout)
pfilout = paste(presult, "TtestCrossactiveBlue.csv", sep = "")
write.csv(dout, pfilout)

# GREEN
ldel = names(lcrossGreen)[which(!is.na(lcrossGreen))]
lGreen[ldel] = NA
dout = computeComparison(lGreen, lcrossGreen)
print (dout)
pfilout = paste(presult, "TtestCrossactiveGreen.csv", sep = "")
write.csv(dout, pfilout)

# RED
ldel = names(lcrossRed)[which(!is.na(lcrossRed))]
lRed[ldel] = NA
dout = computeComparison(lRed, lcrossBlue)
print (dout)
pfilout = paste(presult, "TtestCrossactiveRed.csv", sep = "")
write.csv(dout, pfilout)



# by cell type #
################
# blue
lBlueFree = mergeCol(dAC50, c("hepg2_cell_blue_n", "hek293_cell_blue_n"), delNA = 1)
lBluecell = mergeCol(dAC50, c("hepg2_med_blue_n", "hek293_med_blue_n"), delNA = 1)

dout = computeComparison(lBlueFree, lBluecell)
print (dout)
pfilout = paste(presult, "TtestCellTypeBlue.csv", sep = "")
write.csv(dout, pfilout)

# green
lGreenFree = mergeCol(dAC50, c("hepg2_cell_green_n", "hek293_cell_green_n"), delNA = 1)
lGreencell = mergeCol(dAC50, c("hepg2_med_green_n", "hek293_med_green_n"), delNA = 1)

dout = computeComparison(lGreenFree, lGreencell)
print (dout)
pfilout = paste(presult, "TtestCellTypeGreen.csv", sep = "")
write.csv(dout, pfilout)


# red
lRedFree = mergeCol(dAC50, c("hepg2_cell_red_n", "hek293_cell_red_n"), delNA = 1)
lRedcell = mergeCol(dAC50, c("hepg2_med_red_n", "hek293_med_red_n"), delNA = 1)

dout = computeComparison(lRedFree, lRedcell)
print (dout)
pfilout = paste(presult, "TtestCellTypeRed.csv", sep = "")
write.csv(dout, pfilout)