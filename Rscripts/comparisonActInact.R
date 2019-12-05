#!/usr/bin/env Rscript
source("~/development/Stat/source/testStats.R")

computeTableSignif = function(ddesc, dAC50, prresult){
  
  for(i in seq(1, dim(dAC50)[2])){
    lact = rownames(dAC50[which(dAC50[,i] != "NA"),])
    linact = rownames(dAC50[is.na(dAC50[,i] == "NA"),])
    
    dout = data.frame()
    for(j in seq(1, dim(ddesc)[2])){
      dact = ddesc[lact,j]
      dinact = ddesc[linact,j]
      typeTest = conditionTtest(dact, dinact)
      if (typeTest == 0){
        pval = comparisonTest(dact, dinact, "non-parametric")
      }else{
        pval = comparisonTest(dact, dinact, "parametric")
      }
      dout[j,1] = length(lact)
      dout[j,2] = length(linact)
      dout[j,3] = round(mean(dact, na.rm = TRUE),2)
      dout[j,4] = round(mean(dinact, na.rm = TRUE),2)
      dout[j,5] = pval
      dout[j,6] = signifPvalue(pval)  
      
    }
    rownames(dout) = colnames(ddesc)
    colnames(dout) = c("Nbact", "Nbinact", "Mact", "Minact", "pval", "Signif")
    orderPval = order(dout[,5],decreasing=F)
    dout = dout[orderPval,]
    pfilout = paste(presult, colnames(dAC50)[i], ".csv", sep = "")
    write.csv(dout, pfilout)
  }
}




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
pdesc = args[1]
pAC50 = args[2]
presult = args[3]

#pdesc = "/home/borrela2/interference/spDataAnalysis/Ttest/descClean.csv"
#pAC50 = "/home/borrela2/interference/spDataAnalysis/AC50_all"
#presult = "/home/borrela2/interference/spDataAnalysis/Ttest/"

# AC50
dAC50 = read.csv(pAC50, sep="\t", header = TRUE)
rownames(dAC50) = dAC50[,1]
dAC50 = dAC50[,-1]

# descriptor
ddesc = read.csv(pdesc, sep = ",", header = TRUE)
rownames(ddesc) = ddesc[,1]
ddesc = ddesc[,-1]

# for any conditions
computeTableSignif(ddesc, dAC50, prresult)


ltypeCellBlue = c("hepg2_cell_blue_n", "hepg2_med_blue_n", "hek293_cell_blue_n",  "hek293_med_blue_n") 
ltypeCellGreen = c("hepg2_cell_green_n", "hepg2_med_green_n", "hek293_cell_green_n",  "hek293_med_green_n") 
ltypeCellRed = c("hepg2_cell_red_n", "hepg2_med_red_n", "hek293_cell_red_n",  "hek293_med_red_n") 


lcrossBlue = mergeCol(dAC50, c("hepg2_cell_blue_n", "hepg2_med_blue_n", "hek293_cell_blue_n",  "hek293_med_blue_n"), delNA = 0)
lcrossGreen = mergeCol(dAC50, c("hepg2_cell_green_n", "hepg2_med_green_n", "hek293_cell_green_n",  "hek293_med_green_n"), delNA = 0)
lcrossRed = mergeCol(dAC50, c("hepg2_cell_red_n", "hepg2_med_red_n", "hek293_cell_red_n",  "hek293_med_red_n"), delNA = 0)

dcross = cbind(lcrossBlue, lcrossGreen)
dcross = cbind(dcross, lcrossRed)
colnames(dcross) = c("blue", "green", "red")
dcross = cbind(dcross, mergeCol(dcross, c("blue", "green", "red")))
colnames(dcross) = c("allblue", "allgreen", "allred", "cross")
computeTableSignif(ddesc, dcross, prresult)


dAC50color = cbind(mergeCol(dAC50, ltypeCellBlue), mergeCol(dAC50, ltypeCellGreen))
dAC50color = cbind(dAC50color, mergeCol(dAC50, ltypeCellRed))
colnames(dAC50color) = c("blue", "green", "red")

# for color
computeTableSignif(ddesc, dAC50color, prresult)
