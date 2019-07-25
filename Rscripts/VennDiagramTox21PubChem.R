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
pTox21SMI = args[1]
pTox21AC50 = args[2]
pPubChemSMI = args[3]
pPubChemAC50 = args[4]
prout = args[5]

#pTox21SMI = "/home/borrela2/interference/VennTox21Pubmed/tox21SMI.csv"
#pPubChemSMI = "/home/borrela2/interference/VennTox21Pubmed/PubChemSMI.csv"
#pTox21AC50 = "/home/borrela2/interference/spDataAnalysis/AC50_all"
#pPubChemAC50 = "/home/borrela2/interference/testing/588_resorufin/tableSmi.csv"
#prout = "/home/borrela2/interference/VennTox21Pubmed/"

dAC50Tox21 = read.csv(pTox21AC50, sep = "\t", header = TRUE)
rownames(dAC50Tox21) = dAC50Tox21[,1] 
dSMItox21 = read.csv(pTox21SMI, sep = "\t", header = TRUE)
rownames(dSMItox21) = dSMItox21[,1]

dAC50PubChem = read.csv(pPubChemAC50, sep = "\t", header = TRUE)  
rownames(dAC50PubChem) = dAC50PubChem[,1]

dSMIPubChem = read.csv(pPubChemSMI, sep = "\t", header = TRUE)
rownames(dSMIPubChem) = dSMIPubChem[,1]

xglobal = list()
xglobal$PubChem = dSMIPubChem[,2]
xglobal$tox21 = dSMItox21[,2]
venPlot(xglobal, c("red", "blue"), paste(prout, "Tox21VSPubChem", sep = ""))


lcelltype = c("hepg2_cell_blue_n","hepg2_cell_green_n", "hepg2_cell_red_n", "hek293_cell_blue_n", "hek293_cell_green_n", "hek293_cell_red_n", "Luc_IC50")

for (celltype in lcelltype){
  x = list()
  x$PubChem_active = dAC50PubChem[which(dAC50PubChem[,3] == "Active"),1]
  x$PubChem_active = dSMIPubChem[x$PubChem_active,2]
  
  x$PubChem_inactive = dAC50PubChem[which(dAC50PubChem[,3] == "Inactive"),1]
  x$PubChem_inactive = dSMIPubChem[x$PubChem_inactive,2]
  
  # remove intersection active VS inactive
  linter = intersect(x$PubChem_active, x$PubChem_inactive)
  for(inter in linter){
    x$PubChem_active = x$PubChem_active[-which(inter == x$PubChem_active)]
    x$PubChem_inactive = x$PubChem_inactive[-which(inter == x$PubChem_inactive)]
  }
  
  x$Tox21_active = dAC50Tox21[which(dAC50Tox21[,celltype] != "NA"), 1]
  x$Tox21_active = dSMItox21[x$Tox21_active,2]
  x$Tox21_active = na.omit(x$Tox21_active)
  
  x$Tox21_inactive = dAC50Tox21[which(is.na(dAC50Tox21[,celltype])), 1]
  x$Tox21_inactive = dSMItox21[x$Tox21_inactive,2]
  x$Tox21_inactive = na.omit(x$Tox21_inactive)
  
  linter = intersect(x$Tox21_active, x$Tox21_inactive)
  for(inter in linter){
    x$Tox21_active = x$Tox21_active[-which(inter == x$Tox21_active)]
    x$Tox21_inactive = x$Tox21_inactive[-which(inter == x$Tox21_inactive)]
  }
  venPlot(x, c("blue", "blue", "red", "red"), paste(prout, celltype, sep = ""))
}




