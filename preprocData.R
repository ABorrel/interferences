#!/usr/bin/env Rscript
source ("tool.R")
source("dataManager.R")


################
#     MAIN     #
################

args <- commandArgs(TRUE)
pdesc = args[1]
pdata = args[2] #to take IC50 or class
valcor = args[3]
maxquantile = as.double(args[4])
logAff = args[5]
prout = args[6]

#pdesc = "/home/borrela2/interference/spDataAnalysis/Desc/tableDesc1D2D"
#pdata = "/home/borrela2/interference/spDataAnalysis/tox21-luc-biochem-p1/AC50_sample"
#prout = "/home/borrela2/interference/spDataAnalysis/tox21-luc-biochem-p1/Stat/"
#valcor = 0.9
#maxquantile = 90
#logaff = 1


##############################
# Process descriptors matrix #
##############################

dglobal = openData(pdesc, valcor, prout, c(1))
dglobal = dglobal[[1]]

print("==== Preprocessing ====")
print(paste("Data initial: dim = ", dim(dglobal)[1], dim(dglobal)[2], sep = " "))

rownames(dglobal) = dglobal[,1]
dglobal = dglobal[,-1]

##########
# filter #
##########

dglobal = delnohomogeniousdistribution(dglobal, maxquantile)
print(paste("Data after filtering: dim = ", dim(dglobal)[1], dim(dglobal)[2], sep = " "))

#######################
# order with affinity #
#######################
# Opening
if(pdata != "0"){
  daffinity = read.csv(pdata, sep = "\t", header = TRUE)
  rownames(daffinity) = daffinity[,1]
  daffinity = daffinity[,-1]

  # transform #
  if(logAff == 1){
    daffinity = -log10(daffinity)
  }

  # remove line of NA
  i = 1
  ldel = NULL
  while (i < dim(daffinity)[1]) {
    if(length(which(is.na(daffinity[i,]))) == dim(daffinity)[2]){
      ldel = append(ldel, i)
    }
    i = i + 1
  }
  daffinity = daffinity[-ldel,]

  # merge with data descriptors and remove data remove from the manual curation
  lID = intersect(rownames(daffinity), rownames(dglobal))
  dglobal = dglobal[lID,]
  daffinity = daffinity[lID,]


  # Write table 
  paffout = paste(prout, "IC50Clean.csv", sep = "")
  write.csv(daffinity, paffout, col.names = TRUE, row.names = TRUE)
}
  
pdesout = paste(prout, "descClean.csv", sep = "")
write.csv(dglobal, pdesout, col.names = TRUE, row.names = TRUE)


