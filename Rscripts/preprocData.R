#!/usr/bin/env Rscript
source ("~/development/Rglobal/source/dataManager.R")


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
nbNA = as.integer(args[7])

#pdesc = "/home/borrela2/interference/spDataAnalysis/SOMactiveluc/descActive"
#pdata = "/home/borrela2/interference/spDataAnalysis/SOMactiveluc/AC50Active"
#prout = "/home/borrela2/interference/spDataAnalysis/SOMactiveluc/"
#logAff = 0
  
#pdesc = "/home/borrela2/interference/spDataAnalysis/Desc/tableDesc1D2D"
#pdata = "/home/borrela2/interference/spDataAnalysis/tox21-luc-biochem-p1/AC50_sample"
#prout = "/home/borrela2/interference/spDataAnalysis/tox21-luc-biochem-p1/Stat/"
#valcor = 0.9
#maxquantile = 90
#logaff = 1


##############################
# Process descriptors matrix #
##############################

dglobal = openData(pdesc, valcor, prout,  nbNA)
dglobal = dglobal[[1]]

print("==== Preprocessing ====")
print(paste("Data initial: dim = ", dim(dglobal)[1], dim(dglobal)[2], sep = " "))

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
  if(dim(daffinity)[2] == 2){
    daffinity = cbind(daffinity, daffinity[,2])
    colnames(daffinity)[3] = "Aff"
  }
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
  if(!is.null(ldel)){
    daffinity = daffinity[-ldel,]
  }
  

  # merge with data descriptors and remove data remove from the manual curation
  lID = intersect(rownames(daffinity), rownames(dglobal))
  dglobal = dglobal[lID,]
  daffinity = daffinity[lID,]


  # Write table 
  paffout = paste(prout, "AC50Clean.csv", sep = "")
  write.csv(daffinity, paffout, col.names = TRUE, row.names = TRUE)
}
  
pdesout = paste(prout, "descClean.csv", sep = "")
write.csv(dglobal, pdesout, col.names = TRUE, row.names = TRUE)


