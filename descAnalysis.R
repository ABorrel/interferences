#!/usr/bin/env Rscript
source ("tool.R")
source("cardMatrix.R")
source("PCAplot.R")
source("dendocircular.R")
source("clustering.R")
source("distributions.R")

################
#     MAIN     #
################

args <- commandArgs(TRUE)
pdesc = args[1]
pdata = args[2] #to take affinity
prout = args[3]
valcor = as.double(args[4])
maxQuantile = as.double(args[5])
#logaff = as.integer(args[6])
#plotPCA = as.integer(args[7])
#corMatrix = as.integer(args[8])
#histplot = as.integer(args[9])
#circularDendo = as.integer(args[10])
#optimal_clustering = as.integer(args[11])


#pdesc = "/home/borrela2/interference/luc-biochem/test/desc-test.csv"
#pdata = "/home/borrela2/interference/luc-biochem/test/LD50.csv"
#prout = "/home/borrela2/interference/luc-biochem/test/"

plotPCA = 1
corMatrix = 1
histplot = 1
circularDendo = 1
#valcor = 0.80
#maxQuantile = 85
logaff = 1
optimal_clustering = 0


# Process descriptors matrix #
##############################
dglobal = openData(pdesc, valcor, prout, c(1,2))
dglobal = dglobal[[1]]


rownames(dglobal) = dglobal[,1]
dglobal = dglobal[,-1]
dglobal = dglobal[,-1]

# remove descriptor with a distribution on one quantile
dglobal = delnohomogeniousdistribution(dglobal, maxQuantile)

#######################
# order with affinity #
#######################
# Opening
daffinity = read.csv(pdata, sep = "\t", header = TRUE)
# fusion 
daffinity = mergeIdenticRow(daffinity)

rownames(daffinity) = daffinity[,1]
daffinity = daffinity[,-1]
# remove class + remove NA
#daffinity = daffinity[,-2]

print(dim(daffinity))



if(logaff == 1){
  histAff(daffinity, prout)
  daffinity[,1] = -log10(daffinity[,1])
  histlogAff(daffinity, prout)
}

# remove no selected affinity or bad quality from manual curation
print(dim(dglobal))
lID = intersect(rownames(daffinity), rownames(dglobal))

dglobal = dglobal[lID,]
daffinity = daffinity[lID,]

write.csv(dglobal, paste(prout, "globalTable.csv", sep = ""))

#ord = NULL
#for(i in seq(1,dim(dglobal)[1])){
#  ipos = which(daffinity[,1] == rownames(dglobal)[i])
#  ord = append (ord,ipos)
#}

#daffinity = daffinity[ord,]


#orderaff = order(daffinity[,2],decreasing=T)
#daffinity = daffinity[orderaff,]


#dglobal = dglobal[orderaff,]
#print(rownames(dglobal))
#print(colnames(dglobal))

#print(dim(dglobal))
#print(dim(daffinity))

#print(rownames(dglobal))
#print(rownames(daffinity))



if (corMatrix == 1){
  cardMatrixCor(cor(dglobal), paste(prout, "matrixCor_", valcor, sep = ""), 6)
  
  # matrix cor with pMIC
  lcor = NULL
  ldesc = colnames(dglobal)
  for(desc in ldesc){
    lvaldesc = dglobal[,desc]
    corval = cor(daffinity[,1], lvaldesc, use = "pairwise.complete.obs")
    lcor = append(lcor, corval)
  }
  names(lcor) = ldesc
  
  write.csv(lcor, file = paste(prout, "corAC50-Desc.csv"))
}


if(plotPCA == 1){
  PCAplot(dglobal, daffinity, prout)
}


if (histplot == 1){
  histDataOne(data1 = dglobal, paste(prout, "histDesc_", valcor, ".pdf", sep = ""))
}

if (circularDendo == 1){
  dendogramCircle(dglobal, daffinity, prout)
}


