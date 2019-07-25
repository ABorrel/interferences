#!/usr/bin/env Rscript
source("dendocircular.R")
source("radialPlot.R")
source("clustering.R")
source("~/development/Rglobal/source/dataManager.R")

################
#     MAIN     #
################

args <- commandArgs(TRUE)
pdesc = args[1]
paff = args[2] #to take affinity
pcluster = args[3]
prresult = args[4]
disttype = args[5]
clusterType = args[6]
aggregtype = args[7]
optimalCluster = args[8]


#pdesc = "/home/borrela2/interference/spDataAnalysis/tox21-luc-biochem-p1/Stat/descClean.csv"
#paff = "/home/borrela2/interference/spDataAnalysis/tox21-luc-biochem-p1/Stat/IC50.csv"
#pcluster = "0"
#prresult = "/home/borrela2/interference/spDataAnalysis/tox21-luc-biochem-p1/Stat/clustering/"

#disttype = "euc"
#aggregtype = "ward.D2"#"ward.D2", "complete", "single", "average"
#clusterType = "hclust"#"hclust", "kmeans"
#optimalCluster = "gap_stat"#"silhouette", "wss", "gap_stat"

# verbose #
###########
print(pdesc)
print(paff)
print(pcluster)
print(prresult)
print(disttype)
print(aggregtype)
print(clusterType)
print(optimalCluster)


##############################
# Process descriptors matrix #
##############################
ddesc = read.csv(pdesc, header = TRUE)
rownames(ddesc) = ddesc[,1]
ddesc = ddesc[,-1]

#check SD null for desc, case of secondar clustering
ddesc = delSDNull(ddesc)

if (pcluster == "0"){
  pcluster = optimalCluters(ddesc, prresult, clusterType, optimalCluster, aggregtype)
}

if (paff != 0){
  
  dIC50 = read.csv(paff, header = TRUE)
  rownames(dIC50) = dIC50[,1]
  dIC50 = dIC50[,-1]
  
  dcluster = read.csv(pcluster, header = TRUE, sep = ",")
  rownames(dcluster) = dcluster[,1]
  #dcluster = dcluster[,-1]
  
  dIC50 = read.csv(paff, header = TRUE)
  rownames(dIC50) = dIC50[,1]
  dIC50 = dIC50[,-1]
  
  # M and SD by cluster
  affByCluster(dIC50, dcluster, prresult)
  
  
  # dendogram cluster #
  #####################
  
  dendogramCluster(ddesc, dIC50, dcluster, prresult)
  
  
  # radial plot by cluster #
  ##########################
  
  #lcluster = unique(dcluster[,2])
  prRadialplots = paste(prresult, "Radial/", sep = "")
  dir.create(prRadialplots)
  
  for(clust in lclust){
    dtemp = dIC50[which(dcluster[,2] == clust),]
    radialByCluster(dtemp, paste(prRadialplots, clust, ".svg", sep = ""))
  }
}
