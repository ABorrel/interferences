#!/usr/bin/env Rscript
library(factoextra)

predictClosestCompoundCluster = function(dglobal, distMeth, aggMeth, nbcluster){
  
  dglobal = scale(dglobal)
  ddist = dist(dglobal, method = distMeth)
  
  outclust = hcut(ddist, k = nbcluster, hc_func = "hclust", hc_method = aggMeth)
  clusTest = outclust$cluster["test"]
  
  outcluster = outclust$cluster[which(outclust$cluster == clusTest)]
  
  outcluster = outcluster[-which(names(outcluster)=="test")]
  
  return (outcluster)
}




################
#     MAIN     #
################

args <- commandArgs(TRUE)
pDescAll = args[1]
pDescPred = args[2]
penrichment = args[3]
pcluster = args[4]
distMeth = args[5]
aggMeth = args[6]

#pDescAll = "/home/borrela2/interference/spDataAnalysis/descClean.csv"
#pDescPred = "/home/borrela2/interference/spDataAnalysis/predictions/test/test.txt"
#penrichment = "/home/borrela2/interference/spDataAnalysis/FinalClustering/med_red_n/hek293/Desc/enrichment_cluster.csv"
#pcluster = "/home/borrela2/interference/spDataAnalysis/FinalClustering/med_red_n/hek293/Desc/hek293_med_red_n_cluster.csv"
#distMeth = "euclidian"
#aggMeth = "ward.D2"

#pDescAll = "/home/borrela2/interference/spDataAnalysis/descClean.csv"
#pDescPred = "/home/borrela2/interference/spDataAnalysis/predictions/50-69-1/50-69-1.txt"
#penrichment = "/home/borrela2/interference/spDataAnalysis/FinalClustering/med_red_n/hepg2/Desc-euclidean-ward.D2/enrichment_cluster.csv"
#pcluster = "/home/borrela2/interference/spDataAnalysis/FinalClustering/med_red_n/hepg2/Desc-euclidean-ward.D2/hepg2_med_red_n_cluster.csv"
#distMeth = "euclidian"
#aggMeth = "ward.D2"

dDescAll = read.csv(pDescAll, sep = ",", header = TRUE)
rownames(dDescAll) = dDescAll[,1]
colnames(dDescAll)[1] = "CAS"

dDescPred = read.csv(pDescPred, sep = "\t", header = TRUE)
rownames(dDescPred) = "test"

denrich = read.csv(penrichment, sep = ",", header = TRUE)
rownames(denrich) = denrich[,1]

dcluster = read.csv(pcluster, sep = ",", header = TRUE)
rownames(dcluster) = dcluster[,1]


nbcluster = max(dcluster$Cluster)
dtopredict = dDescPred[,colnames(dDescAll)]
dglobal = rbind(dDescAll, dtopredict)
dglobal = dglobal[,-1]


lCAScluster = predictClosestCompoundCluster(dglobal, distMeth, aggMeth, nbcluster)


lclust = dcluster[lCAScluster, "Cluster"]
lenrich = denrich[lclust, "Enrichment"]
print(mean(lenrich))