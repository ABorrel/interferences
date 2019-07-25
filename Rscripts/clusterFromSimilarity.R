#!/usr/bin/env Rscript
source("dendocircular.R")
source("radialPlot.R")
source("clustering.R")

################
#     MAIN     #
################

args <- commandArgs(TRUE)
pMat = args[1]
clusterType = args[2]
aggregtype = args[3]
optimalCluster = args[4]
prresult = args[5]

#aggregtype = "ward.D2"#"ward.D2", "complete", "single", "average"
#clusterType = "hclust"#"hclust", "kmeans"
#optimalCluster = "silhouette"#, "wss", "gap_stat"
#pMat = "~/interference/FP/FPMol-Tanimoto"
#prresult = "/home/borrela2/interference/spDataAnalysis/FPclusters/FPMol-Tanimoto-hclust_None_wardD2_silhouette/"

# verbose #
###########
print(pMat)
print(prresult)
print(aggregtype)
print(clusterType)
print(optimalCluster)


###################
# Process  matrix #
###################
dmat = read.csv(pMat, header = TRUE, sep = "\t")
colnames(dmat) = rownames(dmat)

# transform to distance
pcluster = optimalClutersDist(dmat, prresult, clusterType, optimalCluster, aggregtype)
