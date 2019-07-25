#!/usr/bin/env Rscript
source("~/development/Rglobal/source/dataManager.R")
source("clustering.R")
source("dendocircular.R")
source("radialPlot.R")



################
#     MAIN     #
################

args <- commandArgs(TRUE)
pdesc = args[1]
pAC50 = args[2]
pcluster = args[3]
prresult = args[4]

#pdesc = "/home/borrela2/interference/spDataAnalysis/clusters/hclust_euc_wardD2_gap_stat/descClean.csv"
#pAC50 = "/home/borrela2/interference/spDataAnalysis/tox21-luc-biochem-p1/AC50_sample"
#pcluster = "/home/borrela2/interference/spDataAnalysis/clusters/hclust_euc_wardD2_gap_stat/clusterMain.csv"
#prresult = "/home/borrela2/interference/spDataAnalysis/tox21-luc-biochem-p1/hclust_euc_wardD2_gap_stat/"



#open ddesc #
#############
ddesc = read.csv(pdesc, header = TRUE)
rownames(ddesc) = ddesc[,1]
ddesc = ddesc[,-1]

#check SD null for desc, case of secondar clustering
ddesc = delSDNull(ddesc)

#open dIC50 #
#############
dIC50 = read.csv(pAC50, header = TRUE, sep = "\t")
rownames(dIC50) = dIC50[,1]
# add line for global
if(dim(dIC50)[2]==2){
  dIC50 = cbind(dIC50, dIC50[,2])
}
# transform to log
dIC50 = dIC50[,-1]
dIC50 = -log10(dIC50)

print(head(dIC50))

dcluster = read.csv(pcluster, header = TRUE, sep = ",")
print(dcluster)
rownames(dcluster) = dcluster[,1]
#dcluster = dcluster[,-1]
colnames(dcluster)[2] = "cluster"

#print(dcluster)
#print("dddd")

# Means and SD by cluster #
###########################


affByCluster(dIC50, dcluster, prresult)


# dendogram cluster #
#####################
dendogramCluster(ddesc, dIC50, dcluster, prresult)


# radial plot by cluster #
##########################

lcluster = unique(dcluster[,2])
prRadialplots = paste(prresult, "Radial/", sep = "")
dir.create(prRadialplots)

for(clust in lcluster){
  dtemp = dIC50[which(dcluster[,2] == clust),]
  radialByCluster(dtemp, paste(prRadialplots, clust, ".svg", sep = ""))
}



