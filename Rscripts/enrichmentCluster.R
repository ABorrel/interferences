#!/usr/bin/env Rscript
library(ggplot2)


################
#     MAIN     #
################

args <- commandArgs(TRUE)
pAC50 = args[1]
pCluster = args[2]
prout = args[3]

#pAC50 = "/home/borrela2/interference/spDataAnalysis/AC50_all"
#pCluster = "/home/borrela2/interference/spDataAnalysis/Descclusters/hclust_euc_wardD2_silhouette/cluster.csv"
#prout = "/home/borrela2/interference/spDataAnalysis/Descclusters/hclust_euc_wardD2_silhouette/"

dAC50 = read.csv(pAC50, sep = "\t", header = TRUE)
rownames(dAC50) = dAC50[,1]

dClust = read.csv(pCluster, sep = ",", header = TRUE)
rownames(dClust) = dClust[,1]

#same size data
lcas = intersect(rownames(dAC50), rownames(dClust))
dAC50 = dAC50[lcas,]
dClust = dClust[lcas,]

#cluster
lcluster = unique(dClust[,2])
lAC50 = colnames(dAC50)
lAC50 = lAC50[-1]

dout = NULL
for(cluster in lcluster){
  print(cluster)
  lpval = NULL
  for (AC50 in lAC50){
    print(AC50)
    #confusion table #
    ##################
    # complete pop
    dtemp = dAC50[,AC50]
    print(dtemp)
    nbInact = length(which(is.na(dtemp)))
    nbAct = length(dtemp) - nbInact
    
    # for cluster
    dclusttemp = dAC50[which(dClust[,"cluster"] == cluster),AC50]
    ninactClust = length(which(is.na(dclusttemp)))
    nactClust = length(dclusttemp) - ninactClust
    FisherMat = matrix(c(nbAct, nbInact, nactClust, ninactClust),
                       nrow = 2,
                       dimnames = list(Pop = c("Act", "Inact"),
                                       Clust = c("Act", "Inact")))
    p = fisher.test(FisherMat, alternative = "less")
    pval = p$p.value
    lpval = append(lpval, log10(p$p.value))
  }
  dout = rbind(dout, lpval)  
}
colnames(dout) = lAC50
rownames(dout) = lcluster

print(dout)

write.csv(dout, paste(prout, "Enrichment", sep = ""))