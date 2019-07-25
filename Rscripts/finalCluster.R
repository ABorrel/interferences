#!/usr/bin/env Rscript
source("dendocircular.R")



################
#     MAIN     #
################

args <- commandArgs(TRUE)
pMat = args[1]
pAC50All = args[2]
pcluster = args[3]
colorConsidered = args[4]
methdist = as.character(args[5])
methagg = args[6]
prout = args[7]

#pMat = "/home/borrela2/interference/FP/MACCS-Tanimoto"
#pAC50All = "/home/borrela2/interference/spDataAnalysis/AC50_all"
#pcluster = "/home/borrela2/interference/spDataAnalysis/FinalClustering/IC50/Luc/MACCS-Tanimoto/Luc_IC50_cluster.csv"
#colorConsidered = "Luc_IC50"
#methdist = "None"
#methagg = "ward.D2"
#prout = "/home/borrela2/interference/spDataAnalysis/FinalClustering/IC50/Luc/MACCS-Tanimoto/"


#pMat = "/home/borrela2/interference/FP/MACCS-Tanimoto"
#pAC50All = "/home/borrela2/interference/spDataAnalysis/AC50_all"
#pcluster = "/home/borrela2/interference/spDataAnalysis/FinalClustering/cell_red_n/hepg2/MACCS-Tanimoto/hepg2_cell_red_n_cluster.csv"
#colorConsidered = "hepg2_cell_red_n"
#methdist = "euclidian"
#methagg = "ward.D2"
#prout = "/home/borrela2/interference/spDataAnalysis/FinalClustering/cell_red_n/hepg2/MACCS-Tanimoto/"


dAC50 = read.csv(pAC50All, sep = "\t", header= TRUE)
rownames(dAC50) = dAC50[,1]

dAC50 = dAC50[,c("CASID", colorConsidered)]

dcluster = read.csv(pcluster, sep = ",", header = TRUE)
rownames(dcluster) = dcluster[,1]
dcluster = dcluster[,-1]

if(methdist == "None"){
  dMat = read.csv(pMat, sep = "\t", header = TRUE)
  colnames(dMat) = rownames(dMat)
  #dMat = dMat[1:100,1:100]
  #ddist = as.dist(dMat)
}else{
  dMat = read.csv(pMat, sep = ",", header = TRUE)
  rownames(dMat) = dMat[,1]
  dMat = dMat[,-1]
  #ddist = dist(dMat, method = methdist)
}

dAC50 = dAC50[rownames(dMat),]

# enrichment
nbClust = length(unique(dcluster[,2]))
nbInact = length(which(is.na(dAC50[,2])))
nbAct = length(dAC50[,2]) - nbInact

lpval = NULL
lcount = NULL
lprob = NULL
for (cluster in seq(1, nbClust)){
  dclusttemp = dAC50[which(dcluster[,2] == cluster),2]
  ninactClust = length(which(is.na(dclusttemp)))
  nactClust = length(dclusttemp) - ninactClust
      
  FisherMat = matrix(c(nbAct, nbInact, nactClust, ninactClust),
                      nrow = 2,
                      dimnames = list(Pop = c("Act", "Inact"),
                                      Clust = c("Act", "Inact")))
 
  #print(FisherMat)     
  p = fisher.test(FisherMat, alternative = "greater")
  pval = p$p.value
  lpval = append(lpval, log10(p$p.value))
  lcount = append(lcount, length(dclusttemp))
  lprob = append(lprob, nactClust/(nactClust + ninactClust))
}
tablepval = cbind(seq(1,nbClust), lpval)
tablepval = cbind(tablepval, lcount)
tablepval = cbind(tablepval, lprob)
rownames(tablepval) = tablepval[,1]
colnames(tablepval) = c("Cluster", "Enrichment", "size", "prob")

write.csv(tablepval, paste(prout, "enrichment_cluster.csv", sep = ""))

tablepval = as.data.frame(tablepval)
Cluster = dcluster[rownames(dAC50),2]
Enrichment = NULL
prob = NULL
for(clust in Cluster){
  Enrichment = append(Enrichment, tablepval$Enrichment[which(tablepval$Cluster == clust)])
  prob = append(prob, tablepval$prob[which(tablepval$Cluster == clust)])
}

dAC50 = cbind(dAC50, Enrichment)
dAC50 = cbind(dAC50, prob)

#dAC50 = dAC50[,-2]
dendoClusteringInter(dMat, dAC50, dcluster, prout, methdist, methagg)

