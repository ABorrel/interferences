#!/usr/bin/env Rscript
library(factoextra)


computeEnrichment = function(lclust, dAC50All, prout){
  
  lAC50 = colnames(dAC50All)
  lAC50 = lAC50[-1]
  
  # reduce just for only normalized data
  lAC50 = c("Luc_IC50", "hepg2_med_blue_n", "hepg2_med_green_n", "hepg2_med_red_n", "hek293_med_blue_n", "hek293_med_green_n",
            "hek293_med_red_n", "hepg2_cell_blue_n", "hepg2_cell_green_n", "hepg2_cell_red_n", "hek293_cell_blue_n", "hek293_cell_green_n",
            "hek293_cell_red_n")  

  for (AC50 in lAC50){
    denrich = NULL
    dAC50 = dAC50All[,AC50]
    nbInact = length(which(is.na(dAC50)))
    nbAct = length(dAC50) - nbInact
    
    #print(AC50)
    #print(nbAct)
    #print(nbInact)
    
    for(dclust in lclust){
      if(!is.null(dclust)){
        nbClust = length(unique(dclust))
        lpval = NULL
        lcount = NULL
        for (cluster in seq(1, nbClust)){
          dclusttemp = dAC50All[which(dclust == cluster),AC50]
          ninactClust = length(which(is.na(dclusttemp)))
          nactClust = length(dclusttemp) - ninactClust
          
          FisherMat = matrix(c(nbAct, nbInact, nactClust, ninactClust),
                             nrow = 2,
                             dimnames = list(Pop = c("Act", "Inact"),
                                             Clust = c("Act", "Inact")))
          
          p = fisher.test(FisherMat, alternative = "greater")
          pval = p$p.value
          lpval = append(lpval, log10(p$p.value))
          lcount = append(lcount, length(dclusttemp))
        }
        enrich = c(nbClust, mean(lpval), sd(lpval), mean(lcount), sd(lcount))
        denrich = rbind(denrich, enrich)
      }
    }
    colnames(denrich) = c("NBclust", "Menrich", "SDenrich", "Mcount", "SDcount")
    write.csv(denrich, paste(prout, AC50, ".csv", sep = ""))
    
    denrich = as.data.frame(denrich)
    p = ggplot(denrich, aes(NBclust, Menrich, color = "Menrich")) + 
      #geom_point(size=1.5, colour="black", shape=21)+
      geom_line()+
      geom_line(aes(NBclust, SDenrich, color = "SD"), denrich)+
      #geom_point(aes(NBclust, SDenrich, color = "SD"), size=1.5, colour="black", shape=21)+
      labs(x = "Number cluster", y = paste("Enrichment"))
    ggsave(paste(prout, AC50, "_enrich.png", sep = ""), width = 8, height = 8, dpi = 300)
    
    p = ggplot(denrich, aes(NBclust, Mcount, color = "M count")) + 
      geom_line()+
      ylim(0,50)+
      geom_line(aes(NBclust, SDcount, color = "SD"), denrich)+
      labs(x = "Number cluster", y = paste("Nb chemicals"))
    ggsave(paste(prout, AC50, "_size.png", sep = ""), width = 8, height = 8, dpi = 300)
    
  }
}



################
#     MAIN     #
################

args <- commandArgs(TRUE)
pMat = args[1]
pAC50All = args[2]
prresult = args[3]
methclustering = args[4]
methdist = args[5]
methagg = args[6]

#pMat = "/home/borrela2/interference/FP/FPMol-Tanimoto"
#pAC50All = "/home/borrela2/interference/spDataAnalysis/AC50_all"
#prresult = "/home/borrela2/interference/spDataAnalysis/FPclusters/FPMol-Tanimoto/"
#methclustering = "hclust"
#methdist = "None"
#methagg = "ward.D2"

#pMat = "/home/borrela2/interference/spDataAnalysis/Descclusters/enrich-index/descClean.csv"
#pAC50All = "/home/borrela2/interference/spDataAnalysis/AC50_all"
#prresult = "/home/borrela2/interference/spDataAnalysis/Descclusters/enrich-index/"
#methclustering = "hclust"
#methdist = "euclidean"
#methagg = "ward.D2"

dAC50 = read.csv(pAC50All, sep = "\t", header= TRUE)
rownames(dAC50) = dAC50[,1]


if(methdist == "None"){
  dMat = read.csv(pMat, sep = "\t", header = TRUE)
  colnames(dMat) = rownames(dMat)
  dMat[is.na(dMat)] = 1
  #dMat = dMat[1:100,1:100]
  ddist = as.dist(dMat)
}else{
  dMat = read.csv(pMat, sep = ",", header = TRUE)
  rownames(dMat) = dMat[,1]
  dMat = dMat[,-1]
  ddist = dist(dMat, method = methdist)
}

imin = 2
imax = dim(dMat)[1] -1
# limit to 2000 clusters
if(imax >= 2000){
  imax = 2000
}
by = 10
i = imin
lclust = list()
hc <- hclust(ddist, method = methagg)
while (i < imax){
  outhcut = cutree(hc, k = i)
  lclust[[i]] = outhcut
  i = i + by
}
#print (lclust)
computeEnrichment(lclust, dAC50, prresult)
