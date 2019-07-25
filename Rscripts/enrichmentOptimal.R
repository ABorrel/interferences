#!/usr/bin/env Rscript
library(factoextra)
library(RootsExtremaInflections)


computeEnrichmentOpt = function(lclust, dAC50, AC50name, prout){
  
  denrich = NULL
  dAC50 = dAC50[,AC50name]
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
        dclusttemp = dAC50[which(dclust == cluster)]
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
  write.csv(denrich, paste(prout, AC50name, ".csv", sep = ""))
    
  denrich = as.data.frame(denrich)
  p = ggplot(denrich, aes(NBclust, Menrich, color = "Menrich")) + 
    #geom_point(size=1.5, colour="black", shape=21)+
    geom_line()+
    geom_line(aes(NBclust, SDenrich, color = "SD"), denrich)+
    #geom_point(aes(NBclust, SDenrich, color = "SD"), size=1.5, colour="black", shape=21)+
    labs(x = "Number cluster", y = paste("Enrichment"))
  ggsave(paste(prout, AC50name, "_enrich.png", sep = ""), width = 8, height = 8, dpi = 300)
    
  p = ggplot(denrich, aes(NBclust, Mcount, color = "M count")) + 
    geom_line()+
    geom_line(aes(NBclust, SDcount, color = "SD"), denrich)+
    labs(x = "Number cluster", y = paste("Nb chemicals"))
  ggsave(paste(prout, AC50name, "_size.png", sep = ""), width = 8, height = 8, dpi = 300)
  
  nbclust = denrich$NBclust[which(denrich$Menrich == max(denrich$Menrich))]
  
  dclustbest = lclust[[nbclust]]
  lcpt = names(dclustbest)
  dw = cbind(lcpt, dclustbest)
  colnames(dw) = c("CAS", "Cluster")
  
  write.csv(dw,paste(prout, AC50name, "_cluster.csv", sep = ""))
  
  return(nbclust)
  
}




################
#     MAIN     #
################

args <- commandArgs(TRUE)
pMat = args[1]
pAC50All = args[2]
pprelimEnrichment = args[3]
prresult = args[4]
methclustering = args[5]
methdist = args[6]
methagg = args[7]


#pMat = "/home/borrela2/interference/FP/FPMol-Tanimoto"
#pAC50All = "/home/borrela2/interference/spDataAnalysis/AC50_all"
#prresult = "/home/borrela2/interference/spDataAnalysis/FPclusters/FPMol-Tanimoto/optclustering/"
#pprelimEnrichment = "/home/borrela2/interference/spDataAnalysis/FPclusters/FPMol-Tanimoto/hek293_med_green_n.csv"
#methclustering = "hclust"
#methdist = "None"
#methagg = "ward.D2"

#pMat = "/home/borrela2/interference/spDataAnalysis/Descclusters/enrichmentIndex/descClean.csv"
#pAC50All = "/home/borrela2/interference/spDataAnalysis/AC50_all"
#prresult = "/home/borrela2/interference/spDataAnalysis/Descclusters/enrichmentIndex/"
#methclustering = "hclust"
#methdist = "euclidean"
#methagg = "ward.D2"

dAC50 = read.csv(pAC50All, sep = "\t", header= TRUE)
rownames(dAC50) = dAC50[,1]


if(methdist == "None"){
  dMat = read.csv(pMat, sep = "\t", header = TRUE)
  colnames(dMat) = rownames(dMat)
  #dMat = dMat[1:100,1:100]
  ddist = as.dist(dMat)
}else{
  dMat = read.csv(pMat, sep = ",", header = TRUE)
  rownames(dMat) = dMat[,1]
  dMat = dMat[,-1]
  ddist = dist(dMat, method = methdist)
}

# define interval best clustering
dtable = read.csv(pprelimEnrichment, sep = ",", header = TRUE )

AC50considered = strsplit(pprelimEnrichment, "/")
AC50considered = AC50considered[[1]][length(AC50considered[[1]])]
AC50considered = strsplit(AC50considered, "[.]")
AC50considered = AC50considered[[1]][1]


nbclusteropt = dtable$NBclust[which(dtable$Menrich == max(dtable$Menrich))]

# case of nb cluster is the max of the last 
if (nbclusteropt > (max(dtable$NBclust) - 100)){
  d = inflexi(dtable$NBclust,as.double(dtable$Menrich),1,length(dtable$NBclust),10,10,plots=FALSE)
  nboptimal = d$finfl[1]
  print(nboptimal)

  imin = nboptimal - 50
  imax = nboptimal
  
}else{
  imin = nbclusteropt - 50
  if(imin < 0){
    imin = 2
  }
  imax = nbclusteropt + 50
}

i = imin
lclust = list()
hc <- hclust(ddist, method = methagg)
while (i < imax){
  outhcut = cutree(hc, k = i)
  lclust[[i]] = outhcut
  i = i + 1
}
#print (lclust)

nbopt = computeEnrichmentOpt(lclust, dAC50, AC50considered, prresult)

