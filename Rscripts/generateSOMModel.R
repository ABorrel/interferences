#!/usr/bin/env Rscript
require(kohonen)
library(ggplot2)

source("elimcor_sansY.R")

###################
# data management #
###################

openData = function (pfilin, valcor, prout, NbmaxNA=10){
  desc = read.csv (pfilin, header = TRUE, sep = "\t", stringsAsFactors = TRUE)
  
  if(dim(desc)[2] ==1){
    # case of fucking ,
    desc = read.csv (pfilin, header = TRUE, sep = ",", stringsAsFactors = TRUE)
  }
  
  if(length(which(duplicated(desc[,1]))) == 0){
    rownames(desc) = desc[,1]
    desc = desc[,-1]  
  }
  
  
  # deleted col with NA
  lcoldel = NULL
  print (dim(desc))
  for (icol in seq(1, dim(desc)[2])){
    #print (desc[,icol])
    if (sum(is.na(as.vector(desc[,icol]))) > NbmaxNA){
      lcoldel = append(lcoldel, icol)
    }
  }
  
  print("=== colnames deleted ===")
  print(colnames(desc)[lcoldel])
  
  if(is.null(lcoldel)){
    desc = na.omit(desc)
  }else{
    desc = desc[,-lcoldel]  
    desc = na.omit(desc)  
  }
  
  
  
  # dell when sd = 0
  sd_desc = apply (desc[,1:(dim(desc)[2])], 2, sd)
  
  #print (sd_desc)
  #print ("--------")
  sd_0 = which (sd_desc == 0.0)
  
  #print ("------------")
  #print (mode(sd_0))
  #print (length (sd_0))
  #print ("------------")
  if (length(sd_0) != 0){
    #print (as.factor (sd_0))
    #desc = desc[,-sd_0]
    desc=subset(desc,select=-sd_0)
  }
  if (valcor != 0){
    out_elimcor = elimcor_sansY (desc, valcor)
    descriptor = out_elimcor$possetap
    
    descriptor = colnames (desc) [descriptor]
    desc = desc[,descriptor]
    #print (dim(desc))
  }
  
  # again with SD null
  sd_desc =apply (desc[,1:(dim(desc)[2])], 2, sd)
  sd_0 = which (sd_desc == 0)
  if (length(sd_0) != 0){
    desc=subset(desc,select=-sd_0)
  }
  return (list((desc),colnames (desc)))
}


delete.na <- function(DF, n=0) {
  DF[rowSums(is.na(DF)) <= n,]
}

is.integer0 <- function(x)
{
  is.integer(x) && length(x) == 0L
}

eucdist = function(x1, x2, y1, y2){
  distout = sqrt((x1-x2)^2 + (y1-y2)^2)
  return (distout)
}



#################
# Color reverse #
#################

colors <- function(n, alpha = 1) {
  rev(heat.colors(n, alpha))
}

coolBlueHotRed <- function(n, alpha = 1) {rainbow(n, end=4/6, alpha=alpha)[n:1]}

#######
# SOM #
#######

generateSOM = function(ddesc, dAC50, xdim, ydim, prout, modelin){
  write.table(dAC50, paste(prout, "dAC50merged", sep = ""), sep = ",", row.names = FALSE)
  ddesc = as.matrix(scale(ddesc))
  if(is.double(modelin)){
    som_grid <- somgrid(xdim=xdim, ydim=ydim, topo="hexagonal")
    som_model <- som(ddesc, 
                     grid=som_grid, 
                     rlen=100, 
                     alpha=c(0.05,0.01), 
                     keep.data = TRUE)  
  }else{
    som_model = modelin
  }
  
  #write cluster
  dclust = som_model$unit.classif
  names(dclust) = rownames(som_model$data[[1]])
  
  
  dclust = dclust[rownames(dAC50)]# take only active chemical
  dclust = cbind(names(dclust), dclust)
  write.table(dclust, paste(prout, "SOMClust", sep = ""), sep = ",", row.names = FALSE)
  
  lAct = rep(0, xdim*ydim)
  names(lAct) = seq(1, xdim*ydim)
  ltable = table(dclust[,2])
  lAct[names(ltable)] = ltable
  
  ltabinit = table(som_model$unit.classif)
  linitial = rep(0, xdim*ydim)
  names(linitial) = seq(1, xdim*ydim)
  linitial[names(ltabinit)] = ltabinit
  
  #lAct = lAct[names(linitial)]
  #lAct[is.na(lAct)] = 0
  #print (linitial)
  #print(lAct)
  lprob = lAct / linitial
  
  
  #png(paste(prout, "_GlobalCount.png", sep = ""))
  #plot(som_model, type = "counts", palette.name = colors, heatkey = TRUE)
  #dev.off()
  
  #dAC50 = dAC50[rownames(ddesc),]
  #print(dAC50)
  # count by cluster of active 
  #lAct = NULL
  #lprob = NULL
  #nbclust = max(som_model$unit.classif)
  #for(clust in seq(1,nbclust)){
  #  dclust = dAC50[which(som_model$unit.classif == clust),]
  #print(dclust)
  #  ninactClust = length(which(is.na(dclust[,2])))
  #  nactClust = length(dclust[,2]) - ninactClust
  #  lAct = append(lAct, nactClust)
  #  lprob = append(lprob, nactClust/(nactClust + ninactClust))
  #}
  
  #names(lAct) = rownames(ddesc)
  write.csv(lAct, paste(prout, "SOMClustAct", sep = ""))
  
  svg(paste(prout, "CountAC50.svg", sep = ""))
  plot(som_model, type = "property", property=lAct, palette.name=coolBlueHotRed, main = "", dpi=300, height = 20, width = 20, bg = "transparent")
  dev.off()
  
  
  lAct[which(lAct == max(lAct))] = max(table(som_model$unit.classif))# have to calibrate based on the max of the original SOM
  svg(paste(prout, "CountAC50_calibrate.svg", sep = ""))
  plot(som_model, type = "property", property=lAct, palette.name=coolBlueHotRed, main = "", dpi=300, height = 20, width = 20, bg = "transparent")
  dev.off()
  
  if(length(unique((na.omit(lprob)))) != 1){
    
    # prob SOM
    svg(paste(prout, "ProbAC50.svg", sep = ""))
    plot(som_model, type = "property", property=lprob, palette.name=coolBlueHotRed, main = "Prob active", dpi=300, height = 20, width = 20, bg = "transparent")
    dev.off()
    
    write.csv(lprob, paste(prout, "SOMClustActProb", sep = ""))
    
    # enrichment
    if (dAC50 != 0){
      nbInact = length(which(is.na(dAC50)))
      nbAct = dim(dAC50)[1] - nbInact
      lpval = NULL
      nbclust = ydim*xdim
      for(clust in seq(1,nbclust)){
        dclust = dAC50[which(som_model$unit.classif == clust),]
        ninactClust = length(which(is.na(dclust[,2])))
        nactClust = dim(dclust)[1] - ninactClust
        FisherMat = matrix(c(nbAct, nbInact, nactClust, ninactClust),
                           nrow = 2,
                           dimnames = list(Pop = c("Act", "Inact"),
                                           Clust = c("Act", "Inact")))
        p = fisher.test(FisherMat, alternative = "greater")
        pval = p$p.value
        lpval = append(lpval, p$p.value)
      }
      print(length(lpval))
      png(paste(prout, "enrich.png", sep = ""))
      plot(som_model, type = "property", property=log10(lpval), palette.name=coolBlueHotRed, main = "Enrichment log(pvalue(Fisher test))", dpi=300, height = 20, width = 20)
      dev.off()
    }
  }
  return (som_model)
}


generateEnrichCurve = function(ddesc, dAC50, nmax, prout){
  
  xydim = as.integer(sqrt(nmax))
  desc = as.matrix(scale(ddesc))
  
  #Confusion table
  nbInact = length(which(is.na(dAC50)))
  nbAct = dim(dAC50)[1] - nbInact
  
  i = 5
  lM = NULL
  lSD = NULL
  lcount = NULL
  li = NULL
  while (i<=nmax) {
    som_grid <- somgrid(xdim=i, ydim=i, topo="hexagonal")
    
    
    som_model <- som(desc, 
                     grid=som_grid, 
                     rlen=100, 
                     alpha=c(0.05,0.01), 
                     keep.data = TRUE)#,
    li = append(li, i)
    i = i + 2
    lclust = unique(som_model$unit.classif)
    lpval = NULL
    lcountClust = NULL
    for(clust in lclust){
      dclust = dAC50[which(som_model$unit.classif == clust),]
      ninactClust = length(which(is.na(dclust[,2])))
      nactClust = length(dclust[,2]) - ninactClust
      FisherMat = matrix(c(nbAct, nbInact, nactClust, ninactClust),
                         nrow = 2,
                         dimnames = list(Pop = c("Inact", "Act"),
                                         Clust = c("Inact", "Act")))
      #print(FisherMat)
      p = fisher.test(FisherMat, alternative = "greater")
      lpval = append(lpval, p$p.value)
      lcountClust = append(lcountClust, ninactClust + nactClust)
    }
    lM = append(lM, mean(log10(lpval)))
    lSD = append(lSD, sd(log10(lpval)))
    lcount = append(lcount, mean(lcountClust))
  }
  
  p = qplot(li, lM, geom = c("point", "smooth"))+
    labs(x = "Size SOM", y = "Mean log10(pvals Fisher's test)")
  #print(p)
  #Fisher’s exact test
  ggsave(paste(prout, "_EnrichCurve.png", sep = ""), dpi=300, height = 6, width = 6)
  
  
  p = qplot(li, lcount, geom = c("point", "smooth"))+
    labs(x = "Size SOM", y = "Mean size cluster")
  #print(p)
  #Fisher’s exact test
  ggsave(paste(prout, "meanSizeCluster.png", sep = ""), dpi=300, height = 6, width = 6)
  
}



generateSOMGlobal = function(ddesc, xdim, ydim, prout, modelin){
  
  pmodel = paste(prout, "SOMmodel.Rdata", sep = "")
  if(file.exists(pmodel)){
    load(pmodel)
    return(som_model) 
  }
  
  ddesc = as.matrix(scale(ddesc))
  som_grid <- somgrid(xdim=xdim, ydim=ydim, topo="hexagonal")
  som_model <- som(ddesc, 
                    grid=som_grid, 
                    rlen=100, 
                    alpha=c(0.05,0.01), 
                    keep.data = TRUE)  
  
  svg(paste(prout, "SOM_count.svg", sep = ""))
  plot(som_model, type = "count", palette.name=coolBlueHotRed, main = "", dpi=300, height = 20, width = 20, bg = "transparent")
  dev.off()
  
  
  #write cluster
  dclust = som_model$unit.classif
  names(dclust) = rownames(ddesc)
  write.csv(dclust, paste(prout, "SOMClust", sep = ""))
  save(som_model, file = pmodel)
  
  return (som_model)
}


specificSOM = function(techno, cutoff, modelG, przscore, prout){
  
  prout = paste(prout, techno, "_" , cutoff, "/", sep = "")
  dir.create(prout)
  
  
  pzscore = paste(przscore, "Zscore_", techno, ".csv", sep = "")
  dzscore = read.csv(pzscore, sep = ",", header = TRUE)
  rownames(dzscore) = dzscore[,1]
  dzscore = dzscore[,-1]
                
  dZselect = dzscore[which(dzscore[,1] > cutoff),]
  if(dim(dZselect)[1] != 0){
    generateSOM(ddesc, dZselect, size, size, prout, modelG)
  }
}



################
#     MAIN     #
################

#args <- commandArgs(TRUE)
#pdesc = args[1]
#prout = args[2]
#size = args[3]
#pmodel = args[4]

pdesc = "./data/tableDesc1D2D"
prout = "./results/SOM/"
size = 15



##############################
# Process descriptors matrix #
##############################
lddesc = openData(pdesc, 0.9, 20)
ddesc = lddesc[[1]]

# => process to the SOM
modelG = generateSOMGlobal(ddesc, size, size, prout, model)


# for technology
cutoff = 5
przscore = "./results/technology/"

techno = "Fluorescence"
specificSOM(techno, cutoff, modelG, przscore, prout)

techno = "Radiometry"
specificSOM(techno, cutoff, modelG, przscore, prout)

techno = "Luminescence"
specificSOM(techno, cutoff, modelG, przscore, prout)

techno = "Spectrophotometry"
specificSOM(techno, cutoff, modelG, przscore, prout)

techno = "Microscopy"
specificSOM(techno, cutoff, modelG, przscore, prout)



# for platform
cutoff = 5
przscore = "./results/platform/"

techno = "APR"
specificSOM(techno, cutoff, modelG, przscore, prout)

techno = "ATG"
specificSOM(techno, cutoff, modelG, przscore, prout)

techno = "BSK"
specificSOM(techno, cutoff, modelG, przscore, prout)

techno = "NVS"
specificSOM(techno, cutoff, modelG, przscore, prout)

techno = "TOX21"
specificSOM(techno, cutoff, modelG, przscore, prout)

techno = "CEETOX"
specificSOM(techno, cutoff, modelG, przscore, prout)

techno = "CLD"
specificSOM(techno, cutoff, modelG, przscore, prout)

