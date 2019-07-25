#!/usr/bin/env Rscript
require(kohonen)
library(ggplot2)



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


mergeAC50 = function(dAC50, lcol){
  
  lCAS = rownames(dAC50)
  imax = dim(dAC50)[1]
  i = 1
  lAC50 = NULL
  while(i <= imax){
    lac50temp = as.double(dAC50[i,lcol]) 
    AC50 = mean(lac50temp, na.rm = TRUE)
    lAC50 = append(lAC50, AC50)
    i = i + 1
  }
  lAC50[lAC50 == "NaN"] = NA
  dout = cbind(lCAS, lAC50)
  rownames(dout) = lCAS
  dout = na.omit(dout)
  return(as.data.frame(dout))
}


################
#     MAIN     #
################

args <- commandArgs(TRUE)
pdesc = args[1]
pAC50 = args[2]
pmodel = args[3]
prout = args[4]

  
#pdesc = "/home/borrela2/interference/spDataAnalysis/SOMactiveluc/descClean.csv"
#pAC50 = "/home/borrela2/interference/spDataAnalysis/SOMactiveluc/AC50Clean.csv"
#pmodel = "/home/borrela2/interference/spDataAnalysis/SOMactiveluc/SOMmodel.Rdata"
#prout = "/home/borrela2/interference/spDataAnalysis/SOMactiveluc/"


#pdesc = "/home/borrela2/interference/spDataAnalysis/SOMactiveluc/descClean.csv"
#pAC50 = "/home/borrela2/interference/spDataAnalysis/SOMactiveluc/AC50Clean.csv"
#pmodel = "/home/borrela2/interference/spDataAnalysis/SOM/SOMmodel.Rdata"
#prout = "/home/borrela2/interference/spDataAnalysis/SOMactiveluc/all/"


#pmodel = "/home/borrela2/interference/spDataAnalysis/SOMactive/SOMmodel.Rdata"
#pdesc = "/home/borrela2/interference/spDataAnalysis/SOMactive/descClean.csv"
#pAC50 = "/home/borrela2/interference/spDataAnalysis/SOMactive/AC50Clean.csv"
#prout = "/home/borrela2/interference/spDataAnalysis/SOMactive/all/"

#pdesc = "/home/borrela2/interference/ToxCast_analysis/SOM/descClean.csv"
#pAC50 = "0"
#prout = "/home/borrela2/interference/ToxCast_analysis/SOM/"

#pdesc = "/home/borrela2/interference/spDataAnalysis/SOMactiveluc/descClean.csv"
#pAC50 = "/home/borrela2/interference/spDataAnalysis/SOMactiveluc/AC50Clean.csv"
#pmodel = "0"
#prout = "/home/borrela2/interference/spDataAnalysis/SOMactiveluc/"

if (pmodel != "0"){
  load(pmodel)  
}else{
  model = as.double(pmodel)
}


##############################
# Process descriptors matrix #
##############################
ddesc = read.csv(pdesc, header = TRUE)
rownames(ddesc) = ddesc[,1]
ddesc = ddesc[,-1]
lchannel = c("hek293_med_blue_n", "hek293_cell_blue_n", "hek293_med_green_n", "hek293_cell_green_n", "hek293_med_red_n", "hek293_cell_red_n", "hepg2_med_blue_n", "hepg2_cell_blue_n","hepg2_med_green_n", "hepg2_cell_green_n","hepg2_med_red_n", "hepg2_cell_red_n")


if(pAC50 == "0"){
  model = generateSOM(ddesc, pAC50, 15, 15, prout, model)
}else{
  dAC50 = read.csv(pAC50, sep = "\t", header = TRUE)
  # control format -> open to coma
  if (dim(dAC50)[2] == 1){
    dAC50 = read.csv(pAC50, sep = ",", header = TRUE)
  }
  
  rownames(dAC50) = dAC50[,1]
  lAC50 = colnames(dAC50)
  lAC50 = lAC50[-1]
  
  # merge all aff
  # case only luc
  if(dim(dAC50)[2] == 3){
    lchannel = c("Luc_IC50")
    dAC50temp = mergeAC50(dAC50, lchannel)
  }else{
    dAC50temp = mergeAC50(dAC50, lchannel)  
  }
  
  model = generateSOM(ddesc, dAC50temp, 15, 15, paste(prout, "allActive",  sep = ""), model)
  
    
  # spe by column
  for (channel in lchannel){
    dir.create(paste(prout, channel, "/", sep = ""))
    dAC50temp = mergeAC50(dAC50, channel)
    #generateEnrichCurve(ddesc, dAC50temp , 50, paste(prout, AC50col, sep = ""))
    generateSOM(ddesc, dAC50temp, 15, 15, paste(prout, channel, "/",  sep = ""), model)
  }
}

# no apply for luciferase
if(lchannel[1] != "Luc_IC50"){
  # by color
  cbluehek293 = c("hek293_med_blue_n", "hek293_cell_blue_n")
  cgreenhek293 = c("hek293_med_green_n", "hek293_cell_green_n")
  credhek293 = c("hek293_med_red_n", "hek293_cell_red_n")
  
  cbluehepg2 = c("hepg2_med_blue_n", "hepg2_cell_blue_n")
  cgreenhepg2 = c("hepg2_med_green_n", "hepg2_cell_green_n")
  credhepg2 = c("hepg2_med_red_n", "hepg2_cell_red_n")
  
  cblue = c("hepg2_med_blue_n", "hepg2_cell_blue_n", "hek293_med_blue_n", "hek293_cell_blue_n")
  cgreen = c("hepg2_med_green_n", "hepg2_cell_green_n", "hek293_med_green_n", "hek293_cell_green_n")
  cred = c("hepg2_med_red_n", "hepg2_cell_red_n", "hek293_med_red_n", "hek293_cell_red_n")
  
  
  dAC50temp = mergeAC50(dAC50, lchannel)
  dir.create(paste(prout, "allcolors", "/", sep = ""))
  generateSOM(ddesc, dAC50temp, 15, 15, paste(prout, "allcolors", "/", sep = ""), model)
  
  dAC50temp = mergeAC50(dAC50, cbluehepg2)
  dir.create(paste(prout, "hepg2_blue", "/", sep = ""))
  generateSOM(ddesc, dAC50temp, 15, 15, paste(prout, "hepg2_blue", "/", sep = ""), model)
  
  dAC50temp = mergeAC50(dAC50, cgreenhepg2)
  dir.create(paste(prout, "hepg2_green", "/", sep = ""))
  generateSOM(ddesc, dAC50temp, 15, 15, paste(prout, "hepg2_green", "/", sep = ""), model)
  
  dAC50temp = mergeAC50(dAC50, credhepg2)
  dir.create(paste(prout, "hepg2_red", "/", sep = ""))
  generateSOM(ddesc, dAC50temp, 15, 15, paste(prout, "hepg2_red", "/", sep = ""), model)
  
  dAC50temp = mergeAC50(dAC50, cbluehek293)
  dir.create(paste(prout, "hek293_blue", "/", sep = ""))
  generateSOM(ddesc, dAC50temp, 15, 15, paste(prout, "hek293_blue", "/", sep = ""), model)
  
  dAC50temp = mergeAC50(dAC50, cgreenhek293)
  dir.create(paste(prout, "hek293_green", "/", sep = ""))
  generateSOM(ddesc, dAC50temp, 15, 15, paste(prout, "hek293_green", "/", sep = ""), model)
  
  dAC50temp = mergeAC50(dAC50, credhek293)
  dir.create(paste(prout, "hek293_red", "/", sep = ""))
  generateSOM(ddesc, dAC50temp, 15, 15, paste(prout, "hek293_red", "/", sep = ""), model)
  
  dAC50temp = mergeAC50(dAC50, cblue)
  dir.create(paste(prout, "blue", "/", sep = ""))
  generateSOM(ddesc, dAC50temp, 15, 15, paste(prout, "blue", "/", sep = ""), model)
  
  dAC50temp = mergeAC50(dAC50, cgreen)
  dir.create(paste(prout, "green", "/", sep = ""))
  generateSOM(ddesc, dAC50temp, 15, 15, paste(prout, "green", "/", sep = ""), model)
  
  dAC50temp = mergeAC50(dAC50, cred)
  dir.create(paste(prout, "red", "/", sep = ""))
  generateSOM(ddesc, dAC50temp, 15, 15, paste(prout, "red", "/", sep = ""), model)
  
}


