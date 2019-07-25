#!/usr/bin/env Rscript
source("~/development/Rglobal/source/dataManager.R")
source("~/development/Rglobal/source/PCAdrawer.R")
source("dendocircular.R")
source("clustering.R")
source("distributions.R")
source("cardMatrix.R")


PCAplot = function (din, daff, prout){
  
  dinScale = scale(din)
  
  data.cor=cor(dinScale)
  data.eigen=eigen(data.cor)
  lambda = data.eigen$values
  var_cap = lambda/sum(lambda)*100
  cp = data.eigen$vectors
  rownames (cp) = colnames (dinScale)
  colnames (cp) = colnames (dinScale)
  data_plot = as.matrix(dinScale)%*%cp
  factor = factorACP (data_plot, cp)
  data_plot = as.data.frame(data_plot[,c(1,2)])
  
  colnames(data_plot) = c("X", "Y")
  
  # projection point
  aff = daff[,1]
  p <- ggplot(data_plot, aes(X,Y)) + 
    geom_point(size=1.5, aes(color = aff)) + 
    scale_color_continuous(name='-log10(AC50)',low='red', high='lightgreen') +
    labs(x =  paste("CP1: ", signif (var_cap[1], 4), "%", sep = "") , y = paste("CP2: ", signif (var_cap[2], 4), "%", sep = ""))
  ggsave(paste(prout, "logAC50_PCA.png", sep = ""), dpi = 300, width = 8, height = 8)
  
  
  #aff = daff[,2]
  #p <- ggplot(data_plot, aes(X,Y)) + 
  #  geom_point(size=1.5, aes(color = aff)) + 
  #  scale_color_discrete(name='Tox class') +
  #  labs(x =  paste("CP1: ", signif (var_cap[1], 4), "%", sep = "") , y = paste("CP2: ", signif (var_cap[2], 4), "%", sep = ""))
  #ggsave(paste(prout, "class_PCA.png", sep = ""), dpi = 300, width = 8, height = 8)
  
  
  # projection descriptors #
  ##########################
  #col.desc = colorDesc(colnames(din))
  col.desc = "black"
  color_arrow = "black"
  print(col.desc)
  
  
  png (paste (prout, "PCA_descriptor.png", sep = ""), 1700, 1500)
  par(mar=c(8,8,8,8))
  plot(data_plot[,1],data_plot[,2], xlab = paste("CP1: ", signif (var_cap[1], 4), "%", sep = ""), ylab = paste("CP2: ", signif (var_cap[2], 4), "%", sep = ""), pch=20, cex.lab = 4, cex.main = 4, cex.axis = 1.75, cex = 4, type = "n")
  #points (data_plot[,1][length (color_point):dim(data_plot)[1]],data_plot[,2][length (color_point):dim(data_plot)[1]], pch=17, cex = 4, col = color_point2)
  abline(h=0,v=0)
  arrows (0,0,cp[,1]*factor,cp[,2]*factor, col = color_arrow, lwd = 4 )
  text (cp[,1]*factor, cp[,2]*factor,rownames(cp), col = color_arrow, cex = 2.5)
  dev.off()
  
  
  svg (file = paste (prout, "PCA_descriptor.svg", sep = ""), 25, 25)
  par(mar=c(8,8,8,8))
  plot(data_plot[,1],data_plot[,2], main = "", xlab = paste("CP1: ", signif (var_cap[1], 4), "%", sep = ""), ylab = paste("CP2: ", signif (var_cap[2], 4), "%", sep = ""), cex.lab = 6, cex.main = 4, cex.axis = 1.75, cex = 6, type = "n")
  #points (data_plot[,1][length (color_point):dim(data_plot)[1]],data_plot[,2][length (color_point):dim(data_plot)[1]], pch=17, cex = 4, col = color_point2)
  abline(h=0,v=0)
  arrows (0,0,cp[,1]*factor,cp[,2]*factor, col = color_arrow, lwd = 3 )
  text (cp[,1]*factor, cp[,2]*factor,rownames(cp), col = color_arrow, cex = 3.5)
  dev.off()
  
}




################
#     MAIN     #
################

args <- commandArgs(TRUE)
pdesc = args[1]
pdata = args[2] #to take affinity
prout = args[3]
valcor = as.double(args[4])
maxQuantile = as.double(args[5])
nbNA = as.integer(args[6])
#logaff = as.integer(args[6])
#plotPCA = as.integer(args[7])
#corMatrix = as.integer(args[8])
#histplot = as.integer(args[9])
#circularDendo = as.integer(args[10])
#optimal_clustering = as.integer(args[11])


#pdesc = "/home/borrela2/interference/luc-biochem/test/desc-test.csv"
#pdata = "/home/borrela2/interference/luc-biochem/test/LD50.csv"
#prout = "/home/borrela2/interference/luc-biochem/test/"

plotPCA = 1
corMatrix = 1
histplot = 1
circularDendo = 1
#valcor = 0.80
#maxQuantile = 85
logaff = 1
optimal_clustering = 0


# Process descriptors matrix #
##############################
dglobal = openDataVexcluded(pdesc, valcor, prout, c(1,2), nbNA)
dglobal = dglobal[[1]]


rownames(dglobal) = dglobal[,1]
dglobal = dglobal[,-1]
dglobal = dglobal[,-1]

# remove descriptor with a distribution on one quantile
dglobal = delnohomogeniousdistribution(dglobal, maxQuantile)

#######################
# order with affinity #
#######################
# Opening
daffinity = read.csv(pdata, sep = "\t", header = TRUE)
# fusion 
daffinity = mergeIdenticRow(daffinity)

rownames(daffinity) = daffinity[,1]
daffinity = daffinity[,-1]
# remove class + remove NA
#daffinity = daffinity[,-2]

print(dim(daffinity))



if(logaff == 1){
  histAff(daffinity, prout)
  daffinity[,1] = -log10(daffinity[,1])
  histlogAff(daffinity, prout)
}

# remove no selected affinity or bad quality from manual curation
print(dim(dglobal))
lID = intersect(rownames(daffinity), rownames(dglobal))

dglobal = dglobal[lID,]
daffinity = daffinity[lID,]

write.csv(dglobal, paste(prout, "globalTable.csv", sep = ""))

#ord = NULL
#for(i in seq(1,dim(dglobal)[1])){
#  ipos = which(daffinity[,1] == rownames(dglobal)[i])
#  ord = append (ord,ipos)
#}

#daffinity = daffinity[ord,]


#orderaff = order(daffinity[,2],decreasing=T)
#daffinity = daffinity[orderaff,]


#dglobal = dglobal[orderaff,]
#print(rownames(dglobal))
#print(colnames(dglobal))

#print(dim(dglobal))
#print(dim(daffinity))

#print(rownames(dglobal))
#print(rownames(daffinity))



if (corMatrix == 1){
  cardMatrixCor(cor(dglobal), paste(prout, "matrixCor_", valcor, sep = ""), 6)
  
  # matrix cor with pMIC
  lcor = NULL
  ldesc = colnames(dglobal)
  for(desc in ldesc){
    lvaldesc = dglobal[,desc]
    corval = cor(daffinity[,1], lvaldesc, use = "pairwise.complete.obs")
    lcor = append(lcor, corval)
  }
  names(lcor) = ldesc
  
  write.csv(lcor, file = paste(prout, "corAC50-Desc.csv"))
}


if(plotPCA == 1){
  PCAplot(dglobal, daffinity, prout)
}


if (histplot == 1){
  histDataOne(data1 = dglobal, paste(prout, "histDesc_", valcor, ".pdf", sep = ""))
}

if (circularDendo == 1){
  dendogramCircle(dglobal, daffinity, prout)
}


