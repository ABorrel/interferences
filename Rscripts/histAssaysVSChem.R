#!/usr/bin/env Rscript
library(ggplot2)

histplot = function(din, xlab, ylab, pout){
  
  din = as.data.frame(din)
  
  p <- ggplot(din, aes(x=test)) + 
    labs(y = ylab, x=xlab)+
    geom_histogram(binwidth=.05, alpha=1, position="identity", color = "grey", fill = "green")
  #geom_boxplot() +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
  #geom_jitter(shape=16, position=position_jitter(0.2))
  
  print(p)
  
  ggsave(paste(pout, ".png", sep = ""), dpi = 300, width = 7, height = 8)
  
  
}


################
#     MAIN     #
################

args <- commandArgs(TRUE)
pAssays = args[1]
pchem = args[2]
prout = args[3]

#pAssays = "/home/borrela2/interference/ToxCast_analysis/prelimAnalysis/histAssays/CountChembyAssays"
#pchem = "/home/borrela2/interference/ToxCast_analysis/prelimAnalysis/histAssays/CountAssaysbyChem"
#prout = "/home/borrela2/interference/ToxCast_analysis/prelimAnalysis/histAssays/"


dassays = read.csv(pAssays, sep = "\t", header = TRUE)
rownames(dassays) = dassays[,1]
dchem = read.csv(pchem, sep = "\t", header = TRUE)

print (dim(dassays))

histplot(dassays, "Number of chemical tested by assays","Count", paste(prout, "histAssays"))
histplot(dchem, "Number of assays tested by chemicals","Count", paste(prout, "histChemicals"))