#!/usr/bin/env Rscript
library(ggplot2)


histActInact = function(din, prout){
  
  dact = rep (2, dim(din)[1])
  dact = as.factor(append(dact, rep(1, dim(din)[1])))
  
  ract = append(din[,2], din[,3])
  
  dplot = cbind(ract, dact)
  colnames(dplot) = c("Ratio", "Act")
  dplot = as.data.frame(dplot)
  
  dplot$Act = as.factor(dplot$Act)
  
  p <- ggplot(dplot, aes(x=Ratio, fill=Act)) + 
    labs(y = "Count of chemicals")+
    geom_histogram(binwidth=.05, alpha=.5, position="identity")
    #geom_boxplot() +
    #geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
    #geom_jitter(shape=16, position=position_jitter(0.2))
  
  ggsave(paste(prout, "_hist.png", sep = ""), dpi = 300, width = 8, height = 7)
  
}


histAct = function(din, prout, xlab){
  
  dact = rep (2, dim(din)[1])
  dact = as.factor(append(dact, rep(1, dim(din)[1])))
  
  ract = din[,2]
  
  dplot = cbind(ract, dact)
  colnames(dplot) = c("Ratio", "Act")
  dplot = as.data.frame(dplot)
  
  dplot$Act = as.factor(dplot$Act)
  
  p <- ggplot(dplot, aes(x=Ratio)) + 
    labs(y = "Count of chemicals", x=xlab)+
    geom_histogram(binwidth=.05, alpha=.5, position="identity", color = "grey", fill = "green")
  #geom_boxplot() +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
  #geom_jitter(shape=16, position=position_jitter(0.2))
  
  ggsave(paste(prout, "_Acthist.png", sep = ""), dpi = 300, width = 7, height = 8)
  
}


################
#     MAIN     #
################

args <- commandArgs(TRUE)
pratio = args[1]

#pratio = "~/interference/ToxCast_analysis/IC50Cluster/detection_technology_type/Fluorescence"



dratio = read.csv(pratio, sep = "\t", header = TRUE)
rownames(dratio) = dratio[,1]

#boxplot
histActInact(dratio[,c("CAS", "ActRatio", "InactRatio")], paste(pratio, "_ratio"))
histAct(dratio[,c("CAS", "ActRatio")], paste(pratio, "_ratioAct"), "Ratio of active assays")
histAct(dratio[,c("CAS", "ZAct")], paste(pratio, "_ratioZscore"), "Z-score ratio of active assays")
histActInact(dratio[,c("CAS", "ZAct", "Zinact")], paste(pratio, "_Zscore"))


