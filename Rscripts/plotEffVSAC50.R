#!/usr/bin/env Rscript
library(ggplot2)
source("~/development/Rglobal/source/dataManager.R")


################
#     MAIN     #
################

args <- commandArgs(TRUE)
pAC50Eff = args[1]


#pAC50Eff = "/home/borrela2/interference/spDataAnalysis/tox21-luc-biochem-p1/efficient/eff"
d = read.csv(pAC50Eff, sep = "\t", header = TRUE)
rownames(d) = d[,1]
d = d[,-1]

i = 1
imax = dim(d)[2]
while(i<imax){
  dtemp = cbind(d[,i], d[,i+1])
  dtemp = na.omit(dtemp)
  delval = which(dtemp[,2] >= 150)
  if (is.integer0(delval) != TRUE){
    dtemp = dtemp[-delval,]
  }
  
  # correlation 
  #png(paste(pAC50Eff, "_", colnames(d)[i], sep =""), 500, 500)
  #plot(log10(dtemp[,1]), dtemp[,2], pch = 20, main = paste("Correlation = ", round(cor(dtemp[,1], dtemp[,2]), digits = 3)), cex = 2, xlab = "log10(AC50)", ylab = "Efficiency")
  #abline(a = 0, b = 1, col = "red", cex = 3)
  #dev.off()

  colnames(dtemp) = c("AC50", "Efficiency")
  dtemmp = as.data.frame(dtemp)
  dtemp = transform(dtemp, Efficiency=as.numeric(as.character(Efficiency)))
  # histogram
  print(head(dtemp))  
  
  ggplot(dtemp, aes(x=Efficiency)) +
    geom_histogram(binwidth=5, alpha=1, position="identity", color = "grey", fill = "green")
    #geom_histogram(aes(y=..density..), position="identity", alpha=0.15)
    #geom_density(alpha=10)
    #theme(text = element_text(size=19))+
    #scale_fill_manual(values=c("#eb9999", "#290000", "#d21919", "#a40000"), labels = c("hek293 cell based", "hek293 cell free", "hepg2 cell based", "hepg2 cell free"))+
    #scale_color_manual(values=c("#eb9999", "#290000", "#d21919", "#a40000"), labels = c("hek293 cell based", "hek293 cell free", "hepg2 cell based", "hepg2 cell free")) +
    #labs(title="",x="log(AC50) (uM)", y = "Density")
  ggsave(paste(pAC50Eff, "_hist_", colnames(d)[i],".png", sep =""), width = 8, height = 7, dpi = 300, bg = "transparent")


 
  i = i + 2
}

