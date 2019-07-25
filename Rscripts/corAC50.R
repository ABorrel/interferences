#!/usr/bin/env Rscript
library(ggplot2)

################
#     MAIN     #
################

args <- commandArgs(TRUE)
pAC50 = args[1]
pCurve = args[2]#classes of interferences 
prout = args[3]

#pAC50 = "/home/borrela2/interference/spDataAnalysis/tox21-spec-hepg2-p1/AC50_sample"
#pCurve = "/home/borrela2/interference/spDataAnalysis/tox21-spec-hepg2-p1/corAC50/curve"
#prout = "/home/borrela2/interference/spDataAnalysis/tox21-spec-hepg2-p1/corAC50/"

dAC50 = read.csv(pAC50, header = TRUE, sep = "\t")
dcurve = read.csv(pCurve, sep = "\t", header = TRUE)

print(head(dAC50))

i = 2
imax = dim(dAC50)[2]-1
print(imax)
while(i <= imax-1){
  j = i + 1
  while(j <= imax){
    mattemp = cbind(dAC50[i], dAC50[j])
    mattemp = na.omit(mattemp)
    mattemp = -log10(mattemp)
    corval = cor(mattemp[1], mattemp[2])
    colnames(mattemp) = c("AC50_1","AC50_2")
    
    p = ggplot(mattemp, aes(AC50_1, AC50_2))+
      geom_point(size=1.5, col="black", shape=21) + 
      geom_text(x=-2, y=0, label = paste("r=",round(corval,2), sep = ""), size = 8)+
      theme(axis.text.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 25, hjust = 0.5, vjust =0.1))+
      labs(x = expression("-log AC50"), y =expression("-log AC50"))
    #xlim (c(-2.5, 2.5)) +
    #geom_segment(aes(x = -2.5, y = -2.5, xend = 2.5, yend = 2.5)) + 
    #ylim (c(-2.5, 2.5)) 
    #print(p)
    ggsave(paste(prout, colnames(dAC50[i]), "-",colnames(dAC50[j]), ".png", sep=""), width = 10,height = 10, dpi = 300)
    j = j + 1
    
  }
  i = i+ 1
}

