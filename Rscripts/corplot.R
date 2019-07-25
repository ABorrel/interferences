#!/usr/bin/env Rscript
library(ggplot2)

################
#     MAIN     #
################

args <- commandArgs(TRUE)
pfilin = args[1]

#pfilin = "/home/borrela2/interference/spDataAnalysis/cor/tox21-spec-hek293-p1_tox21-spec-hepg2-p1"

# open affinity file #
######################
din = read.csv(pfilin, sep = "\t", header = TRUE)
rownames(din) = din[,1]
din = na.omit(din)

corval = cor(din[,2],din[,3])

p = ggplot(din, aes(AC50_1, AC50_2))+
  geom_point(size=1.5, col="black", shape=21) + 
  geom_text(x=10, y=95, label = paste("r=",round(corval,2), sep = ""), size = 8)+
  theme(axis.text.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 25, hjust = 0.5, vjust =0.1))+
  labs(x = expression("AC50"), y =expression("AC50"))
  #xlim (c(-2.5, 2.5)) +
  #geom_segment(aes(x = -2.5, y = -2.5, xend = 2.5, yend = 2.5)) + 
  #ylim (c(-2.5, 2.5)) 
#print(p)
ggsave(paste(pfilin, ".png",sep=""), width = 6,height = 6, dpi = 300)






