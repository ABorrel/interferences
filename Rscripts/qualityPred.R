#!/usr/bin/env Rscript
library(ggplot2)
library(plyr)


################
#     MAIN     #
################

args <- commandArgs(TRUE)
pprob = args[1]
prout = args[2]

#pprob = "/home/borrela2/interference/testing/587_TexasRed/HepG2/cell_red_n/LDAclass/sumProb"
#prout = "/home/borrela2/interference/testing/587_TexasRed/HepG2/cell_red_n/LDAclass/"

#pprob = "/home/borrela2/interference/spDataAnalysis/tox21-luc-biochem-p1/QSARclass/Prob/IC50_train"
#prout = "/home/borrela2/interference/spDataAnalysis/tox21-luc-biochem-p1/QSARclass/Prob/"

#pprob = "./../../trash/IC50_train"
#prout = "./../../"



d = read.csv(pprob, sep = "\t")
if (dim(d)[2] == 1){
  d = read.csv(pprob, sep = ",")  
}
rownames(d) = d[,1]
d = d[,-1]

d[which(d[,"Real"] == 0),"Real"] = NA
dact = na.omit(d)
dinact = d[which(is.na(d[,3])),]


# Use geom_pointrange
ggplot(dact, aes(x=Real, y=Mpred)) + 
  geom_pointrange(aes(ymin=Mpred-SDpred, ymax=Mpred+SDpred))

ggsave(paste(pprob, ".png", sep = ""),  width = 8, height = 8, dpi = 300, bg="transparent")



lact = rep("act", dim(d)[1])
lact[which(is.na(d[,3]))] = "inact"

dhist = d
dhist[,3] = lact

mu <- ddply(dhist, "Real", summarise, grp.mean=mean(Mpred))


ggplot(dhist, aes(x=Mpred, color=Real, fill=Real)) +
    #geom_histogram(aes(y=..density..), position="identity", alpha=0.15)+
    geom_density(alpha=0.3)+
    geom_vline(data=mu, aes(xintercept=grp.mean, color=Real),
               linetype=c("dashed", "solid"))+
    #xlim(0,120)+
    theme(text = element_text(size=19))+
    scale_fill_manual(values=c("#eb9999", "#290000"), labels = c("active", "inactive"))+
    #scale_color_manual(values=c("#eb9999", "#290000", "#d21919", "#a40000"), labels = c("hek293 cell based", "hek2$
    labs(title="",x="Prob", y = "Density")

  ggsave(paste(prout, "hist.png", sep = ""),  width = 8, height = 7, dpi = 300, bg="transparent")




h = c("NB active", "NB inact", "TP", "FN", "TN", "FP", "Mact prob", "SDact prob", "Minact prob", "SDinact prob")

nbact = dim(dact)[1]
nbinact = 
TP = length(which(dact[,1] >= 0.5))
FN = length(which(dact[,1] < 0.5))
Mact = mean(dact[,1])
SDact = sd(dact[,1])

TN = length(which(dinact[,1] < 0.5))
FP = length(which(dinact[,1] >= 0.5))
Minact = mean(dinact[,1])
SDinact = sd(dinact[,1])

dsum = rbind(h, c(nbact, nbinact, TP, FN, TN, FP, Mact,SDact, Minact, SDinact))
write.csv(dsum, file = paste(pprob, ".csv", sep = ""))


