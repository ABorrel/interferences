#!/usr/bin/env Rscript
library(ggplot2)
library(plyr)


# luciferase #
##############

multipleHistluc = function(d, prout){

  ggplot(d, aes(x=-log10(set1))) + 
      geom_histogram(aes(y=..density..),     # Histogram with density instead of count on y-axis
                     binwidth=.5,
                     colour="black", fill="white") +
    labs(x = "-log10(AC50)", y = "Frequencies") + 
  xlim (c(-4, 4))+
    ylim(c(0, 0.90))+
    theme(axis.text.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 25, hjust = 0.5, vjust =0.1))+
      geom_density(alpha=.2, fill="#FF6666")
  ggsave(paste(prout, "Set1.png", sep = ""), dpi = 300, width = 8, height = 7)
  
  
  ggplot(d, aes(x=-log10(set2))) + 
    geom_histogram(aes(y=..density..),     # Histogram with density instead of count on y-axis
                   binwidth=.5,
                   colour="black", fill="white") +
    labs(x = "-log10(AC50)", y = "Frequencies") + 
    xlim (c(-4, 4))+
    ylim(c(0, 0.90))+
    theme(axis.text.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 25, hjust = 0.5, vjust =0.1))+
    geom_density(alpha=.2, fill="#007F00")
  ggsave(paste(prout, "Set2.png", sep = ""), dpi = 300, width = 8, height = 7)
  
  
  ggplot(d, aes(x=-log10(set3))) + 
    geom_histogram(aes(y=..density..),     # Histogram with density instead of count on y-axis
                   binwidth=.5,
                   colour="black", fill="white") +
    labs(x = "-log10(AC50)", y = "Frequencies") + 
    xlim (c(-4, 4))+
    ylim(c(0, 0.90))+
    theme(axis.text.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 25, hjust = 0.5, vjust =0.1))+
    geom_density(alpha=.2, fill="#0e12df")
  ggsave(paste(prout, "Set3.png", sep = ""), dpi = 300, width = 8, height = 7)
  
  
}


combineluchist = function(d, prout){
  
  # combine plot
  dset1 = cbind(d$set1, rep("set1", length(d$set1)))
  dset2 = cbind(d$set2, rep("set2", length(d$set2)))
  dset3 = cbind(d$set3, rep("set3", length(d$set3)))
  
  din = rbind(dset1, dset2)
  din = rbind(din, dset3)
  din = na.omit(din)
  
  colnames(din) = c("AC50", "set")
  din[,1] = -log10(as.double(din[,1]))
  din = transform(din, AC50=as.numeric(as.character(AC50)))
  
  
  #din = data.frame(din)
  mu <- ddply(din, "set", summarise, grp.mean=mean(AC50))
  print(mu)
  
  ggplot(din, aes(x=AC50, color=set, fill=set)) +
    geom_histogram(aes(y=..density..), position="identity", alpha=0.5)+
    geom_density(alpha=0.6)+
    geom_vline(data=mu, aes(xintercept=grp.mean, color=set),
               linetype="dashed")+
    #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
    #scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
    labs(title="Histogram plot",x="-log10(AC50)", y = "Density")+
    theme_classic()
  
  ggsave(paste(prout, "combine.png", sep = ""), dpi = 300, width = 8, height = 7)
  
}




# hepg2 #
#########


multipleHisthepg2 = function(d, prout){
  
  # med red #
  ggplot(d, aes(x=-log10(med_red))) + 
    geom_histogram(aes(y=..density..),     # Histogram with density instead of count on y-axis
                   binwidth=.5,
                   colour="black", fill="white") +
    labs(x = "-log10(AC50)", y = "Frequencies") + 
    xlim (c(-4, 4))+
    ylim(c(0, 1.0))+
    theme(axis.text.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 25, hjust = 0.5, vjust =0.1))+
    geom_density(alpha=.2, fill="#c70000")
  ggsave(paste(prout, "med_red.png", sep = ""), dpi = 300, width = 8, height = 7)
  
  # med blue #
  ggplot(d, aes(x=-log10(med_blue))) + 
    geom_histogram(aes(y=..density..),     # Histogram with density instead of count on y-axis
                   binwidth=.5,
                   colour="black", fill="white") +
    labs(x = "-log10(AC50)", y = "Frequencies") + 
    xlim (c(-4, 4))+
    ylim(c(0, 1.0))+
    theme(axis.text.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 25, hjust = 0.5, vjust =0.1))+
    geom_density(alpha=.2, fill="#0003c7")
  ggsave(paste(prout, "med_blue.png", sep = ""), dpi = 300, width = 8, height = 7)
  
  # med_green
  ggplot(d, aes(x=-log10(med_green))) + 
    geom_histogram(aes(y=..density..),     # Histogram with density instead of count on y-axis
                   binwidth=.5,
                   colour="black", fill="white") +
    labs(x = "-log10(AC50)", y = "Frequencies") + 
    xlim (c(-4, 4))+
    ylim(c(0, 0.90))+
    theme(axis.text.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 25, hjust = 0.5, vjust =0.1))+
    geom_density(alpha=.2, fill="#0e12df")
  ggsave(paste(prout, "med_green.png", sep = ""), dpi = 300, width = 8, height = 7)
  
  
  
  
  # med_blue_n
  ggplot(d, aes(x=-log10(med_blue_n))) + 
    geom_histogram(aes(y=..density..),     # Histogram with density instead of count on y-axis
                   binwidth=.5,
                   colour="black", fill="white") +
    labs(x = "-log10(AC50)", y = "Frequencies") + 
    xlim (c(-4, 4))+
    ylim(c(0, 0.90))+
    theme(axis.text.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 25, hjust = 0.5, vjust =0.1))+
    geom_density(alpha=.2, fill="#0e12df")
  ggsave(paste(prout, "med_blue_n.png", sep = ""), dpi = 300, width = 8, height = 7)
  
  
  
  # cell_green
  ggplot(d, aes(x=-log10(cell_green))) + 
    geom_histogram(aes(y=..density..),     # Histogram with density instead of count on y-axis
                   binwidth=.5,
                   colour="black", fill="white") +
    labs(x = "-log10(AC50)", y = "Frequencies") + 
    xlim (c(-4, 4))+
    ylim(c(0, 1.0))+
    theme(axis.text.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 25, hjust = 0.5, vjust =0.1))+
    geom_density(alpha=.2, fill="#7af4c3")
  ggsave(paste(prout, "cell_green.png", sep = ""), dpi = 300, width = 8, height = 7)
  
  
  
  # cell_blue
  ggplot(d, aes(x=-log10(cell_blue))) + 
    geom_histogram(aes(y=..density..),     # Histogram with density instead of count on y-axis
                   binwidth=.5,
                   colour="black", fill="white") +
    labs(x = "-log10(AC50)", y = "Frequencies") + 
    xlim (c(-4, 4))+
    ylim(c(0, 0.90))+
    theme(axis.text.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 25, hjust = 0.5, vjust =0.1))+
    geom_density(alpha=.2, fill="#0e12df")
  ggsave(paste(prout, "cell_blue.png", sep = ""), dpi = 300, width = 8, height = 7)
  
  
  
  
  # cell_red
  ggplot(d, aes(x=-log10(cell_red))) + 
    geom_histogram(aes(y=..density..),     # Histogram with density instead of count on y-axis
                   binwidth=.5,
                   colour="black", fill="white") +
    labs(x = "-log10(AC50)", y = "Frequencies") + 
    xlim (c(-4, 4))+
    ylim(c(0, 0.90))+
    theme(axis.text.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 25, hjust = 0.5, vjust =0.1))+
    geom_density(alpha=.2, fill="#e3230d")
  ggsave(paste(prout, "cell_red.png", sep = ""), dpi = 300, width = 8, height = 7)
  
  
  
  
  # cell_blue_n
  ggplot(d, aes(x=-log10(cell_blue_n))) + 
    geom_histogram(aes(y=..density..),     # Histogram with density instead of count on y-axis
                   binwidth=.5,
                   colour="black", fill="white") +
    labs(x = "-log10(AC50)", y = "Frequencies") + 
    xlim (c(-4, 4))+
    ylim(c(0, 0.90))+
    theme(axis.text.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 25, hjust = 0.5, vjust =0.1))+
    geom_density(alpha=.2, fill="#0e12df")
  ggsave(paste(prout, "cell_blue_n.png", sep = ""), dpi = 300, width = 8, height = 7)
  
  
  
  
}


combinehepg2hist = function(d, prout){
  
  # combine plot
  dmed_red = cbind(d$med_red, rep("med_red", length(d$med_red)))
  dmed_blue = cbind(d$med_blue, rep("med_blue", length(d$med_blue)))
  dmed_green = cbind(d$med_green, rep("med_green", length(d$med_green)))
  dcell_green = cbind(d$cell_green, rep("cell_green", length(d$cell_green)))
  dcell_blue = cbind(d$cell_blue, rep("cell_blue", length(d$cell_blue)))
  dcell_red = cbind(d$cell_red, rep("cell_red", length(d$cell_red)))
  dcell_blue_n = cbind(d$cell_blue_n, rep("cell_blue_n", length(d$cell_blue_n)))
  dmed_blue_n = cbind(d$med_blue_n, rep("med_blue_n", length(d$med_blue_n)))
  
  
  din = rbind(dmed_red, dmed_blue)
  din = rbind(din, dmed_green)
  din = rbind(din, dcell_green)
  din = rbind(din, dcell_blue)
  din = rbind(din, dcell_red)
  din = rbind(din, dcell_blue_n)
  din = rbind(din, dmed_blue_n)
  
  din = na.omit(din)
  
  colnames(din) = c("AC50", "Fluo")
  din[,1] = -log10(as.double(din[,1]))
  din = transform(din, AC50=as.numeric(as.character(AC50)))
  
  
  #din = data.frame(din)
  mu <- ddply(din, "Fluo", summarise, grp.mean=mean(AC50))
  print(mu)
  
  ggplot(din, aes(x=AC50, color=Fluo, fill=Fluo)) +
    geom_histogram(aes(y=..density..), position="identity", alpha=0.5)+
    geom_density(alpha=0.6)+
    geom_vline(data=mu, aes(xintercept=grp.mean, color=Fluo),
               linetype="dashed")+
    #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
    #scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
    labs(title="Histogram plot",x="-log10(AC50)", y = "Density")+
    theme_classic()
  
  ggsave(paste(prout, "combine.png", sep = ""), dpi = 300, width = 8, height = 7)
  
}





################
#     MAIN     #
################

args <- commandArgs(TRUE)
pAC50 = args[1]
prout = args[2]
nameassays = args[3]


#pAC50 = "/home/borrela2/interference/spDataAnalysis/tox21-spec-hepg2-p1/AC50_sample"
#prout = "/home/borrela2/interference/spDataAnalysis/tox21-spec-hepg2-p1/histIC50/"
#nameassays = "hepg2"

dAC50 = read.csv(pAC50, sep = "\t", header = TRUE)

if(nameassays == "luc"){
  multipleHistluc(dAC50, prout)
  combineluchist(dAC50, prout)
}else{
  multipleHisthepg2(dAC50, prout)
  combinehepg2hist(dAC50, prout)
}


