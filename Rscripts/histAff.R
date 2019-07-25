#!/usr/bin/env Rscript
library(ggplot2)
library(plyr)

combinered = function(d, prout){
  
  # combine plot
  dmed_hepg2_cell_red = cbind(d$hepg2_cell_red_n, rep("hepg2_cell_red_n", length(d$hepg2_cell_red_n)))
  dmed_hepg2_med_red = cbind(d$hepg2_med_red_n, rep("hepg2_med_red_n", length(d$hepg2_med_red_n)))
  dmed_hek293_cell_red = cbind(d$hek293_cell_red_n, rep("hek293_cell_red_n", length(d$hek293_cell_red_n)))
  dmed_hek293_med_red = cbind(d$hek293_med_red_n, rep("hek293_med_red_n", length(d$hek293_med_red_n)))
  
  din = rbind(dmed_hepg2_cell_red, dmed_hepg2_med_red)
  din = rbind(din, dmed_hek293_cell_red)
  din = rbind(din, dmed_hek293_med_red)
  din = na.omit(din)
  
  colnames(din) = c("AC50", "Assays")
  #din[,1] = log10(as.double(din[,1]))
  din = transform(din, AC50=as.numeric(as.character(AC50)))
  
  
  #din = data.frame(din)
  mu <- ddply(din, "Assays", summarise, grp.mean=mean(AC50))
  #print(mu)
  
  ggplot(din, aes(x=AC50, color=Assays, fill=Assays)) +
    #geom_histogram(aes(y=..density..), position="identity", alpha=0.15)+
    geom_density(alpha=0.3)+
    geom_vline(data=mu, aes(xintercept=grp.mean, color=Assays),
               linetype=c("dashed", "solid", "dashed", "solid"))+
    xlim(0,120)+
    theme(text = element_text(size=19))+
    scale_fill_manual(values=c("#eb9999", "#290000", "#d21919", "#a40000"), labels = c("Hek293 cell based", "Hek293 cell free", "HepG2 cell based", "HepG2 cell free"))+
    scale_color_manual(values=c("#eb9999", "#290000", "#d21919", "#a40000"), labels = c("Hek293 cell based", "Hek293 cell free", "HepG2 cell based", "HepG2 cell free")) + 
    labs(title="",x="AC50 (uM)", y = "Density")
  
  ggsave(paste(prout, "red.png", sep = ""),  width = 8, height = 7, dpi = 300, bg="transparent")
  
  
  din$AC50 = log10(din$AC50)
  mu <- ddply(din, "Assays", summarise, grp.mean=mean(AC50))
  ggplot(din, aes(x=AC50, color=Assays, fill=Assays)) +
    #geom_histogram(aes(y=..density..), position="identity", alpha=0.15)+
    geom_density(alpha=0.3)+
    geom_vline(data=mu, aes(xintercept=grp.mean, color=Assays),
               linetype=c("dashed", "solid", "dashed", "solid"))+
    xlim(-1.5,3)+
    theme(text = element_text(size=19))+
    scale_fill_manual(values=c("#eb9999", "#290000", "#d21919", "#a40000"), labels = c("Hek293 cell based", "Hek293 cell free", "HepG2 cell based", "HepG2 cell free"))+
    scale_color_manual(values=c("#eb9999", "#290000", "#d21919", "#a40000"), labels = c("Hek293 cell based", "Hek293 cell free", "HepG2 cell based", "HepG2 cell free")) + 
    labs(title="",x="log(AC50) (uM)", y = "Density")
  
  ggsave(paste(prout, "logred.png", sep = ""),  width = 8, height = 7, dpi = 300, bg="transparent")
  
  
}




combineblue = function(d, prout){
  
  # combine plot
  dmed_hepg2_cell_blue = cbind(d$hepg2_cell_blue_n, rep("hepg2_cell_blue_n", length(d$hepg2_cell_blue_n)))
  dmed_hepg2_med_blue = cbind(d$hepg2_med_blue_n, rep("hepg2_med_blue_n", length(d$hepg2_med_blue_n)))
  dmed_hek293_cell_blue = cbind(d$hek293_cell_blue_n, rep("hek293_cell_blue_n", length(d$hek293_cell_blue_n)))
  dmed_hek293_med_blue = cbind(d$hek293_med_blue_n, rep("hek293_med_blue_n", length(d$hek293_med_blue_n)))
  
  din = rbind(dmed_hepg2_cell_blue, dmed_hepg2_med_blue)
  din = rbind(din, dmed_hek293_cell_blue)
  din = rbind(din, dmed_hek293_med_blue)
  din = na.omit(din)
  
  colnames(din) = c("AC50", "Assays")
  #din[,1] = log10(as.double(din[,1]))
  din = transform(din, AC50=as.numeric(as.character(AC50)))
  
  
  #din = data.frame(din)
  mu <- ddply(din, "Assays", summarise, grp.mean=mean(AC50))
  #print(mu)
  
  ggplot(din, aes(x=AC50, color=Assays, fill=Assays)) +
    #geom_histogram(aes(y=..density..), position="identity", alpha=0.15)+
    geom_density(alpha=0.3)+
    geom_vline(data=mu, aes(xintercept=grp.mean, color=Assays),
               linetype=c("dashed", "solid", "dashed", "solid"))+
    xlim(0,120)+
    theme(text = element_text(size=19))+
    scale_fill_manual(values=c("#cde2ff", "#032c66", "#6aa8ff", "#044299"), labels = c("Hek293 cell based", "Hek293 cell free", "HepG2 cell based", "HepG2 cell free"))+
    scale_color_manual(values=c("#cde2ff", "#032c66", "#6aa8ff", "#044299"), labels = c("Hek293 cell based", "Hek293 cell free", "HepG2 cell based", "HepG2 cell free")) + 
    labs(title="",x="AC50 (uM)", y = "Density")
  
  ggsave(paste(prout, "blue.png", sep = ""),  width = 8, height = 7, dpi = 300, bg="transparent")
  
  
  din$AC50 = log10(din$AC50)
  mu <- ddply(din, "Assays", summarise, grp.mean=mean(AC50))
  ggplot(din, aes(x=AC50, color=Assays, fill=Assays)) +
    #geom_histogram(aes(y=..density..), position="identity", alpha=0.15)+
    geom_density(alpha=0.3)+
    geom_vline(data=mu, aes(xintercept=grp.mean, color=Assays),
               linetype=c("dashed", "solid", "dashed", "solid"))+
    xlim(-1.5,3)+
    theme(text = element_text(size=19))+
    scale_fill_manual(values=c("#cde2ff", "#032c66", "#6aa8ff", "#044299"), labels = c("Hek293 cell based", "Hek293 cell free", "HepG2 cell based", "HepG2 cell free"))+
    scale_color_manual(values=c("#cde2ff", "#032c66", "#6aa8ff", "#044299"), labels = c("Hek293 cell based", "Hek293 cell free", "HepG2 cell based", "HepG2 cell free")) + 
    labs(title="",x="log(AC50) (uM)", y = "Density")
  
  ggsave(paste(prout, "logblue.png", sep = ""),  width = 8, height = 7, dpi = 300, bg="transparent")
  
}


combinegreen = function(d, prout){
  
  # combine plot
  dmed_hepg2_cell_green = cbind(d$hepg2_cell_green_n, rep("hepg2_cell_green_n", length(d$hepg2_cell_green_n)))
  dmed_hepg2_med_green = cbind(d$hepg2_med_green_n, rep("hepg2_med_green_n", length(d$hepg2_med_green_n)))
  dmed_hek293_cell_green = cbind(d$hek293_cell_green_n, rep("hek293_cell_green_n", length(d$hek293_cell_green_n)))
  dmed_hek293_med_green = cbind(d$hek293_med_green_n, rep("hek293_med_green_n", length(d$hek293_med_green_n)))
  
  din = rbind(dmed_hepg2_cell_green, dmed_hepg2_med_green)
  din = rbind(din, dmed_hek293_cell_green)
  din = rbind(din, dmed_hek293_med_green)
  din = na.omit(din)
  
  colnames(din) = c("AC50", "Assays")
  #din[,1] = log10(as.double(din[,1]))
  din = transform(din, AC50=as.numeric(as.character(AC50)))
  
  
  #din = data.frame(din)
  mu <- ddply(din, "Assays", summarise, grp.mean=mean(AC50))
  #print(mu)
  
  ggplot(din, aes(x=AC50, color=Assays, fill=Assays)) +
    #geom_histogram(aes(y=..density..), position="identity", alpha=0.15)+
    geom_density(alpha=0.3)+
    geom_vline(data=mu, aes(xintercept=grp.mean, color=Assays),
               linetype=c("dashed", "solid", "dashed", "solid"))+
    xlim(0,120)+
    theme(text = element_text(size=19))+
    scale_fill_manual(values=c("#c4f0b2", "#0b2900", "#62d732", "#298f00"), labels = c("Hek293 cell based", "Hek293 cell free", "HepG2 cell based", "Hepg2 cell free"))+
    scale_color_manual(values=c("#c4f0b2", "#0b2900", "#62d732", "#298f00"), labels = c("Hek293 cell based", "Hek293 cell free", "HepG2 cell based", "HepG2 cell free")) + 
    labs(title="",x="AC50 (uM)", y = "Density")
  
  ggsave(paste(prout, "green.png", sep = ""),  width = 8, height = 7, dpi = 300, bg="transparent")
  
  
  din$AC50 = log10(din$AC50)
  mu <- ddply(din, "Assays", summarise, grp.mean=mean(AC50))
  ggplot(din, aes(x=AC50, color=Assays, fill=Assays)) +
    #geom_histogram(aes(y=..density..), position="identity", alpha=0.15)+
    geom_density(alpha=0.3)+
    geom_vline(data=mu, aes(xintercept=grp.mean, color=Assays),
               linetype=c("dashed", "solid", "dashed", "solid"))+
    xlim(-1.5,3)+
    theme(text = element_text(size=19))+
    scale_fill_manual(values=c("#c4f0b2", "#0b2900", "#62d732", "#298f00"), labels = c("Hek293 cell based", "Hek293 cell free", "HepG2 cell based", "HepG2 cell free"))+
    scale_color_manual(values=c("#c4f0b2", "#0b2900", "#62d732", "#298f00"), labels = c("Hek293 cell based", "Hek293 cell free", "HepG2 cell based", "HepG2 cell free")) + 
    labs(title="",x="log(AC50) (uM)", y = "Density")
  
  ggsave(paste(prout, "loggreen.png", sep = ""),  width = 8, height = 7, dpi = 300, bg="transparent")
  
  
}


combineluc = function(d, prout){
  
  # combine plot
  dluc =  cbind(d$Luc_IC50, rep("Luc_IC50", length(d$Luc_IC50)))
  #dluc = dluc[-which(d[,1]>200),]
  
  din = dluc
  colnames(din) = c("IC50", "Assays")
  #din[,1] = log10(as.double(din[,1]))
  din = transform(din, IC50=as.numeric(as.character(IC50)))
  din = na.omit(din)
  
  din = data.frame(din)
  mu <- ddply(din, "Assays", summarise, grp.mean=mean(IC50))
  
  print(din)
  
  ggplot(din, aes(x=IC50, color=Assays, fill=Assays)) +
    #geom_histogram(aes(y=..density..), position="identity", alpha=0.15)+
    geom_density(alpha=0.3)+
    xlim(0,120)+
    theme(text = element_text(size=19))+
    geom_vline(data=mu, aes(xintercept=grp.mean, color=Assays),
               linetype=c("dashed"))+
    scale_fill_manual(values=c("#585858"), labels = c("IC50"))+
    scale_color_manual(values=c("#585858"), labels = c("IC50")) + 
    labs(title="",x="IC50 (uM)", y = "Density")
  
  ggsave(paste(prout, "luc.png", sep = ""),  width = 8, height = 7, dpi = 300, bg="transparent")
  
  
  din$IC50 = log10(din$IC50)
  mu <- ddply(din, "Assays", summarise, grp.mean=mean(IC50))
  ggplot(din, aes(x=IC50, color=Assays, fill=Assays)) +
    #geom_histogram(aes(y=..density..), position="identity", alpha=0.15)+
    geom_density(alpha=0.3)+
    xlim(-1.5,3)+
    theme(text = element_text(size=19))+
    geom_vline(data=mu, aes(xintercept=grp.mean, color=Assays),
               linetype=c("dashed"))+
    scale_fill_manual(values=c("#585858"), labels = c("IC50"))+
    scale_color_manual(values=c("#585858"), labels = c("IC50")) + 
    labs(title="",x="log(IC50) (uM)", y = "Density")
  
  ggsave(paste(prout, "logluc.png", sep = ""),  width = 8, height = 7, dpi = 300, bg="transparent")
  
}





################
#     MAIN     #
################

args <- commandArgs(TRUE)
pAC50all = args[1]
prout = args[2]

#pAC50all = "/home/borrela2/interference/spDataAnalysis/AC50_all"
#prout = "/home/borrela2/interference/spDataAnalysis/hist/"

dAC50 = read.csv(pAC50all, sep = "\t", header = TRUE)
rownames(dAC50) = dAC50[,1]

#dred = dAC50[,c("CASID", "hepg2_cell_red_n", "hepg2_med_red_n", "hek293_cell_red_n", "hek293_med_red_n")]
combinered(dAC50, prout)
combinegreen(dAC50, prout)
combineblue(dAC50, prout)
combineluc(dAC50, prout)
