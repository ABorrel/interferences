#!/usr/bin/env Rscript
library(ggplot2)
library(gridExtra)
library(grid)


generateReponseCurves = function(dac50, nbNA, praw, prout){
  
  # case where we have 3 IC50
  #pdf(paste(prout, "IC50_", nbNA,".pdf", sep = ""), 10,15)
  
  lplot = NULL
  iplot = 0
  for(i in seq(1, dim(dac50)[1])){
  #for(i in seq(1, 100)){
    if (length(which(is.na(dac50[i,]))) == nbNA){
      
      nameCpd = dac50[i,1]
      datraw = read.table(paste(praw, dac50[i,1], sep = ""), sep = "\t", header = TRUE)
      
      # table AC50 #
      ##############
      tableAC50 = t(dac50[i,])
      rownames(tableAC50) = colnames(dac50[i,])
      colnames(tableAC50) = "AC50"
      # del first line
      #tableAC50 = tableAC50[-which(rownames(tableAC50) == "CAS"),]
      
      
      # format table for plot
      tableAC50 = tableGrob(tableAC50, theme = ttheme_default(base_size = 7))
      

      # Legend #
      ##########
      plegend = ggplot(data = datraw, 
          aes(x = log10(CONC), y = DATA, color = Fluorophores))+
         geom_line()
      # legend
      tmp <- ggplot_gtable(ggplot_build(plegend))
      leg <- which(sapply(tmp$grobs, function(x) x$name) ==  "guide-box")
      legend <- tmp$grobs[[leg]]
      

      # plot #
      ########
      tp = ggplot(data = datraw, 
        aes(x = log10(CONC), y = DATA, color = Fluorophores)) +
        theme(legend.position="none")+
        ggtitle(nameCpd)+
        geom_line()
      
      
      tp2 = grid.arrange(tp, arrangeGrob(legend, tableAC50),  layout_matrix = matrix(c(1, 1, 2), ncol = 3))
      iplot = iplot + 1
      lplot[[iplot]] = tp2
      #if (length(lplot) == 12){
      #  #tg = grid.arrange(lplot[[1]], lplot[[2]], lplot[[3]], lplot[[4]], lplot[[5]], lplot[[6]], lplot[[7]], lplot[[8]],lplot[[9]], lplot[[10]],lplot[[11]], lplot[[12]], nrow = 4, ncol=3)
      #  tg = grid.arrange(grobs = lplot, nrow=4, ncol=3)
      #  print (tg)
      #  lplot = NULL
      #  iplot = 0
      #}
    }
  }
  if(dim(tableAC50)[1] > 4){
    ml = marrangeGrob(lplot, nrow = 3, ncol = 3)
  } else{
    ml = marrangeGrob(lplot, nrow = 4, ncol = 3)
  }
  ggsave(paste(prout, "IC50_", nbNA,".pdf", sep = ""), ml, width = 12, height = 15)
  #if(!is.null(lplot)){
  #  tg = grid.arrange(grobs = lplot, nrow=4, ncol=3)
  #  print (tg)
  #}
  #dev.off()
}



################
#     MAIN     #
################

args <- commandArgs(TRUE)
praw = args[1]
pAC50 = args[2]
prout = args[3]

#praw = "./../../spDataAnalysis/responseCurve/"
#pAC50 = "./../../spDataAnalysis/AC50_sample"
#prout = "./../../spDataAnalysis/"


dac50 = read.csv(pAC50, sep = "\t", header = TRUE)
rownames(dac50) = dac50[,1]

if(!is.integer(which(colnames(dac50) == "control"))){
  dac50 = dac50[,-which(colnames(dac50) == "control")]
}


if(!is.integer(which(rownames(dac50) == ""))){
  dac50 = dac50[-which(rownames(dac50) == ""),]
}


for(i in seq(0,dim(dac50)[2]-1)){
  generateReponseCurves(dac50, i, praw, prout)
}
