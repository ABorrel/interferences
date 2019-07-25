#!/usr/bin/env Rscript
library(ggplot2)
library(gridExtra)
library(grid)

source("~/development/Rglobal/source/dataManager.R")

generateReponseCurves = function(dac50, nbNA, praw, prout){
  
  # case where we have 3 IC50
  #pdf(paste(prout, "IC50_", nbNA,".pdf", sep = ""), 10,15)
  
  lplot = NULL
  iplot = 0
  for(i in seq(1, dim(dac50)[1])){
    #i = 1403
  #for(i in seq(1, 100)){
    if (length(which(is.na(dac50[i,]))) == nbNA){
      nameCpd = dac50[i,1]
      datraw = read.table(paste(praw, dac50[i,1], sep = ""), sep = "\t", header = TRUE)
      # reduce datraw
      datrawfiltered = NULL
      for(j in seq(1, dim(datraw)[1])){
        if (!is.integer0(which(colnames(dac50) == datraw[j,c("Fluorophores")]))){
          datrawfiltered = rbind(datrawfiltered, datraw[j,])
        }
      }
      colnames(datrawfiltered) = colnames(datraw)
      datraw = datrawfiltered
      
      
      # table AC50 #
      ##############
      tableAC50 = t(as.data.frame(dac50[i,]))
      rownames(tableAC50) = colnames(dac50[i,])
      ccurves = NULL
      for(nameTest in rownames(tableAC50)){
        if(nameTest == "CAS"){
          ccurves = append(ccurves, "")
        }else{
          ccurves = append(ccurves, datraw[which(datraw$Fluorophores == nameTest)[1],c("CurveType")])
        }
      }      
      
      tableAC50 = cbind(tableAC50, ccurves)
      colnames(tableAC50) = c("AC50", "Curve")
      # del first line
      #tableAC50 = tableAC50[-which(rownames(tableAC50) == "CAS"),]
      
      
      # format table for plot
      tableAC50 = tableGrob(tableAC50, theme = ttheme_default(base_size = 12))

      
      if(length(unique(datraw$Fluorophores)) == 3){
        
        # Legend #
        ##########
        plegend = ggplot(data = datraw, 
                         aes(x = log10(CONC), y = DATA, color = Fluorophores))+
          scale_color_manual(values=c("#282828", "#282828", "#282828"))+
          scale_linetype_manual(values=c("solid", "dotted", "longdash" ))+
          theme(legend.text = element_text(size = 12), legend.title = element_blank())+
          geom_line(aes(linetype = Fluorophores))
        # legend
        tmp <- ggplot_gtable(ggplot_build(plegend))
        leg <- which(sapply(tmp$grobs, function(x) x$name) ==  "guide-box")
        legend <- tmp$grobs[[leg]]
        
        # plot #
        ########
        tp = ggplot(data = datraw, 
                    aes(x = log10(CONC), y = DATA, color = Fluorophores)) +
          scale_color_manual(values=c("#282828", "#282828", "#282828"))+
          scale_linetype_manual(values=c("solid", "dotted", "longdash" ))+
          labs(x = "Concentration (Log(uM))", y = "Response (%)") +
          theme(axis.text.y = element_text(size = 12, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 12, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 14, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 14, hjust = 0.5, vjust =0.1))+
          theme(legend.position="none")+
          ggtitle(nameCpd)+
          geom_line(aes(linetype = Fluorophores))
        
      }else{
        
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
      }
      
      tp2 = grid.arrange(tp, arrangeGrob(legend, tableAC50),  layout_matrix = matrix(c(1, 1, 2), ncol = 3))
      
      if(nameCpd == "775304-57-9"){
        ggsave(paste(prout, "ataluren.png", sep = ""), tp2, width = 7, height = 6, dpi = 300)
        #return()
      }
      
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
  if(!is.null(lplot)){
    if(dim(dac50)[2] > 4){
      ml = marrangeGrob(lplot, nrow = 3, ncol = 3)
    } else{
      ml = marrangeGrob(lplot, nrow = 4, ncol = 3)
    }
    ggsave(paste(prout, "IC50_", nbNA,".pdf", sep = ""), ml, width = 21, height = 18)
    #if(!is.null(lplot)){
    #  tg = grid.arrange(grobs = lplot, nrow=4, ncol=3)
    #  print (tg)
    #}
    #dev.off()
  }
}

generateReponseCurvesSpFluo = function(dac50, chead, praw, prout, fluo){
  
  # case where we have 3 IC50
  #pdf(paste(prout, "IC50_", nbNA,".pdf", sep = ""), 10,15)
  
  lplot = NULL
  iplot = 0
  icol = NULL
  for (h in chead){
    icol = append(icol, which(colnames(dac50) == h))
  }
  icol = append(1,icol)
  dac50 = dac50[,icol]
  #print(head(dac50))
  for(i in seq(1, dim(dac50)[1])){
  #for(i in seq(1, 20)){
    nameCpd = dac50[i,1]
    datraw = read.table(paste(praw, dac50[i,1], sep = ""), sep = "\t", header = TRUE)
    #reduce matrix
    ikeep = NULL
    for(h in chead){
      ikeep = append(ikeep, which(datraw[,c("Fluorophores")] == h))
    }
    datraw = datraw[ikeep,]
    # table AC50 #
    ##############
    tableAC50 = t(dac50[i,])
    rownames(tableAC50) = colnames(dac50[i,])
    ccurves = NULL
    for(nameTest in rownames(tableAC50)){
      if(nameTest == "CAS"){
        ccurves = append(ccurves, "")
      }else{
        ccurves = append(ccurves, datraw[which(datraw$Fluorophores == nameTest)[1],c("CurveType")])
      }
    } 
    
    tableAC50 = cbind(tableAC50, ccurves)
    colnames(tableAC50) = c("AC50", "Curve")
    # del first line
    #tableAC50 = tableAC50[-which(rownames(tableAC50) == "CAS"),]
      
      
    # format table for plot
    tableAC50 = tableGrob(tableAC50, theme = ttheme_default(base_size = 6.5))
      
    #print(head(tableAC50))
      
    if(fluo == "blue"){
      # Legend #
      ##########
      plegend = ggplot(data = datraw, 
                       aes(x = log10(CONC), y = DATA, color = Fluorophores))+
        scale_color_manual(values=c("#0000ff", "#00c5ff", "#0000ff", "#00c5ff"))+
        scale_linetype_manual(values=c("solid", "solid", "dotted","dotted" ))+
        geom_line(aes(linetype = Fluorophores))
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
        scale_color_manual(values=c("#0000ff", "#00c5ff", "#0000ff", "#00c5ff"))+
        scale_linetype_manual(values=c("solid", "solid", "dotted","dotted" ))+
        labs(x = "Concentration (LogM)", y = "Response (%)") +
        geom_line(aes(linetype = Fluorophores))
      
      
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
    } else if (fluo == "red"){
      plegend = ggplot(data = datraw, 
        aes(x = log10(CONC), y = DATA, color = Fluorophores))+
        scale_color_manual(values=c("#ff0000", "#c70808", "#ff0000", "#c70808"))+
        scale_linetype_manual(values=c("solid", "solid", "dotted","dotted"))+
        geom_line(aes(linetype = Fluorophores))
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
        scale_color_manual(values=c("#ff0000", "#c70808", "#ff0000", "#c70808"))+
        scale_linetype_manual(values=c("solid", "solid", "dotted","dotted" ))+
        labs(x = "Concentration (LogM)", y = "Response (%)") +
        geom_line(aes(linetype = Fluorophores))
      
      
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
    }else if (fluo == "green"){
      plegend = ggplot(data = datraw, 
                       aes(x = log10(CONC), y = DATA, color = Fluorophores))+
        scale_color_manual(values=c("#008000", "#21dd21", "#008000", "#21dd21"))+
        scale_linetype_manual(values=c("solid", "solid", "dotted","dotted" ))+
        geom_line(aes(linetype = Fluorophores))
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
        scale_color_manual(values=c("#008000", "#21dd21", "#008000", "#21dd21"))+
        scale_linetype_manual(values=c("solid", "solid", "dotted","dotted" ))+
        labs(x = "Concentration (LogM)", y = "Response (%)") +
        geom_line(aes(linetype = Fluorophores))
      
      
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
  if(!is.null(lplot)){
    if(dim(tableAC50)[1] > 4){
      ml = marrangeGrob(lplot, nrow = 3, ncol = 3)
    } else{
      ml = marrangeGrob(lplot, nrow = 4, ncol = 3)
    }
    ggsave(paste(prout, "IC50.pdf", sep = ""), ml, width = 21, height = 18)
    #if(!is.null(lplot)){
    #  tg = grid.arrange(grobs = lplot, nrow=4, ncol=3)
    #  print (tg)
    #}
    #dev.off()
  }
}


################
#     MAIN     #
################

args <- commandArgs(TRUE)
praw = args[1]
pAC50 = args[2]
prout = args[3]

#praw = "/home/borrela2/interference/spDataAnalysis/tox21-spec-hek293-p1/responseCurve/"
#pAC50 = "/home/borrela2/interference/spDataAnalysis/tox21-spec-hek293-p1/AC50_sample"
#prout = "/home/borrela2/interference/spDataAnalysis/tox21-spec-hek293-p1/"

#praw = "/home/borrela2/interference/spDataAnalysis/tox21-luc-biochem-p1/responseCurve/"
#pAC50 = "/home/borrela2/interference/spDataAnalysis/tox21-luc-biochem-p1/AC50_sample"
#prout = "/home/borrela2/interference/spDataAnalysis/tox21-luc-biochem-p1/"


dac50 = read.csv(pAC50, sep = "\t", header = TRUE)
rownames(dac50) = dac50[,1]

#remove control
if(!is.integer(which(colnames(dac50) == "control"))){
  dac50 = dac50[,-which(colnames(dac50) == "control")]
}


if(!is.integer(which(rownames(dac50) == ""))){
  dac50 = dac50[-which(rownames(dac50) == ""),]
}


if(dim(dac50)[2] >= 10){
  #dac50temp = dac50[,c("CAS", "med_green", "med_blue", "med_red", "cell_green", "cell_blue", "cell_red")]
  #for(nbNA in seq(0,dim(dac50temp)[2]-1)){
  #  generateReponseCurves(dac50temp, nbNA, praw, prout)
  #}
  
  generateReponseCurvesSpFluo(dac50, c("med_blue_n", "cell_blue_n", "cell_blue", "med_blue"), praw, paste(prout, "blue_", sep = ""), "blue")
  generateReponseCurvesSpFluo(dac50, c("med_green_n", "cell_green_n", "med_green", "cell_green"), praw, paste(prout, "green_", sep = ""), "green")
  generateReponseCurvesSpFluo(dac50, c("med_red_n", "cell_red_n", "med_red", "cell_red"), praw, paste(prout, "red_", sep = ""), "red")
  
  
  dac50temp = dac50[,c("CAS", "med_green_n", "med_blue_n", "med_red_n", "cell_green_n", "cell_blue_n", "cell_red_n")]
  for(nbNA in seq(0,dim(dac50temp)[2]-1)){
    generateReponseCurves(dac50temp, nbNA, praw, paste(prout, "N_", sep = ""))
  }
}else{
  for(nbNA in seq(0,dim(dac50)[2]-1)){
    generateReponseCurves(dac50, nbNA, praw, prout)  
    }
}