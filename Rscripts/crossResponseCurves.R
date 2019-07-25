#!/usr/bin/env Rscript
library(ggplot2)
library(gridExtra)
library(grid)


generateReponseCurvesSpFluoCross = function(pac501, pac502, chead, praw1, praw2, prout, fluo){
  
  # case where we have 3 IC50
  #pdf(paste(prout, "IC50_", nbNA,".pdf", sep = ""), 10,15)
  
  lplot = NULL
  iplot = 0
  
  dac501 = openAC50(pac501, "HepG2", chead)
  dac502 = openAC50(pac502, "Hek293", chead)
  
  
  lcompound = unique(dac501[,1], dac502[,1]) 
  # loop on compound
  for(nameCpd in lcompound){
    #nameCpd = "632-69-9"#red
    #nameCpd = "396-01-0" #blue
    #nameCpd = "2321-07-5" #green
    
    #for(i in seq(1, 100)){
    #for(i in seq(1, 20)){
    print(nameCpd)
    datraw1 = uploadRawData(nameCpd, chead, "HepG2", praw1)
    datraw2 = uploadRawData(nameCpd, chead, "Hek293", praw2)
    datraw = rbind(datraw1, datraw2)
    
    
    # table AC50 #
    ##############
    dac50 = cbind(dac501[c(nameCpd),],dac502[c(nameCpd),])
    dac50 = dac50[,unique(colnames(dac50))]
    tableAC50 = t(dac50)
    rownames(tableAC50) = colnames(dac50)
    # add curve in the table
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
    
    # case of fluorescine
    if(nameCpd == "2321-07-5" & fluo == "green"){
      datraw = na.omit(datraw)
      datraw$Fluorophores = as.character(datraw$Fluorophores)
      datraw[datraw == "HepG2-med_green_n"] = "HepG2 cell free"
      datraw[datraw == "Hek293-med_green_n"] = "Hek293 cell free"
      datraw[datraw == "HepG2-cell_green_n"] = "HepG2 cell based"
      datraw[datraw == "Hek293-cell_green_n"] = "Hek293 cell based"
      datraw[is.na(datraw)] = "NA"
      datraw = as.data.frame(datraw)
      datraw$CONC = as.double(datraw$CONC)
      datraw$DATA = as.double(datraw$DATA)
      
      rownames(tableAC50) = c("CAS", "HepG2 cell free", "HepG2 cell based", "Hek293 cell free", "Hek293 cell based")
    }
    
    if(nameCpd == "396-01-0" & fluo == "blue"){
      datraw = na.omit(datraw)
      datraw$Fluorophores = as.character(datraw$Fluorophores)
      datraw[datraw == "HepG2-med_blue_n"] = "HepG2 cell free"
      datraw[datraw == "Hek293-med_blue_n"] = "Hek293 cell free"
      datraw[datraw == "HepG2-cell_blue_n"] = "HepG2 cell based"
      datraw[datraw == "Hek293-cell_blue_n"] = "Hek293 cell based"
      datraw[is.na(datraw)] = "NA"
      datraw = as.data.frame(datraw)
      datraw$CONC = as.double(datraw$CONC)
      datraw$DATA = as.double(datraw$DATA)
      
      rownames(tableAC50) = c("CAS", "HepG2 cell free", "HepG2 cell based", "Hek293 cell free", "Hek293 cell based")
    }
    
    if(nameCpd == "632-69-9" & fluo == "red"){
      datraw = na.omit(datraw)
      datraw$Fluorophores = as.character(datraw$Fluorophores)
      datraw[datraw == "HepG2-med_red_n"] = "HepG2 cell free"
      datraw[datraw == "Hek293-med_red_n"] = "Hek293 cell free"
      datraw[datraw == "HepG2-cell_red_n"] = "HepG2 cell based"
      datraw[datraw == "Hek293-cell_red_n"] = "Hek293 cell based"
      datraw[is.na(datraw)] = "NA"
      datraw = as.data.frame(datraw)
      datraw$CONC = as.double(datraw$CONC)
      datraw$DATA = as.double(datraw$DATA)
      
      rownames(tableAC50) = c("CAS", "HepG2 cell free", "HepG2 cell based", "Hek293 cell free", "Hek293 cell based")
    }
    
    
    # extract only active
    i = 2
    imax = dim(tableAC50)[1]
    flagAct = 0
    while(i <= imax){
      #print(tableAC50[i,1])
      #print(tableAC50[i,2])
      if(!is.na(tableAC50[i,1]) && tableAC50[i,2] != "4"){
        flagAct = 1
        i = imax
      }
      i = i + 1
    }
    
    if(flagAct == 1){
      # format table for plot
      tableAC50 = tableGrob(tableAC50, theme = ttheme_default(base_size = 10))
      
      if(fluo == "blue"){
        # Legend #
        ##########
        plegend = ggplot(data = datraw, 
                         aes(x = log10(CONC), y = DATA, color = Fluorophores))+
          scale_color_manual(values=c("#0000ff", "#00c5ff", "#0000ff", "#00c5ff"))+
          scale_linetype_manual(values=c("solid", "solid", "dotted","dotted" ))+
          theme(legend.text = element_text(size = 12), legend.title = element_blank())+
          geom_line(aes(linetype = Fluorophores))
        # legend
        tmp <- ggplot_gtable(ggplot_build(plegend))
        leg <- which(sapply(tmp$grobs, function(x) x$name) ==  "guide-box")
        legend <- tmp$grobs[[leg]]
        
        
        # plot #
        ########
        print(datraw)
        tp = ggplot(data = datraw, 
                    aes(x = log10(CONC), y = DATA, color = Fluorophores)) +
          theme(legend.position="none")+
          ggtitle(nameCpd)+
          scale_color_manual(values=c("#0000ff", "#00c5ff", "#0000ff", "#00c5ff"))+
          scale_linetype_manual(values=c("solid", "solid", "dotted","dotted" ))+
          labs(x = "Concentration (Log(uM))", y = "Response (%)") + 
          theme(axis.text.y = element_text(size = 12, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 12, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 14, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 14, hjust = 0.5, vjust =0.1))+
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
        
        if(nameCpd == "396-01-0"){
          ggsave(paste(prout, "tiramterene.png", sep = ""), tp2, width = 8, height = 6, dpi = 300)
        }
        
        
      } else if (fluo == "red"){
        plegend = ggplot(data = datraw, 
                         aes(x = log10(CONC), y = DATA, color = Fluorophores))+
          scale_color_manual(values=c("#ff0000", "#c70808","#ff0000", "#c70808"))+
          scale_linetype_manual(values=c("solid", "solid", "dotted","dotted" ))+
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
          theme(legend.position="none")+
          ggtitle(nameCpd)+
          scale_color_manual(values=c("#ff0000", "#c70808","#ff0000", "#c70808"))+
          scale_linetype_manual(values=c("solid", "solid", "dotted","dotted" ))+
          labs(x = "Concentration (Log(uM))", y = "Response (%)") +
          theme(axis.text.y = element_text(size = 12, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 12, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 14, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 14, hjust = 0.5, vjust =0.1))+
          geom_line(aes(linetype = Fluorophores))
        
        
        tp2 = grid.arrange(tp, arrangeGrob(legend, tableAC50),  layout_matrix = matrix(c(1, 1, 2), ncol = 3))
        
        if(nameCpd == "632-69-9"){
          ggsave(paste(prout, "rosebangal_analogue.png", sep = ""), tp2, width = 8, height = 6, dpi = 300)
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
      }else if (fluo == "green"){
        
        plegend = ggplot(data = datraw, 
                         aes(x = log10(CONC), y = DATA, color = Fluorophores))+
          scale_color_manual(values=c("#008000", "#21dd21","#008000", "#21dd21"))+
          scale_linetype_manual(values=c("solid", "solid", "dotted","dotted" ))+
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
          theme(legend.position="none")+
          ggtitle(nameCpd)+
          scale_color_manual(values=c("#008000", "#21dd21","#008000", "#21dd21"))+
          scale_linetype_manual(values=c("solid", "solid", "dotted","dotted" ))+
          labs(x = "Concentration (Log(uM))", y = "Response (%)") +
          theme(axis.text.y = element_text(size = 12, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 12, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 14, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 14, hjust = 0.5, vjust =0.1))+
          geom_line(aes(linetype = Fluorophores))
        
        
        tp2 = grid.arrange(tp, arrangeGrob(legend, tableAC50),  layout_matrix = matrix(c(1, 1, 2), ncol = 3))
        iplot = iplot + 1
        
        if(nameCpd == "2321-07-5"){
          ggsave(paste(prout, "luciferin.png", sep = ""), tp2, width = 8, height = 6, dpi = 300)
        }
        

     
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
  }
  if(dim(tableAC50)[1] > 4){
    ml = marrangeGrob(lplot, nrow = 3, ncol = 3)
  } else{
    ml = marrangeGrob(lplot, nrow = 4, ncol = 3)
  }
  ggsave(paste(prout, "IC50.pdf", sep = ""), ml, width = 26, height = 15)
  #if(!is.null(lplot)){
  #  tg = grid.arrange(grobs = lplot, nrow=4, ncol=3)
  #  print (tg)
  #}
  #dev.off()
}



uploadRawData = function(nameCpd, chead, prefix, praw){
  
  print(paste(praw, nameCpd, sep = ""))
  
  datraw = read.table(paste(praw, nameCpd, sep = ""), sep = "\t", header = TRUE)
  #reduce matrix
  ikeep = NULL
  for(h in chead){
    ikeep = append(ikeep, which(datraw[,c("Fluorophores")] == h))
  }
  datraw = datraw[ikeep,]
  datraw$Fluorophores <- as.character(datraw$Fluorophores)
  
  
  i = 1
  imax = dim(datraw)[1]
  for(i in seq(1,imax)){
    datraw$Fluorophores[i] = paste(prefix, "-", datraw$Fluorophores[i], sep = "")
    
  }
  datraw$Fluorophores <- as.factor(datraw$Fluorophores)
  return (datraw)
}
  
  

openAC50 = function(pin, prefix, chead){
  
  dac50 = read.csv(pin, sep = "\t", header = TRUE)
  rownames(dac50) = dac50[,1]
  
  #remove control
  if(!is.integer(which(colnames(dac50) == "control"))){
    dac50 = dac50[,-which(colnames(dac50) == "control")]
  }
  
  
  if(!is.integer(which(rownames(dac50) == ""))){
    dac50 = dac50[-which(rownames(dac50) == ""),]
  }
  
  
  icol = NULL
  for (h in chead){
    icol = append(icol, which(colnames(dac50) == h))
  }
  
  lcolname = NULL
  for (h in chead){
    lcolname = append(lcolname, paste(prefix, "-", h, sep = ""))
  }
  lcolname = append("CAS", lcolname)
  
  icol = append(1,icol)
  dac50 = dac50[,icol]
  
  colnames(dac50) = lcolname
  
  return(dac50)
}


################
#     MAIN     #
################

args <- commandArgs(TRUE)
prawHepg2 = args[1]
prawHek293 = args[2]
pAC50Hepg2 = args[3]
pAC50Hek293 = args[4]
prout = args[5]


#prawHepg2 = "/home/borrela2/interference/spDataAnalysis/tox21-spec-hepg2-p1/responseCurve/"
#prawHek293 = "/home/borrela2/interference/spDataAnalysis/tox21-spec-hek293-p1/responseCurve/"
#pAC50Hepg2 = "/home/borrela2/interference/spDataAnalysis/tox21-spec-hepg2-p1/AC50_sample"
#pAC50Hek293  = "/home/borrela2/interference/spDataAnalysis/tox21-spec-hek293-p1/AC50_sample"
#prout = "/home/borrela2/interference/spDataAnalysis/crossCurvesResponse/"
  

# AC50
generateReponseCurvesSpFluoCross(pAC50Hepg2, pAC50Hek293, c("cell_blue_n", "med_blue_n"), prawHepg2, prawHek293, paste(prout, "blue_n_", sep = ""), "blue")
generateReponseCurvesSpFluoCross(pAC50Hepg2, pAC50Hek293, c("med_green_n", "cell_green_n"), prawHepg2, prawHek293, paste(prout, "green_n_", sep = ""), "green")
generateReponseCurvesSpFluoCross(pAC50Hepg2, pAC50Hek293, c("med_red_n", "cell_red_n"), prawHepg2, prawHek293, paste(prout, "red_n_", sep = ""), "red")