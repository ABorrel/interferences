#!/usr/bin/env Rscript
library(ggplot2)
library(arcdiagram)
source("arcBands.R")
source("arcDiagram.R")
source("cardMatrix.R")

piesDraw = function(dAssays, prout){
  
  
  # generate for all colnames
  lcharacAssay = colnames(dAssays)
  
  for(characAssay in lcharacAssay){
    dtemp = rownames(dAssays)
    dtemp = cbind(dtemp, as.character(dAssays[,characAssay]))
    colnames(dtemp) = c("ID", "Classes")
    print(dtemp)
    dtemp = as.data.frame(dtemp)
    if(length(unique(dtemp[,2]))<100){
      
      p <- ggplot(dtemp, aes(x=factor(1), fill=Classes))+
        coord_polar(theta = "y")
      
      p <- p +
        geom_bar(color='black')+
        guides(fill=guide_legend(override.aes=list(colour=NA)))+
        theme(axis.text.x=element_text(color='black'))
      
      
      ggsave(filename = paste(prout, "pie_", characAssay, ".png", sep = ""), dpi = 300, width = 10, height = 7)
    }
  }
}


histAC50 = function(dAssays, dAC50, typesplit){
  
  # prepare data
  ltypesplit = unique(dAssays[,which(colnames(dAssays) == typesplit)])
  ltypesplit = ltypesplit[-which(ltypesplit == "")]
  
  dcount = NULL
  nbchemical = dim(dAC50)[1]
  nbAssay = dim(dAC50)[2]
  
  
  lAC50 = list()
  nbelemt = length(ltypesplit)
  for(i in seq(1,nbelemt)){
    lAC50[[i]] = vector() 
  }
  
  names(lAC50) = ltypesplit
  
  i = 2# remove chemical ID
  while(i <= nbchemical){
    j = 2
    print(i)
    print(nbchemical)
    while(j <= nbAssay){
      if (is.na(dAC50[i,j])){
        j = j + 1
      }else{
        nameAssay = colnames(dAC50)[j]
        typeAssay = dAssays[which(rownames(dAssays) == nameAssay),typesplit]
        if(typeAssay != ""){
          lAC50[[which(names(lAC50) == typeAssay)]] = append(lAC50[[which(names(lAC50) == typeAssay)]], dAC50[i,j]) 
        }
      }
      j = j + 1
    }
    i = i + 1
  }
  
  
  for(namek in names(lAC50)){
    if (is.null(lAC50[[which(names(lAC50) == namek)]])){
      print(namek)
    }else{
      AC50s = as.data.frame(lAC50[[which(names(lAC50) == namek)]])
      colnames(AC50s) = "AC50"
      ggplot(AC50s, aes(x=AC50)) + 
        geom_histogram(aes(y=..density..),     # Histogram with density instead of count on y-axis
                       binwidth=.5,
                       colour="black", fill="white") +
        labs(x = "AC50", y = "Frequencies") + 
        #xlim (c(-2.5, 2.5))+
        #ylim(c(0, 0.65))+
        theme(axis.text.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 25, hjust = 0.5, vjust =0.1))+
        geom_density(alpha=.2, fill="#FF6666")
      ggsave(paste(prout, namek, "_AC50.png", sep = ""), dpi = 300, width = 8, height = 7)
      
      AC50s = -log10(AC50s)
      ggplot(AC50s, aes(x=AC50)) + 
        geom_histogram(aes(y=..density..),     # Histogram with density instead of count on y-axis
                       binwidth=.5,
                       colour="black", fill="white") +
        labs(x = "-log(AC50)", y = "Frequencies") + 
        #xlim (c(-2.5, 2.5))+
        #ylim(c(0, 0.65))+
        theme(axis.text.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 25, hjust = 0.5, vjust =0.1))+
        geom_density(alpha=.2, fill="#FF6666")
      ggsave(paste(prout, namek, "_logAC50.png", sep = ""), dpi = 300, width = 8, height = 7)
    }
  }
}



corrAssaysGlobal = function(dAC50, prout){
  
  dAC50 = dAC50[,-1] #remove ID col
  tabout = data.frame()
  # pdf for correlation
  i = 1
  imax = dim(dAC50)[2]
  while(i <= imax){
    j = i
    while(j <= imax){
      tempvect = cbind(dAC50[,i], dAC50[,j])
      tempvect = na.omit(tempvect)
      valcor = cor(tempvect[,1], tempvect[,2])
      tabout[i,j] = valcor
      j = j + 1
    }
    i = i + 1
  }
  colnames(tabout) = colnames(dAC50)
  rownames(tabout) = colnames(dAC50)
  print (tabout)
  write.csv(tabout, paste(prout, "corTable.csv", sep = ""))
}


corrAssaysType = function(dAC50, dAssays, prout, analysis){
  
  dAC50 = dAC50[,-1] #remove ID col
  #print(colnames(dAC50))
  # prepare data
  lanalysis = unique(dAssays[,which(colnames(dAssays) == analysis)])
  #lanalysis[which(lanalysis == "")] = "UNK"
  #print(lanalysis)
  
  pdf(paste(prout, analysis, ".pdf", sep = ""), width = 25, height = 25)
  for (analysis_sub in lanalysis){
    lassay = rownames(dAssays)[which(dAssays[,analysis] == analysis_sub)]
    print (length(lassay))
    dAC50temp = NULL
    lassaytemp = NULL
    for(assay in lassay){
      if (!is.integer0(which(colnames(dAC50) == assay))){
        dAC50temp = cbind(dAC50temp, dAC50[,which(colnames(dAC50) == assay)])
        lassaytemp = append(lassaytemp, assay)
      }
    }
    lassay = lassaytemp# to remove col no selected
    if(length(lassay) != 0){
      print(length(lassay))
      print (dim(dAC50temp))
      tabout = data.frame()
      # pdf for correlation
      i = 1
      imax = dim(dAC50temp)[2]
      while(i <= imax){
        j = i
        while(j <= imax){
          tempvect = cbind(dAC50temp[,i], dAC50temp[,j])
          tempvect = na.omit(tempvect)
          valcor = cor(tempvect[,1], tempvect[,2])
          tabout[i,j] = valcor
          j = j + 1
        }
        i = i + 1
      }
      print (dim(tabout))
      colnames(tabout) = lassay
      rownames(tabout) = lassay
      if (dim(tabout)[1] != 1){
        cardMatrixCorpdf(tabout, 6, analysis_sub)  
      }
      #print (tabout)
      }
  }
dev.off()  
}



ArcDiagramDraw = function(dAssays, dAC50, cutoff, typeArcDia, prout){
  
  # prepare data
  ltypeArc = unique(dAssays[,which(colnames(dAssays) == typeArcDia)])
  ltypeArc = ltypeArc[-which(ltypeArc == "")]
  print(ltypeArc)
  
  dcount = NULL
  nbchemical = dim(dAC50)[1]
  nbAssay = dim(dAC50)[2]
  
  ## define structure ##
  ######################
  # 1. count active inactive by analysis
  for(typedia in ltypeArc){
    countChem = c(0,0)
    dcount = rbind(dcount, countChem)
  }
  rownames(dcount) = ltypeArc
  colnames(dcount) = c("Active", "Inactive")

  #2. combination edge
  i = 1
  edAct = NULL
  edInact = NULL
  edContract = NULL

  while(i < length(ltypeArc)){
    j = i + 1
    while(j <= length(ltypeArc)){
      cadd = c(as.character(ltypeArc[i]), as.character(ltypeArc[j]))
      edAct = rbind(edAct, cadd)
      edInact = rbind(edInact, cadd)
      edContract = rbind(edAct, cadd)
      j = j + 1
    }
    i = i + 1
  }
  # active
  rownames(edAct) = seq(1, dim(edAct)[1])
  edAct = as.data.frame(edAct)
  edAct = cbind(edAct, rep(0, dim(edAct)[1]))
  
  # inactive
  rownames(edInact) = seq(1, dim(edInact)[1])
  edInact = as.data.frame(edInact)
  edInact = cbind(edInact, rep(0, dim(edInact)[1]))
  
  # contract
  rownames(edContract) = seq(1, dim(edContract)[1])
  edContract = as.data.frame(edContract)
  edContract = cbind(edContract, rep(0, dim(edContract)[1]))
  
  
  # analyse data  
  i = 2# remove chemical ID
  while(i <= nbchemical){
    j = 2
    lactive = rep(0,length(ltypeArc))
    names(lactive) = ltypeArc
    while(j <= nbAssay){
      if (is.na(dAC50[i,j])){
        j = j + 1
      }else{
        nameAssay = colnames(dAC50)[j]
        typeAssay = dAssays[which(rownames(dAssays) == nameAssay),typeArcDia]
        #print(paste(nameAssay, typeAssay))
        #print(dAssays[nameAssay, typeArcDia])
        
        typedia = dAssays[nameAssay, typeArcDia]
        if(typedia != ""){
          if(dAC50[i,j] <= cutoff){
            dcount[which(rownames(dcount) == typedia),1] = dcount[which(rownames(dcount) == typedia),1] + 1
            lactive[typedia] = 1
          }else{
            dcount[which(rownames(dcount) == typedia),2] = dcount[which(rownames(dcount) == typedia),2] + 1
          }
        }
      }
      j = j + 1
    }
    i = i + 1
    
    
    k = 1
    while(k < length(lactive)){
      l = k + 1
      while(l <= length(lactive)){
        if (lactive[k] == 1 && lactive[l] == 1){
          #print(which(edAct[,2] == names(lactive)[k]))
          #print(which(edAct[,1] == names(lactive)[k]))
          #print("****")
          if(is.null(which(edAct[,1] == names(lactive)[k]))){
            iedAct = intersect(which(edAct[,2] == names(lactive)[k]), which(edAct[,1] == names(lactive)[l]))
          }else{
            iedAct = intersect(which(edAct[,1] == names(lactive)[k]), which(edAct[,2] == names(lactive)[l]))
          }
          edAct[iedAct,3] = edAct[iedAct,3] + 1
          
        }else if(lactive[k] == 0 && lactive[l] == 0){
          if(is.null(which(edInact[,1] == names(lactive)[k]))){
            iedInact = intersect(which(edInact[,2] == names(lactive)[k]), which(edInact[,1] == names(lactive)[l]))
          }else{
            iedInact = intersect(which(edInact[,1] == names(lactive)[k]), which(edInact[,2] == names(lactive)[l]))
          }
          edInact[iedInact,3] = edInact[iedInact,3] + 1
          
        }else{
          if(is.null(which(edContract[,1] == names(lactive)[k]))){
            iedContract = intersect(which(edContract[,2] == names(lactive)[k]), which(edContract[,1] == names(lactive)[l]))
          }else{
            iedContract = intersect(which(edContract[,1] == names(lactive)[k]), which(edContract[,2] == names(lactive)[l]))
          }
          edContract[iedContract,3] = edContract[iedContract,3] + 1
        }
        l = l + 1
      }
      k = k + 1
    }
  }
  #print (dcount)
  print(edAct)
  print(edInact)
  print(edContract)
  #col.bands = c("#4ECDC4", "#6d849c")
  
  lsum = dcount[,1] + dcount[,2]
  lwdsInact = edInact[,3]/lsum
  lwdsAct = edAct[,3]/lsum
  lwdsContract = edContract[,3]/lsum
  
  edActs = as.matrix(edAct[,-3])
  edInact = as.matrix(edInact[,-3])
  edContract = as.matrix(edContract[,-3])
  #dcount = dcount/100
  # draw arc diagram
 
  png(paste(prout, "arcDiagramAct_", typeArcDia, ".png", sep = ""), 2000, 1000)
  arcDiagram(edActs, 
           lwd=lwdsAct, mar=c(10,1,4,1))
  dev.off()

  png(paste(prout, "arcDiagramInact_", typeArcDia, ".png", sep = ""), 2000, 1000)
  arcDiagram(edInact, 
             lwd=lwdsInact, mar=c(10,1,4,1))
  dev.off()
  
  png(paste(prout, "arcDiagramContract_", typeArcDia, ".png", sep = ""), 2000, 1000)
  arcDiagram(edContract, 
             lwd=lwdsContract, mar=c(10,1,4,1))
  dev.off()
}




################
#     MAIN     #
################

args <- commandArgs(TRUE)
pAssays = args[1]
pAC50 = args[2]
pChem = args[3]
prout = args[4]


# path by pass
#pAssays = "~/interference/data/ToxCast_assay_master_2017-06-14.csv"
#pAC50 = "~/interference/data/ac50_Matrix_170614.csv"
#pChem = "~/interference/data/toxcast_chemicals_2017-06-14.csv"
#prout = "~/interference/ToxCast_analysis/"

distributionAnalysis = 1
arcdiagramAnalysis = 0
histAC50Analysis = 1
corAssay = 1

# open data
dAssays = read.csv(pAssays, sep = ",", header = TRUE)
rownames(dAssays) = dAssays[,1]
dAC50 = read.csv(pAC50, sep = ",", header = TRUE)
rownames(dAC50) = dAC50[,1]


# pie chart
if (distributionAnalysis == 1){
  prPie = paste(prout, "pieAssays/", sep = "")
  dir.create(prPie)
  piesDraw(dAssays, prPie)
}


#lanalysis = c("detection_technology_type", "detection_technology_type_sub", "detection_technology", "assay_source_name", "assay_design_type_sub")
lanalysis = c("detection_technology_type")
for(analysis in lanalysis){
  if(arcdiagramAnalysis == 1){
    pradcDia = paste(prout, "ArcDia/", sep = "")
    dir.create(pradcDia)
    ArcDiagramDraw(dAssays,dAC50, 10, analysis, prout) 
  }
  if(histAC50Analysis == 1){
    phistAC50 = paste(prout, "histAC50/", sep = "")
    dir.create(phistAC50)
    histAC50(dAssays, dAC50, analysis)
  }
}


if(corAssay == 1){
  #corrAssaysGlobal(dAC50, prout)
  pcorAssays = paste(prout, "corAssays/", sep = "")
  dir.create(pcorAssays)
  
  for(analysis in lanalysis){
    corrAssaysType(dAC50, dAssays, prout, analysis)
  }
}