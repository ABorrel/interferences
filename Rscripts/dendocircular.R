#!/usr/bin/env Rscript
library(ape)
library(phangorn)
library(ggtree)
library(factoextra)




dendogramCircle = function(ddes, daff, prout){
  
  #calibrate affinity for color
  daff = as.data.frame(daff)
  daff = cbind(rownames(daff), daff)
  
 
  matTrans1 <- scale(ddes)
  d <- dist(matTrans1, method = "euc")
  tupgma2 <- upgma(d, method="ward.D2")
  legend.title = "AC50"
  
  pfilout = paste(prout, "AC50_dendo.png", sep = "")
  t4 <- ggtree(tupgma2, layout="circular", size=1)
  t4 <- t4 %<+% daff + geom_text(aes(color=AC50, label=label, angle=angle), hjust=-0.5, size=1) +
    geom_tippoint(aes(color=AC50), alpha=0.75, size=0.5)+
    scale_color_continuous(low='red', high='lightgreen') +
    theme(legend.position="right")+
    theme(plot.margin = unit(c(0,0,0,0), "cm")) +
    geom_treescale(x = 1, y = 1, width = NULL, offset = NULL,
                   color = "white", linesize = 1E-100, fontsize = 1E-100)
  ggsave(pfilout, dpi=300, height = 10, width = 14)
  
  
  pfilout = paste(prout, "class_dendo.png", sep = "")
  t4 <- ggtree(tupgma2, layout="circular", size=1)
  t4 <- t4 %<+% daff + geom_text(aes(color=ASSAY_OUTCOME, label=label, angle=angle), hjust=-0.5, size=1) +
    geom_tippoint(aes(color=ASSAY_OUTCOME), alpha=0.75, size=0.5)+
    scale_color_discrete() +
    theme(legend.position="right")+
    theme(plot.margin = unit(c(0,0,0,0), "cm")) +
    geom_treescale(x = 1, y = 1, width = NULL, offset = NULL,
                   color = "white", linesize = 1E-100, fontsize = 1E-100)
  ggsave(pfilout, dpi=300, height = 10, width = 14)
  
  
  
}





dendogramCluster = function(ddes, daff, dcluster, prout){
  
  
  #calibrate affinity for color
  #minMatrix = min(daff)
  #maxMatrix = max(daff)
  
  #daff[which(daff == max(daff, na.rm = TRUE), arr.ind = TRUE)] = 3.5
  #daff[which(daff == min(daff, na.rm = TRUE), arr.ind = TRUE)] = -3.5
  
  daff = rbind(daff, daff[dim(daff)[1],])
  daff[dim(daff)[1],1] = -3.5
  daff = rbind(daff, daff[dim(daff)[1],])
  daff[dim(daff)[1],1] = 3.5
  
  
  ddes = rbind(ddes, ddes[dim(ddes)[1],])
  #ddes[dim(ddes)[1],1] = -3.5
  ddes = rbind(ddes, ddes[dim(ddes)[1],])
  #ddes[dim(ddes)[1],1] = 3.5
  #print(daff)
  #for(i in seq(1, dim(daff)[2])){
  #  daff[which(daff[,i] == max(daff[,i])), i] = maxMatrix
  #  daff[which(daff[,i] == min(daff[,i])), i] = minMatrix
  #}
  daff = as.data.frame(daff)
  
  
  matTrans1 <- scale(ddes)
  d <- dist(matTrans1, method = "euc")
  tupgma2 <- upgma(d, method="ward.D2")
  tupgma2 <- groupOTU(tupgma2, max(as.vector(dcluster[,2])))

  pfilout = paste(prout, "dendo_cluster_name.png", sep = "")
    
  t1 <- ggtree(tupgma2, layout="circular", size=0.8, aes(color=cluster))
  t1 <- t1 %<+% dcluster + geom_text(aes(color=cluster, label = cluster, angle=angle,  fontface="bold"), hjust=-1.5, size=1.2) +
    geom_tippoint(aes(color=cluster), alpha=0.75, size=1)+
    #scale_colour_gradientn(colours=rainbow(20))+
    theme(plot.margin = unit(c(0,0,0,0), "cm")) + 
    geom_treescale(x = 40, y = 40, width = NULL, offset = NULL,
                   color = "white", linesize = 1E-100, fontsize = 1E-100)
  
  
  #daff = daff[,-1]
  t2 = gheatmap(t1, daff, font.size = 2, offset = 3, width = 0.5, colnames_offset_x = 2, colnames_offset_y = 0 ,low='red', high='lightgreen')
    #scale_color_continuous(low='red', high='lightgreen') +
    #scale_colour_gradientn(colours = c("#DF541E", "#D4621A", "#CA7116", "#BF7F12", "#AA9D0B", "#A0AB07", "#95BA03", "#8BC900"),
    #                       values=c(-4.0,-3.0, -2.0, -1.0, 0, 1.0, 2.0, 3.0,4.0))
    #theme(legend.position="right")+

    #theme_tree()
  #print (t2)
  open_tree(t2, 15) %>% rotate_tree(15)
  ggsave(pfilout, dpi=300, height = 11, width = 11)
  
  
  
  #pfilout = paste(prout, "dendo_cluster.png", sep = "")
  
  #t1 <- ggtree(tupgma2, layout="circular", size=0.8, aes(color=cluster))
  #t1 <- t1 %<+% dcluster +
  #  geom_tippoint(aes(color=cluster), alpha=0.75, size=1)+
    #theme(legend.position="right")+
  #  theme(plot.margin = unit(c(0,0,0,0), "cm")) + 
  #  geom_treescale(x = 40, y = 40, width = NULL, offset = NULL,
  #                 color = "white", linesize = 1E-100, fontsize = 1E-100)
  
  #daff = daff[,-1]
  #t2 = gheatmap(t1, daff, font.size =0 , offset = 3, width = 0.2, colnames_offset_x = 2, colnames_offset_y = -0.5, low = "red", high = "lightgreen") +
    #scale_color_continuous(low='red', high='lightgreen') +
  #  theme_tree()
  #print (t2)
  #open_tree(t2, 15) %>% rotate_tree(15)
  #ggsave(pfilout, dpi=300, height = 11, width = 11)
}




dendoClusteringInter = function(ddes, daff, dcluster, prout, methdist, methagg){
  
  #daff = rbind(daff, daff[dim(daff)[1],])
  #daff[dim(daff)[1],1] = -3.5
  #aff = rbind(daff, daff[dim(daff)[1],])
  #aff[dim(daff)[1],1] = 3.5
  #ddes = rbind(ddes, ddes[dim(ddes)[1],])
  #ddes = rbind(ddes, ddes[dim(ddes)[1],])
  
  #ddes = ddes[1:100, 1:100]
  #outclust = hcut(ddes, k = max(as.vector(dcluster[,2])), hc_method = "ward.D2")
  #fviz_dend(outclust, show_labels = FALSE, type = "circular")
  #ggsave(paste(prout, "dendov1_", metcluster, "_",  metagregation, "_", metOptNB, ".png", sep = ""), dpi=300, height = 12, width = 13)
  
  
  daff = as.data.frame(daff)
  
  if(methdist == "None"){
    d = as.dist(ddes)
  }else{
    matTrans1 <- scale(ddes)
    d <- dist(matTrans1, method = methdist)    
  }

  tupgma2 <- upgma(d, method=methagg)
  tupgma2 <- groupOTU(tupgma2, max(as.vector(dcluster[,2])))
  
  pfilout = paste(prout, "dendo_cluster.png", sep = "")
  
  t1 <- ggtree(tupgma2, layout="circular", size=0.5, aes(color=Cluster))
  t1 <- t1 %<+% dcluster + geom_text(aes(color=Cluster, label = Cluster, angle=angle,  fontface="bold"), hjust=-0.5, size=1.2) +
    geom_tippoint(aes(color=Cluster), alpha=0.75, size=1)+
    scale_colour_gradientn(colours=rainbow(max(as.vector(dcluster[,2]))))+
    theme(plot.margin = unit(c(0,0,0,0), "cm")) + 
    geom_treescale(x = -10, y = 0, width = 3, offset = 0,
                   color = "white", linesize = 1E-100, fontsize = 1E-100)

    #t2 = gheatmap(t1, daff, font.size = 2, offset = 0, width = 0.5, colnames_offset_x = 2, colnames_offset_y = 0 ,low='red', high='lightgreen')
  open_tree(t1, 10) %>% rotate_tree(10)
  ggsave(pfilout, dpi=300, height = 20, width = 20)
  
  
  pfilout = paste(prout, "dendo_enrich.png", sep = "")
  
  t1 <- ggtree(tupgma2, layout="circular", size=0.8)
  t1 <- t1 %<+% daff + geom_text(aes(color=Enrichment, label = label,  angle=angle,  fontface="bold"), hjust=-0.5, size=1.2) +
    geom_tippoint(aes(color=Enrichment), alpha=0.75, size=1)+
    scale_color_continuous(low='red', high='lightgreen') +
    theme(plot.margin = unit(c(0,0,0,0), "cm")) + 
    theme(legend.position="right") +
    geom_treescale(x = -10, y = 0.1, width = 2, offset = 0,
                   color = "white", linesize = 1E-100, fontsize = 1E-100)
  
  open_tree(t1, 10) %>% rotate_tree(10)
  ggsave(pfilout, dpi=300, height = 15, width = 15)
  
  
  
  
  pfilout = paste(prout, "dendo_prob.png", sep = "")
  
  t1 <- ggtree(tupgma2, layout="circular", size=0.8)
  t1 <- t1 %<+% daff + geom_text(aes(color=prob, label = label,  angle=angle,  fontface="bold"), hjust=-0.5, size=1.2) +
    geom_tippoint(aes(color=prob), alpha=0.75, size=1)+
    scale_color_continuous(low='red', high='lightgreen') +
    theme(plot.margin = unit(c(0,0,0,0), "cm")) + 
    theme(legend.position="right") +
    geom_treescale(x = -10, y = 0.1, width = 2, offset = 0,
                   color = "white", linesize = 1E-100, fontsize = 1E-100)
  
  open_tree(t1, 10) %>% rotate_tree(10)
  ggsave(pfilout, dpi=300, height = 15, width = 15)
  
  
}



dendogramInterfer = function(ddes, dinter, nameDendo, prout){
  
  dinter = as.data.frame(dinter)  
  
  matTrans1 <- scale(ddes)
  d <- dist(matTrans1, method = "euc")
  tupgma2 <- upgma(d, method="ward.D2")
  
  pfilout = paste(prout, "dendo_", nameDendo, ".png", sep = "")
  t1 <- ggtree(tupgma2, layout="circular", size=0.8)
  t1 <- t1 %<+% dinter + geom_text(aes(color=Type, label = label, angle=angle,  fontface="bold"), hjust=-1, size=1.2) +
    geom_tippoint(aes(color=Type), alpha=0.75, size=1)
    #scale_colour_gradientn(colours=rainbow(20))+
  #  theme(plot.margin = unit(c(0,0,0,0), "cm")) + 
  #  geom_treescale(x = 40, y = 40, width = NULL, offset = NULL,
  #                 color = "white", linesize = 1E-100, fontsize = 1E-100)
  
  open_tree(t1, 15) %>% rotate_tree(15)
  ggsave(pfilout, dpi=300, height = 11, width = 11)
  
  
  pfilout = paste(prout, "dendoZscore_", nameDendo, ".png", sep = "")
  t1 <- ggtree(tupgma2, layout="circular", size=0.8)
  t1 <- t1 %<+% dinter +
  #scale_colour_gradientn(colours=rainbow(20))+
    scale_color_continuous(low='red', high='lightgreen') +
    theme(legend.position="right")+
    theme(plot.margin = unit(c(0,0,0,0), "cm")) 
    #geom_treescale(x = 40, y = 40, width = NULL, offset = NULL,
    #               color = "white", linesize = 1E-100, fontsize = 1E-100)
  
  open_tree(t1, 15) %>% rotate_tree(15)
  ggsave(pfilout, dpi=300, height = 11, width = 11)
  
  
}



dendoAtive = function(ddesc, daff, methDist, methAgg, prout, logAff = 0){
  
  
  if(logAff == 0){
    
    if(dim(daff)[2] == 1){
      dmean = daff
      dmean = cbind(rownames(dmean), dmean)
      dmean = as.data.frame(dmean)
      colnames(dmean) = c("ID", "Maff")
      dmean = dmean[rownames(ddesc),]
    }else{
      Maff = rowMeans(daff, na.rm = TRUE)  
      #calibrate for color
      Maff[which(Maff == max(Maff, na.rm = TRUE))] = max(daff, na.rm = TRUE)
      Maff[which(Maff == min(Maff, na.rm = TRUE))] = min(daff, na.rm = TRUE)
      Maff[which(Maff == "NaN")] = NA
      dmean = as.data.frame(Maff)
      dmean = cbind(rownames(dmean), dmean)
      dmean = as.data.frame(dmean)
    }
    
    
    matTrans1 <- scale(ddesc)
    d <- dist(matTrans1, method = methDist)
    tupgma2 <- upgma(d, method=methAgg)
    
    t1 <- ggtree(tupgma2, layout="circular", size=0.8)
    t1 <- t1 %<+% dmean +  geom_tippoint(aes(color=Maff), alpha=0.75, size=3)+
      scale_color_continuous(low='red', high='lightgreen')+
      theme(plot.margin = unit(c(0,0,0,0), "cm")) + 
      geom_treescale(x = 40, y = 40, width = NULL, offset = NULL,
                     color = "white", linesize = 1E-100, fontsize = 1E-100)
    daffplot = as.data.frame(daff)
    t2 = gheatmap(t1, daffplot, font.size = 0, offset = 9, width = 0.2, colnames_offset_x = 2, colnames_offset_y = 0, low = "red", high = "lightgreen")
    open_tree(t2, 15) %>% rotate_tree(15)
    ggsave(paste(prout, "_allactive.png", sep = ""), dpi=300, height = 11, width = 11)
    
    
    # with text inside
    t1 <- ggtree(tupgma2, layout="circular", size=0.8)
    t1 <- t1 %<+% dmean +  geom_tippoint(aes(color=Maff), alpha=0.75, size=3)+
      geom_text(aes(color=Maff, label = label, angle=angle,  fontface="bold"), hjust=-1, size=0.7) +
      scale_color_continuous(low='red', high='lightgreen')+
      theme(plot.margin = unit(c(0,0,0,0), "cm")) + 
      geom_treescale(x = 40, y = 40, width = NULL, offset = NULL,
                     color = "white", linesize = 1E-100, fontsize = 1E-100)
    daffplot = as.data.frame(daff)
    t2 = gheatmap(t1, daffplot, font.size = 0, offset = 9, width = 0.2, colnames_offset_x = 2, colnames_offset_y = 0, low = "red", high = "lightgreen")
    open_tree(t2, 15) %>% rotate_tree(15)
    ggsave(paste(prout, "_allactive_text.png", sep = ""), dpi=300, height = 11, width = 11)
    
    # reduce the matrix of descriptor to consider only active chemicals #
    #####################################################################
    dmean = na.omit(dmean)
    lID = rownames(dmean)
    ddesc = ddesc[lID,]
    
    matTrans1 <- scale(ddesc)
    d <- dist(matTrans1, method = methDist)
    tupgma2 <- upgma(d, method=methAgg)
    
    t1 <- ggtree(tupgma2, layout="circular", size=0.8)
    t1 <- t1 %<+% dmean +  geom_tippoint(aes(color=Maff), alpha=0.75, size=3)+
      scale_color_continuous(low='red', high='lightgreen')+
      theme(plot.margin = unit(c(0,0,0,0), "cm")) + 
      geom_treescale(x = 40, y = 40, width = NULL, offset = NULL,
                     color = "white", linesize = 1E-100, fontsize = 1E-100)
    daffplot = as.data.frame(daff)
    t2 = gheatmap(t1, daffplot, font.size = 0, offset = 9, width = 0.2, colnames_offset_x = 2, colnames_offset_y = 0, low = "red", high = "lightgreen")
    open_tree(t2, 15) %>% rotate_tree(15)
    ggsave(paste(prout, "_spe.png", sep = ""), dpi=300, height = 11, width = 11)
    
    
    # with text inside
    t1 <- ggtree(tupgma2, layout="circular", size=0.8)
    t1 <- t1 %<+% dmean +  geom_tippoint(aes(color=Maff), alpha=0.75, size=3)+
      geom_text(aes(color=Maff, label = label, angle=angle,  fontface="bold"), hjust=-1, size=0.7) +
      scale_color_continuous(low='red', high='lightgreen')+
      theme(plot.margin = unit(c(0,0,0,0), "cm")) + 
      geom_treescale(x = 40, y = 40, width = NULL, offset = NULL,
                     color = "white", linesize = 1E-100, fontsize = 1E-100)
    daffplot = as.data.frame(daff)
    t2 = gheatmap(t1, daffplot, font.size = 0, offset = 9, width = 0.2, colnames_offset_x = 2, colnames_offset_y = 0, low = "red", high = "lightgreen")
    open_tree(t2, 15) %>% rotate_tree(15)
    ggsave(paste(prout, "_spe_text.png", sep = ""), dpi=300, height = 11, width = 11)
  
    
  }else{
    
    daff = log10(daff)
    
    if(dim(daff)[2] == 1){
      dmean = daff
      dmean = cbind(rownames(dmean), dmean)
      dmean = as.data.frame(dmean)
      colnames(dmean) = c("ID", "Maff")
      dmean = dmean[rownames(ddesc),]
    }else{
      Maff = rowMeans(daff, na.rm = TRUE)  
      #calibrate for color
      Maff[which(Maff == max(Maff, na.rm = TRUE))] = max(daff, na.rm = TRUE)
      Maff[which(Maff == min(Maff, na.rm = TRUE))] = min(daff, na.rm = TRUE)
      Maff[which(Maff == "NaN")] = NA
      dmean = as.data.frame(Maff)
      dmean = cbind(rownames(dmean), dmean)
      dmean = as.data.frame(dmean)
    }
    
    
    matTrans1 <- scale(ddesc)
    d <- dist(matTrans1, method = methDist)
    tupgma2 <- upgma(d, method=methAgg)
    
    t1 <- ggtree(tupgma2, layout="circular", size=0.8)
    t1 <- t1 %<+% dmean +  geom_tippoint(aes(color=Maff), alpha=0.75, size=3)+
      scale_color_continuous(low='red', high='lightgreen')+
      theme(plot.margin = unit(c(0,0,0,0), "cm")) + 
      geom_treescale(x = 40, y = 40, width = NULL, offset = NULL,
                     color = "white", linesize = 1E-100, fontsize = 1E-100)
    daffplot = as.data.frame(daff)
    t2 = gheatmap(t1, daffplot, font.size = 0, offset = 9, width = 0.2, colnames_offset_x = 2, colnames_offset_y = 0, low = "red", high = "lightgreen")
    open_tree(t2, 15) %>% rotate_tree(15)
    ggsave(paste(prout, "_log_allactive.png", sep = ""), dpi=300, height = 11, width = 11)
    
    
    # with text inside
    t1 <- ggtree(tupgma2, layout="circular", size=0.8)
    t1 <- t1 %<+% dmean +  geom_tippoint(aes(color=Maff), alpha=0.75, size=3)+
      geom_text(aes(color=Maff, label = label, angle=angle,  fontface="bold"), hjust=-1, size=0.7) +
      scale_color_continuous(low='red', high='lightgreen')+
      theme(plot.margin = unit(c(0,0,0,0), "cm")) + 
      geom_treescale(x = 40, y = 40, width = NULL, offset = NULL,
                     color = "white", linesize = 1E-100, fontsize = 1E-100)
    daffplot = as.data.frame(daff)
    t2 = gheatmap(t1, daffplot, font.size = 0, offset = 9, width = 0.2, colnames_offset_x = 2, colnames_offset_y = 0, low = "red", high = "lightgreen")
    open_tree(t2, 15) %>% rotate_tree(15)
    ggsave(paste(prout, "_log_allactive_text.png", sep = ""), dpi=300, height = 11, width = 11)
    
    # reduce the matrix of descriptor to consider only active chemicals #
    #####################################################################
    dmean = na.omit(dmean)
    lID = rownames(dmean)
    ddesc = ddesc[lID,]
    
    matTrans1 <- scale(ddesc)
    d <- dist(matTrans1, method = methDist)
    tupgma2 <- upgma(d, method=methAgg)
    
    t1 <- ggtree(tupgma2, layout="circular", size=0.8)
    t1 <- t1 %<+% dmean +  geom_tippoint(aes(color=Maff), alpha=0.75, size=3)+
      scale_color_continuous(low='red', high='lightgreen')+
      theme(plot.margin = unit(c(0,0,0,0), "cm")) + 
      geom_treescale(x = 40, y = 40, width = NULL, offset = NULL,
                     color = "white", linesize = 1E-100, fontsize = 1E-100)
    daffplot = as.data.frame(daff)
    t2 = gheatmap(t1, daffplot, font.size = 0, offset = 9, width = 0.2, colnames_offset_x = 2, colnames_offset_y = 0, low = "red", high = "lightgreen")
    open_tree(t2, 15) %>% rotate_tree(15)
    ggsave(paste(prout, "_log_spe.png", sep = ""), dpi=300, height = 11, width = 11)
    
    
    # with text inside
    t1 <- ggtree(tupgma2, layout="circular", size=0.8)
    t1 <- t1 %<+% dmean +  geom_tippoint(aes(color=Maff), alpha=0.75, size=3)+
      geom_text(aes(color=Maff, label = label, angle=angle,  fontface="bold"), hjust=-1, size=0.7) +
      scale_color_continuous(low='red', high='lightgreen')+
      theme(plot.margin = unit(c(0,0,0,0), "cm")) + 
      geom_treescale(x = 40, y = 40, width = NULL, offset = NULL,
                     color = "white", linesize = 1E-100, fontsize = 1E-100)
    daffplot = as.data.frame(daff)
    t2 = gheatmap(t1, daffplot, font.size = 0, offset = 9, width = 0.2, colnames_offset_x = 2, colnames_offset_y = 0, low = "red", high = "lightgreen")
    open_tree(t2, 15) %>% rotate_tree(15)
    ggsave(paste(prout, "_log_spe_text.png", sep = ""), dpi=300, height = 11, width = 11)
    
  }
  
  
}



dendoAtiveDiss = function(ddis, daff, methAgg, prout){
  
  daff = daff[rownames(ddis),]
  
  Maff = rowMeans(daff, na.rm = TRUE)
  #calibrate for color
  Maff[which(Maff == max(Maff, na.rm = TRUE))] = max(daff, na.rm = TRUE)
  Maff[which(Maff == min(Maff, na.rm = TRUE))] = min(daff, na.rm = TRUE)
  Maff[which(Maff == "NaN")] = NA
  dmean = as.data.frame(Maff)
  dmean = cbind(rownames(dmean), dmean)
  dmean = as.data.frame(dmean)
  
  ddis = as.dist(ddis)
  tupgma2 <- upgma(ddis, method=methAgg)
  
  
  t1 <- ggtree(tupgma2, layout="circular", size=0.8)
  t1 <- t1 %<+% dmean +  geom_tippoint(aes(color=Maff), alpha=0.75, size=3)+
    scale_color_continuous(low='red', high='lightgreen')+
    theme(plot.margin = unit(c(0,0,0,0), "cm")) + 
    geom_treescale(x = 40, y = 40, width = NULL, offset = NULL,
                   color = "white", linesize = 1E-100, fontsize = 1E-100)
  
  daffplot = as.data.frame(daff)
  
  t2 = gheatmap(t1, daffplot, font.size = 0, offset = 9, width = 0.2, colnames_offset_x = 2, colnames_offset_y = 0, low = "red", high = "lightgreen")
  
  open_tree(t2, 15) %>% rotate_tree(15)
  
  ggsave(paste(prout, "dendoAff.png", sep = ""), dpi=300, height = 11, width = 11)
}



