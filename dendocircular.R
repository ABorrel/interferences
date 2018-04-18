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
  minMatrix = min(daff)
  maxMatrix = max(daff)
  
  for(i in seq(1, dim(daff)[2])){
    daff[which(daff[,i] == max(daff[,i])), i] = maxMatrix
    daff[which(daff[,i] == min(daff[,i])), i] = minMatrix
  }
  #daff = daff[,c("Escherichia.coli", "Pseudomonas.aeruginosa",  "Staphylococcus.aureus" , "Streptococcus.pneumoniae")]
  
  daff = as.data.frame(daff)
  daff = cbind(rownames(daff), daff)
  dcluster[,2] = as.character(dcluster[,2])
  
  print(daff)
  
  matTrans1 <- scale(ddes)
  d <- dist(matTrans1, method = "euc")
  tupgma2 <- upgma(d, method="ward.D2")
  tupgma2 <- groupOTU(tupgma2, max(dcluster[,2]))

  pfilout = paste(prout, "dendo_cluster_name.png", sep = "")
    
  t1 <- ggtree(tupgma2, layout="circular", size=0.8, aes(color=cluster))
  t1 <- t1 %<+% dcluster + geom_text(aes(color=cluster, label = cluster, angle=angle,  fontface="bold"), hjust=-1.5, size=1.2) +
    geom_tippoint(aes(color=cluster), alpha=0.75, size=1)+
    #theme(legend.position="right")+
    theme(plot.margin = unit(c(0,0,0,0), "cm")) + 
    geom_treescale(x = 40, y = 40, width = NULL, offset = NULL,
                   color = "white", linesize = 1E-100, fontsize = 1E-100)
  
  
  daff = daff[,-1]
  t2 = gheatmap(t1, daff, font.size = 2, offset = 3, width = 0.5, colnames_offset_x = 2, colnames_offset_y = 0, low = "red", high = "lightgreen") +
    #scale_color_continuous(low='red', high='lightgreen') +
    theme_tree()
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




