#!/usr/bin/env Rscript
library(ggplot2)
library(svglite)

is.integer0 <- function(x)
{
  is.integer(x) && length(x) == 0L
}


################
#     MAIN     #
################

args <- commandArgs(TRUE)
pCAS = args[1]
#pCAS = "/home/borrela2/interference/spDataAnalysis/tox21-spec-hek293-p1/curveType/all"

dCAS = read.csv(pCAS, sep = "\t", header = TRUE)
p = ggplot(dCAS, aes(x=factor(Curves)))+
  geom_bar(stat="count", width=0.7, fill="steelblue")+
  labs(x="Curve type", y = "Count active")+
  theme(axis.text.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 20, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 25, hjust = 0.5, vjust =0.1))

#print(p)
ggsave(paste(pCAS, ".png", sep = ""), width = 10,height = 10, dpi = 300)

if (regexpr("-luc-", pCAS)[[1]] != -1){
  xmax = 125
  dCAS = dCAS[-which(dCAS$Curves == "-2.3"),]
  xleg = "IC50 (uM)"
}else{
  xmax = 80
  xleg = "AC50 (uM)"
}

print(xmax)
#dCAS$Curves = as.factor(dCAS$Curves)

lcurveall = as.factor(c(-2.4, -2.3, -2.2, -2.1, -1.4, -1.2, -1.1, 1.1, 1.2, 1.3, 1.4, 2.1, 2.2, 2.3, 2.4))
lcurve = unique(dCAS$Curves)
for (curvetype in lcurveall){
  print (which(lcurve == curvetype))
  if (is.integer0(which(lcurve == curvetype)) == TRUE){
    add = c("XXXX", curvetype, 80)
    names(add) = c("CASID", "Curves", "Aff")
    
    #add = as.matrix(add)
    print (add)
    dCAS = rbind(dCAS, t(add))
    dCAS = rbind(dCAS, t(add))
    
  }
  
}

dCAS$Curves = as.factor(dCAS$Curves)
dCAS$Aff = as.double(dCAS$Aff)

p = ggplot(dCAS, aes(x=Aff, fill=Curves))+
  geom_density(alpha = 0.8, position = "stack") + 
  theme(text = element_text(size=19))+
  xlim(0,xmax)+
  labs(title="",x=xleg, y = "Density")
ggsave(paste(pCAS, "_density.svg", sep = ""), width = 8, height = 7, dpi = 300, bg="transparent")


# for the log10
dCAS$Aff = log10(dCAS$Aff)
if (xleg == "IC50"){
  xleg = "log10(IC50) (uM)"
}else{
  xleg = "log10(AC50) (uM)"
}

p = ggplot(dCAS, aes(x=Aff, fill=Curves))+
  geom_density(alpha = 0.8, position = "stack") + 
  theme(text = element_text(size=19))+
  xlim(-1.5,3)+
  labs(title="",x=xleg, y = "Density")
ggsave(paste(pCAS, "_log10_density.svg", sep = ""), width = 8, height = 7, dpi = 300, bg="transparent")
