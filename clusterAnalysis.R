#!/usr/bin/env Rscript
source ("tool.R")
source("dendocircular.R")
source("radialPlot.R")


################
#     MAIN     #
################

args <- commandArgs(TRUE)
pdesc = args[1]
paff = args[2] #to take affinity
pcluster = args[3]
prout = args[4]
valcor = as.double(args[5])
maxQuantile = as.double(args[6])
logaff = as.integer(args[7])


#pdesc = "/home/aborrel/fluoroquinolones/results/desc/desc_compound.csv"
#paff = "/home/aborrel/fluoroquinolones/MIC_currated.csv"
#pcluster = "/home/aborrel/fluoroquinolones/results/desc_analysis/0.8/Table_hclust_ward.D2_gap_stat.csv"
#prout = "/home/aborrel/fluoroquinolones/results/CrossClusterAnalysis/hclust_ward.D2_gap_stat/"

#valcor = 0.80
#maxQuantile = 85
#logaff = 1


##############################
# Process descriptors matrix #
##############################
dglobal = openData(pdesc, valcor, prout, c(1,2))
dglobal = dglobal[[1]]

rownames(dglobal) = dglobal[,1]
dglobal = dglobal[,-1]
dglobal = dglobal[,-1]

# remove descriptor with a distribution on one quantile
dglobal = delnohomogeniousdistribution(dglobal, maxQuantile)

###################
# Affinity Matrix #
###################
# Opening
daffinity = read.csv(paff, sep = ",", header = TRUE)
rownames(daffinity) = daffinity[,1]
daffinity = daffinity[,-1]
print(dim(daffinity))

# transform
if(logaff == 1){
  daffinity = -log10(daffinity)
}

# merge with data descriptors
lID = intersect(rownames(daffinity), rownames(dglobal))

dglobal = dglobal[lID,]
daffinity = daffinity[lID,]


##################
# cluster matrix #
##################

dcluster = read.csv(pcluster)
rownames(dcluster) = dcluster[,1]
dcluster = dcluster[lID,]


print (dim(dcluster))
print (dim(dglobal))
print (dim(daffinity))

# dendogram cluster #
#####################

#dendogramCluster(dglobal, daffinity, dcluster, prout)


# radial plot by cluster #
##########################

lcluster = unique(dcluster[,2])
daffinity = daffinity[,c("Escherichia.coli", "Pseudomonas.aeruginosa",  "Staphylococcus.aureus" , "Streptococcus.pneumoniae")]
for(cluster in lcluster){
  dtemp = daffinity[which(dcluster[,2] == cluster),]
  radialByCluster(dtemp, paste(prout, cluster, ".svg", sep = ""))
}
