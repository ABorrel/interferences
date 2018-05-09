from os import path

import assayResults
import pathFolder
import analyseDB
import loadDB
import runExternalSoft
import clusteringDB
import QSARModel


#########
# MAIN  #
#########


pSDFToxCast= "/home/borrela2/ToxCast_release20151019/DSSTox_ToxCastRelease_20151019.sdf"
pSDFTox21 = "/home/borrela2/Tox21/TOX21SL.sdf"

pluc = "/home/borrela2/interference/data/luc/tox21-luc-biochem-p1/tox21-luc-biochem-p1.txt"
phek293 = "/home/borrela2/interference/data/luc/tox21-spec-hek293-p1/tox21-spec-hek293-p1.txt"
phepg2 = "/home/borrela2/interference/data/luc/tox21-spec-hepg2-p1/tox21-spec-hepg2-p1.txt"

prMain = "/home/borrela2/interference/"
prresults = "/home/borrela2/interference/spDataAnalysis/"

# log folders
prlog = "/home/borrela2/interference/spDataAnalysis/log/"
pathFolder.createFolder(prlog)

# load assays
cluc = assayResults.assays(pluc, prresults, prlog)
chepg2 = assayResults.assays(phepg2, prresults, prlog)
chek293 = assayResults.assays(phek293, prresults, prlog)

# plot correlation # -> not used
####################
#cluc.cor3assays(chepg2, chek293)


# plot response curves #
########################
#cluc.responseCurves(drawn=1)
#chek293.responseCurves(drawn=1)
#chepg2.responseCurves(drawn=1)


# barplot curve type #
######################
#prbarplot = chepg2.proutSP + "curveType/"
#pathFolder.createFolder(prbarplot)
#chepg2.barplotCurveClass(prbarplot)

#prbarplot = cluc.proutSP + "curveType/"
#pathFolder.createFolder(prbarplot)
#cluc.barplotCurveClass(prbarplot)

#prbarplot = chek293.proutSP + "curveType/"
#pathFolder.createFolder(prbarplot)
#chek293.barplotCurveClass(prbarplot)


# IC50 hist #
#############
#cluc.AC50Distribution()
#chek293.AC50Distribution()
#chepg2.AC50Distribution()



# MOLECULAR DESCRIPTOR #
########################

# Create folder with smi and sdf files
#cluc.extractChemical(pSDFTox21)
#chek293.extractChemical(pSDFTox21)
#chepg2.extractChemical(pSDFTox21)


prSMI = prMain + "SMI/"
# compute descriptors

prDesc = prMain + "Desc/"
pathFolder.createFolder(prDesc)

prlogDesc = prMain + "log/"
pathFolder.createFolder(prlogDesc)

prPNG = prMain + "PNG/"
pathFolder.createFolder(prPNG)

cDesc = analyseDB.Descriptors(prSMI, prDesc, prPNG, prresults, prlogDesc)
cDesc.computeDesc()
#cDesc.generatePNG()


#####################
####  ANALYSIS   ####
#####################

#############################
# set seetings for analysis

corval = 0.9
maxQuantile = 90
splitratio = 0.15
nbCV = 10


##### analysis for all Chemical ######
# Double clustering on the main data #
######################################
disttype = "euc"
aggregtype = "ward.D2"#"ward.D2", "complete", "single", "average"
clusterType = "hclust"#"hclust", "kmeans"
optimalCluster = "gap_stat"#"silhouette", "wss", "gap_stat"

#prcluster = prresults + "clusters/"
#pathFolder.createFolder(prcluster)
#cclust = clusteringDB.clustering(cDesc.pdesc1D2D, prcluster, cDesc.prPNG, corval, maxQuantile)
#cclust.createMainClustering(disttype, aggregtype, clusterType, optimalCluster)



### for luc  ###
################
################
# prep
cluc.writeAC50()
cluc.combineAC50()
pranalysis = cluc.proutSP + "Stat/"
pathFolder.createFolder(pranalysis)
cDesc.setConstantPreproc(cluc.pAC50, corval, maxQuantile, pranalysis)

#ranking chemical based on AC50
#cDesc.rankingAC50()
prRank = cluc.proutSP + "ranking/"
pathFolder.createFolder(prRank)
cluc.rankingTop(100, prPNG, prRank)

prRank = cluc.proutSP + "rankinggood/"
pathFolder.createFolder(prRank)
cluc.rankingTop(100, prPNG, prRank, 1)

# clustering -> independant #
#############################
disttype = "euc"
aggregtype = "ward.D2"#"ward.D2", "complete", "single", "average"
clusterType = "hclust"#"hclust", "kmeans"
optimalCluster = "wss"#"silhouette", "wss", "gap_stat"
#cDesc.clustering(disttype, aggregtype, clusterType, optimalCluster)

optimalCluster = "silhouette"
#cDesc.clustering(disttype, aggregtype, clusterType, optimalCluster)

optimalCluster = "gap_stat"
#cDesc.clustering(disttype, aggregtype, clusterType, optimalCluster)


# QSAR modeling #
#################
prQSAR = cluc.proutSP + "QSAR/"
pathFolder.createFolder(prQSAR)

#cModelluc = QSARModel.Model(cDesc.pdesc1D2D, cluc.pAC50, corval, maxQuantile, splitratio, nbCV, prQSAR)
#cModelluc.prepData()
#cModelluc.buildQSARReg()


# for HEPG #
############
############

# prep
chepg2.writeAC50()
pranalysis = chepg2.proutSP + "Stat/"
pathFolder.createFolder(pranalysis)
cDesc.setConstantPreproc(chepg2.pAC50, corval, maxQuantile, pranalysis)

# cor with different AC50 available
###################################
#chepg2.corAC50()
# rank AC50
prRank = chepg2.proutSP + "ranking/"
pathFolder.createFolder(prRank)
chepg2.rankingTop(100, prPNG, prRank)

prRank = chepg2.proutSP + "rankinggood/"
pathFolder.createFolder(prRank)
chepg2.rankingTop(100, prPNG, prRank, 1)
#cDesc.rankingAC50()

# clustering
disttype = "euc"
aggregtype = "ward.D2"#"ward.D2", "complete", "single", "average"
clusterType = "hclust"#"hclust", "kmeans"
optimalCluster = "wss"#"silhouette", "wss", "gap_stat"
#cDesc.clustering(disttype, aggregtype, clusterType, optimalCluster)

optimalCluster = "silhouette"
#cDesc.clustering(disttype, aggregtype, clusterType, optimalCluster)

optimalCluster = "gap_stat"
#cDesc.clustering(disttype, aggregtype, clusterType, optimalCluster)




# for HEK #
############

# prep
chek293.writeAC50()
chek293.corAC50()
pranalysis = chek293.proutSP + "Stat/"
pathFolder.createFolder(pranalysis)
cDesc.setConstantPreproc(chek293.pAC50, corval, maxQuantile, pranalysis)

#rank by AC50
#cDesc.rankingAC50()
prRank = chek293.proutSP + "ranking/"
pathFolder.createFolder(prRank)
chek293.rankingTop(100, prPNG, prRank)

prRank = chek293.proutSP + "rankinggood/"
pathFolder.createFolder(prRank)
chek293.rankingTop(100, prPNG, prRank, 1)

# clustering
disttype = "euc"
aggregtype = "ward.D2"#"ward.D2", "complete", "single", "average"
clusterType = "hclust"#"hclust", "kmeans"
optimalCluster = "wss"#"silhouette", "wss", "gap_stat"
#cDesc.clustering(disttype, aggregtype, clusterType, optimalCluster)

optimalCluster = "silhouette"
#cDesc.clustering(disttype, aggregtype, clusterType, optimalCluster)

optimalCluster = "gap_stat"
#cDesc.clustering(disttype, aggregtype, clusterType, optimalCluster)




# apply main cluster #
######################

#cclust.applyMainClusters(cluc.pAC50, cluc.proutSP)
#cclust.applyMainClusters(chek293.pAC50, chek293.proutSP)
#cclust.applyMainClusters(chepg2.pAC50, chepg2.proutSP)



#################
####   SOM   ####
#################
#### for luc  ####
##################
prSOM = cluc.proutSP + "SOM/"
pathFolder.createFolder(prSOM)
clusteringDB.createSOM(cDesc.pdesc1D2D, cluc.pAC50, corval, maxQuantile, prSOM)

#### for chepg2  ####
#####################
prSOM = chepg2.proutSP + "SOM/"
pathFolder.createFolder(prSOM)
#clusteringDB.createSOM(cDesc.pdesc1D2D, chepg2.pAC50, corval, maxQuantile, prSOM)


#### for HEK293  ####
#####################
prSOM = chek293.proutSP + "SOM/"
pathFolder.createFolder(prSOM)
#clusteringDB.createSOM(cDesc.pdesc1D2D, chek293.pAC50, corval, maxQuantile, prSOM)



###################################
## cross clustering for all data ##
###################################

# cross color red/green/blue with luc

#clust.corelAllAssays(cluc, chepg2, chek293)