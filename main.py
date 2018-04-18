from idlelib import run
from os import path

import assayResults
import pathFolder
import analyseDB
import loadDB
import runExternalSoft
import clusteringDB


#########
# MAIN  #
#########


pSDFToxCast= "/home/borrela2/ToxCast_release20151019/DSSTox_ToxCastRelease_20151019.sdf"
pSDFTox21 = "/home/borrela2/Tox21/TOX21SL.sdf"

pluc = "/home/borrela2/interference/data/luc/tox21-luc-biochem-p1/tox21-luc-biochem-p1.txt"
phek293 = "/home/borrela2/interference/data/luc/tox21-spec-hek293-p1/tox21-spec-hek293-p1.txt"
phepg2 = "/home/borrela2/interference/data/luc/tox21-spec-hepg2-p1/tox21-spec-hepg2-p1.txt"

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
#cluc.responseCurves()
#chek293.responseCurves()
#chepg2.responseCurves()


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


prSMI = prresults + "SMI/"
# compute descriptors

prDesc = prresults + "Desc/"
pathFolder.createFolder(prDesc)

prlogDesc = prDesc + "log/"
pathFolder.createFolder(prlogDesc)

cDesc = analyseDB.Descriptors(prSMI, prDesc, prlogDesc)
cDesc.computeDesc()


#####################
####  ANALYSIS   ####
#####################

#############################
# set seetings for analysis

corval = 0.9
maxQuantile = 90


##### analysis for all Chemical ######
# Double clustering on the main data #
######################################
disttype = "euc"
aggregtype = "ward.D2"#"ward.D2", "complete", "single", "average"
clusterType = "hclust"#"hclust", "kmeans"
optimalCluster = "gap_stat"#"silhouette", "wss", "gap_stat"

prcluster = prresults + "clusters/"
pathFolder.createFolder(prcluster)
cclust = clusteringDB.clustering(cDesc.pdesc1D2D, prcluster, corval, maxQuantile)
cclust.createMainClustering(disttype, aggregtype, clusterType, optimalCluster)

prSOM = prresults + "SOM/"
pathFolder.createFolder(prSOM)
#clusteringDB.createSOM(cDesc.pdesc1D2D, corval, maxQuantile, prSOM)


### for luc  ###
################
################
# prep
cluc.writeAC50()
pranalysis = cluc.proutSP + "Stat/"
pathFolder.createFolder(pranalysis)
out = cDesc.setConstantPreproc(cluc.pAC50, corval, maxQuantile, pranalysis)
#cDesc.rankingAC50()

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


# apply main cluster #
######################

cclust.applyMainClusters(cluc.pAC50, cluc.proutSP)



# for HEPG #
############

# prep
chepg2.writeAC50()
pranalysis = chepg2.proutSP + "Stat/"
pathFolder.createFolder(pranalysis)
out = cDesc.setConstantPreproc(chepg2.pAC50, corval, maxQuantile, pranalysis)
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
pranalysis = chek293.proutSP + "Stat/"
pathFolder.createFolder(pranalysis)
out = cDesc.setConstantPreproc(chek293.pAC50, corval, maxQuantile, pranalysis)
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





###################################
## cross clustering for all data ##
###################################

