from os import path

import assayResults
import pathFolder
import analyseDB
import loadDB
import runExternalSoft
import clusteringDB
import QSARModel
import MCS

#########
# MAIN  #
#########

pSDFToxCast = "/home/borrela2/ToxCast_release20151019/DSSTox_ToxCastRelease_20151019.sdf"
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


# merge AC50/IC50 from different assays #
#########################################

pAC50All = assayResults.mergeAssays(cluc, chepg2, chek293)


# plot correlation #
####################
#cluc.cor3assays(chepg2, chek293)


# plot response curves #
########################
#cluc.responseCurves(drawn=1)
#chek293.responseCurves(drawn=1)
#chepg2.responseCurves(drawn=1)

# cross curve by color and type of assays #
############################################

prCrossCurve = prresults + "crossCurvesResponse/"
pathFolder.createFolder(prCrossCurve)
#chepg2.crossResponseCurves(chek293, prCrossCurve)

# barplot curve type #
######################
prbarplot = chepg2.proutSP + "curveType/"
pathFolder.createFolder(prbarplot)
#chepg2.barplotCurveClass(prbarplot)

prbarplot = cluc.proutSP + "curveType/"
pathFolder.createFolder(prbarplot)
#cluc.barplotCurveClass(prbarplot)

prbarplot = chek293.proutSP + "curveType/"
pathFolder.createFolder(prbarplot)
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


# set FP for different parameters #
###################################

# FP type
lAllMetrics = ['Tanimoto', "Dice", "Cosine", "Sokal", "Russel", "RogotGoldberg", "AllBit", "Kulczynski",
                       "McConnaughey", "Asymmetric", "BraunBlanquet"]
AffFP = ["FPMol", "FPMACCS", "FPpairs", "FPTorsion", "FPMorgan"]

prFP = prMain + "FP/"
pathFolder.createFolder(prFP)
#cDesc.computeFPMatrix(prFP, "FPMol", 'Sokal')
#cDesc.computeFPMatrix(prFP, "FPMol", 'Tanimoto')
#cDesc.computeFPMatrix(prFP, "FPMACCS", 'Tanimoto')
cDesc.computeFPMatrix(prFP, "FPMACCS", 'Dice')
#cDesc.computeFPMatrix(prFP, "FPpairs", 'Dice')
#cDesc.computeFPMatrix(prFP, "FPTorsion", 'Dice')
#cDesc.computeFPMatrix(prFP, "FPpairs", 'Dice')
#cDesc.computeFPMatrix(prFP, "FPMorgan", 'Dice')
#############
###  MCS  ###
#############
# define matrix MCS
#prMCS = pathFolder.createFolder(prMain + "MCS/")
#mcs = MCS.MCSMatrix(cDesc.prSMIclean, prMCS)
#mcs.computeMatrixMCS()

#####################
####  ANALYSIS   ####
#####################

#############################
# set seetings for analysis

corval = 0.9
maxQuantile = 90
splitratio = 0.15
nbCV = 10

######################################
#    clustering for descriptor       #
##### analysis for all Chemical ######
# Double clustering on the main data #
######################################
distMeth = "euclidean"
aggregtype = "ward.D2"#"ward.D2", "complete", "single", "average"
clusterType = "hclust"#"hclust", "kmeans"
optimalCluster = "silhouette"#"gap_stat", "silhouette", "wss", "gap_stat"
doubleclustering = 0


prcluster = prresults + "Descclusters/"
pathFolder.createFolder(prcluster)
cclust = clusteringDB.clustering(cDesc, prcluster, corval, maxQuantile, distmeth=distMeth, aggregtype=aggregtype, clusterType=clusterType, optimalCluster=optimalCluster)
#cclust.createMainClustering(doublecluster=doubleclustering)
#cclust.enrichmentCluster(pAC50All)
#cclust.enrichmentIndex(pAC50All, FP=0)
#cclust.optimalClusteringForEnrich(FP=0)
#cclust.visualizeOptimalClustering(prresults, FP=0)

#####################
# clustering for FP #
#####################

# for dissimilarity matrix frey, mcclain, cindex, silhouette and dunn
prcluster = prresults + "FPclusters/"
optimalCluster = "silhouette"
aggregtype = "ward.D2"# take a look to be sure it is the best aggregation methods
distMeth = None
pathFolder.createFolder(prcluster)
cclust = clusteringDB.clustering(cDesc, prcluster, corval, maxQuantile, distmeth=distMeth, aggregtype=aggregtype, clusterType=clusterType, optimalCluster=optimalCluster)
#cclust.createMainClustering(doublecluster=doubleclustering)
#cclust.enrichmentCluster(pAC50All)
cclust.enrichmentIndex(pAC50All, FP=1)
cclust.optimalClusteringForEnrich(FP=1)
cclust.visualizeOptimalClustering(prresults, FP=1)
www
###############
# MAIN - SOM  #
###############
prSOM = prresults + "SOM/"
pathFolder.createFolder(prSOM)
#cDesc.setConstantPreproc("0", corval, maxQuantile, prSOM)
#pmodelSOM = cDesc.MainSOM(15)

### for luc  ###
################
################
# prep
cluc.writeAC50()
cluc.combineAC50()
pranalysis = cluc.proutSP + "Stat/"
pathFolder.createFolder(pranalysis)
#cDesc.setConstantPreproc(cluc.pAC50, corval, maxQuantile, pranalysis)
#cluc.summarize(pranalysis)

#ranking chemical based on AC50
#cDesc.rankingAC50()
prRank = cluc.proutSP + "ranking/"
pathFolder.createFolder(prRank)
#cluc.rankingTop(100, prPNG, prRank)

prRank = cluc.proutSP + "rankinggood/"
pathFolder.createFolder(prRank)
#cluc.rankingTop(100, prPNG, prRank, 1)


# for HEPG #
############
############

# prep
chepg2.writeAC50()
pranalysis = chepg2.proutSP + "Stat/"
pathFolder.createFolder(pranalysis)
cDesc.setConstantPreproc(chepg2.pAC50, corval, maxQuantile, pranalysis)
#chepg2.summarize(pranalysis)
#prVenn = pathFolder.createFolder(pranalysis + "Venn/")
#chepg2.drawVennPlot(prVenn, prPNG)

# cor with different AC50 available
###################################
#chepg2.corAC50()
# rank AC50
prRank = chepg2.proutSP + "ranking/"
pathFolder.createFolder(prRank)
#chepg2.rankingTop(100, prPNG, prRank)

prRank = chepg2.proutSP + "rankinggood/"
pathFolder.createFolder(prRank)
#chepg2.rankingTop(100, prPNG, prRank, 1)
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

# prep and preliminary analysis
chek293.writeAC50()
#chek293.corAC50()
pranalysis = chek293.proutSP + "Stat/"
pathFolder.createFolder(pranalysis)
#cDesc.setConstantPreproc(chek293.pAC50, corval, maxQuantile, pranalysis)
#chek293.summarize(pranalysis)
#prVeen = pathFolder.createFolder(pranalysis + "Venn/")
#chek293.drawVennPlot(prVeen, prPNG)

#rank by AC50
cDesc.rankingAC50()
prRank = chek293.proutSP + "ranking/"
pathFolder.createFolder(prRank)
#chek293.rankingTop(100, prPNG, prRank)

prRank = chek293.proutSP + "rankinggood/"
pathFolder.createFolder(prRank)
#chek293.rankingTop(100, prPNG, prRank, 1)

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
#clusteringDB.createSOM(cDesc.pdesc1D2D, cluc.pAC50, corval, maxQuantile, pmodelSOM, prSOM)

#### for chepg2  ####
#####################
prSOM = chepg2.proutSP + "SOM/"
pathFolder.createFolder(prSOM)
#clusteringDB.createSOM(cDesc.pdesc1D2D, chepg2.pAC50, corval, maxQuantile, pmodelSOM, prSOM)

#### for HEK293  ####
#####################
prSOM = chek293.proutSP + "SOM/"
pathFolder.createFolder(prSOM)
#clusteringDB.createSOM(cDesc.pdesc1D2D, chek293.pAC50, corval, maxQuantile, pmodelSOM, prSOM)


####################
#   PCA analysis   #
####################

#### for chepg2  ####
#####################
prPCA = chepg2.proutSP + "PCA/"
pathFolder.createFolder(prPCA)
#chepg2.createPCA(cDesc.pdesc1D2D, chepg2.pAC50, corval, maxQuantile, prPCA)

#### for chek293 ####
#####################
prPCA = chek293.proutSP + "PCA/"
pathFolder.createFolder(prPCA)
#chek293.createPCA(cDesc.pdesc1D2D, chek293.pAC50, corval, maxQuantile, prPCA)



##################
# MDS analysis   #
##################

#### for chepg2  ####
#####################
prMDS = chepg2.proutSP + "MDS/"
pathFolder.createFolder(prMDS)
#chepg2.createMDS(cDesc.pdesc1D2D, chepg2.pAC50, corval, maxQuantile, prMDS)

#### for chek293 ####
#####################
prMDS = chek293.proutSP + "MDS/"
pathFolder.createFolder(prMDS)
#chek293.createMDS(cDesc.pdesc1D2D, chek293.pAC50, corval, maxQuantile, prMDS)


##### QSAR modeling ######
##########################

# for luc #
###########
# ----> Regression
typeQSAR = "Reg"
prQSAR = cluc.proutSP + "QSARReg/"
pathFolder.createFolder(prQSAR)

#cluc.combineAC50()
#cModelluc = QSARModel.Model(cDesc.pdesc1D2D, cluc.pAC50, typeQSAR, corval, maxQuantile, splitratio, nbCV, prQSAR)
#cModelluc.prepData()
#cModelluc.buildQSARReg()


# ----> Classification
typeQSAR = "class"
prQSAR = cluc.proutSP + "QSARClass/"
pathFolder.createFolder(prQSAR)

#cluc.combineAC50()
#cModelluc = QSARModel.Model(cDesc.pdesc1D2D, cluc.pAC50, typeQSAR, corval, maxQuantile, splitratio, nbCV, prQSAR)
#cModelluc.prepData()
#cModelluc.buildQSARClass()

# for HEPG2 #
#############
# --> classification
typeQSAR = "class"
prQSARClass = chepg2.proutSP + "QSARclass/"
pathFolder.createFolder(prQSARClass)

#cModelHEPG2 = QSARModel.Model(cDesc.pdesc1D2D, chepg2.pAC50, typeQSAR, corval, maxQuantile, splitratio, nbCV, prQSARClass)
#cModelHEPG2.prepData()
#cModelHEPG2.buildQSARClass()


# ---> regression
typeQSAR = "Reg"
prQSARReg = pathFolder.createFolder(chepg2.proutSP + "QSARreg/")

#cModelHEPG2 = QSARModel.Model(cDesc.pdesc1D2D, chepg2.pAC50, typeQSAR, corval, maxQuantile, splitratio, nbCV, prQSARClass)
#cModelHEPG2.prepData()
#cModelHEPG2.buildQSARClass()


###################################
## cross clustering for all data ##
###################################

# cross color red/green/blue with luc

#clust.corelAllAssays(cluc, chepg2, chek293)

##########################
# Venn diagram by color  #
##########################

prCrossVenn = pathFolder.createFolder(prresults + "CrossVenn/")
#analyseDB.VennCross(cluc, chepg2, chek293, prPNG, prCrossVenn)


#################
#  cross PCA    #
#################
prCrossPCA = pathFolder.createFolder(prresults + "CrossPCA/")
#analyseDB.PCACross(cDesc.pdesc1D2D, chepg2.pAC50, chek293.pAC50, corval, maxQuantile, prCrossPCA)
