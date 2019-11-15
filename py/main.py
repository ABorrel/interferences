from os import path, listdir
from scipy.linalg import lu

import assayResults
import pathFolder
import analyseDB
import loadDB
import runExternalSoft
import clusteringDB
import QSARModel
#import MCS
import cytox

#########
# MAIN  #
#########

prMain = "/home/borrela2/interference2/"
PRTESTING = prMain + "testing/"
#prMain = "c://Users/Aborrel/Desktop/NIEHS-work/interference/"

pSDFToxCast = "/home/borrela2/ToxCast_release20151019/DSSTox_ToxCastRelease_20151019.sdf"
pSDFTox21 = "/home/borrela2/data/Tox21/TOX21SL.sdf"
prCytox = prMain + "data/Judson_2016-cytotox/"

pluc = prMain + "data/luc/tox21-luc-biochem-p1/tox21-luc-biochem-p1.txt"
phek293 = prMain +"data/luc/tox21-spec-hek293-p1/tox21-spec-hek293-p1.txt"
phepg2 = prMain + "data/luc/tox21-spec-hepg2-p1/tox21-spec-hepg2-p1.txt"
prresults = pathFolder.createFolder(prMain + "spDataAnalysis/")

# log folders
prlog = prMain + "spDataAnalysis/log/"
pathFolder.createFolder(prlog)

curvecutoff = 3
effcutoff = 30
# load assays
cluc = assayResults.assays(pluc, curvecutoff, effcutoff, curvePositive=1, curveNegative=1, prcytox=prCytox, prout=prresults, prlog=prlog)
chepg2 = assayResults.assays(phepg2, curvecutoff, effcutoff, curvePositive=1, curveNegative=0, prcytox=prCytox, prout=prresults, prlog=prlog)
chek293 = assayResults.assays(phek293, curvecutoff, effcutoff, curvePositive=1, curveNegative=0, prcytox=prCytox, prout=prresults, prlog=prlog)

# png
prPNG = prMain + "PNG/"
pathFolder.createFolder(prPNG)

# filter cytox => from judson 2016 #
####################################

# for luc
#=> no filter
#cluc.writeAC50(filtercurvefit=0, filterefficacy=0, filterburst=0, combine=0)
#cluc.summarize(pathFolder.createFolder(cluc.proutSP + "Stat/"))
#=> combine
#cluc.writeAC50(filtercurvefit=0, filterefficacy=0, filterburst=0, combine=1)
#cluc.summarize(pathFolder.createFolder(cluc.proutSP + "Stat/"))
#=> curve filter
#cluc.writeAC50(filtercurvefit=1, filterefficacy=0, filterburst=0, combine=0)
#cluc.summarize(pathFolder.createFolder(cluc.proutSP + "Stat/"))
#=> efficacy filter
#cluc.writeAC50(filtercurvefit=1, filterefficacy=1, filterburst=0, combine=0)
#cluc.summarize(pathFolder.createFolder(cluc.proutSP + "Stat/"))
#=> combining
cluc.writeAC50(filtercurvefit=1, filterefficacy=1, filterburst=0, combine=1)
#cluc.summarize(pathFolder.createFolder(cluc.proutSP + "Stat/"))

#======> curve filter and burst // not applied
#cluc.writeAC50(filtercurvefit=1, filterefficacy=1, filterburst=1, combine=1)
#cluc.summarize(pathFolder.createFolder(cluc.proutSP + "Stat/"))

# for hepg2
#=> no filter
#chepg2.writeAC50(filtercurvefit=0, filterefficacy=0, filterburst=0, combine=0)
#chepg2.summarize(pathFolder.createFolder(chepg2.proutSP + "Stat/"))
#=> curve filter
#chepg2.writeAC50(filtercurvefit=1, filterefficacy=0, filterburst=0, combine=0)
#chepg2.summarize(pathFolder.createFolder(chepg2.proutSP + "Stat/"))
#=> efficacy
#chepg2.writeAC50(filtercurvefit=1, filterefficacy=1, filterburst=0, combine=0)
#chepg2.summarize(pathFolder.createFolder(chepg2.proutSP + "Stat/"))
#=> curve filter and burst
chepg2.writeAC50(filtercurvefit=1, filterefficacy=1, filterburst=1, combine=0)
#chepg2.summarize(pathFolder.createFolder(chepg2.proutSP + "Stat/"))

# for hek293
#=> no filter
#chek293.writeAC50(filtercurvefit=0, filterefficacy=0, filterburst=0, combine=0)
#chek293.summarize(pathFolder.createFolder(chek293.proutSP + "Stat/"))
#=> curve filter
#chek293.writeAC50(filtercurvefit=1, filterefficacy=0, filterburst=0, combine=0)
#chek293.summarize(pathFolder.createFolder(chek293.proutSP + "Stat/"))
#=> efficacy
#chek293.writeAC50(filtercurvefit=1, filterefficacy=1, filterburst=0, combine=0)
#chek293.summarize(pathFolder.createFolder(chek293.proutSP + "Stat/"))
#=> curve filter and burst
chek293.writeAC50(filtercurvefit=1, filterefficacy=1, filterburst=1, combine=0)
#chek293.summarize(pathFolder.createFolder(chek293.proutSP + "Stat/"))

###############
# analyse DYE #
###############
#import dyeAnalysis

#pDye = prMain + "data/Tox21_color_chemicals.csv"
#prDye = pathFolder.createFolder(prMain + "DYE_analysis/")
#cDYE = dyeAnalysis.DYE(pDye, [chepg2, chek293], prDye)
#cDYE.crossDyeAssays()


####################
# analyse spectrum #
####################
#import photoChemDB

#pphotochem = prMain + "data/PhotochemCAD3/PCAD3CompdDatabase2018/2018_03PCAD3.db"
#prSpectrum = pathFolder.createFolder(prMain + "SPECTRUM_analysis/")

#cSpectrum = photoChemDB.photoChem(pphotochem, [chepg2, chek293], prSpectrum)
#cSpectrum.crossSpectrumAssays()


# merge AC50/IC50 from different assays #
#########################################

pAC50All = assayResults.mergeAssays(cluc, chepg2, chek293)
#prhist = pathFolder.createFolder(prresults + "hist/")
#assayResults.histogramAC50(pAC50All, prhist)

# plot correlation #
####################
#cluc.cor3assays(chepg2, chek293)

# plot response curves #
########################

####### => no tested
#cluc.responseCurves(drawn=1)
#chek293.responseCurves(drawn=1)
#chepg2.responseCurves(drawn=1)


# cross curve by color and type of assays #
############################################
#prCrossCurve = prresults + "crossCurvesResponse/"
#pathFolder.createFolder(prCrossCurve)
#chepg2.crossResponseCurves(chek293, prCrossCurve)

#ddd

# extract chemical with type of curve
######################################
#prChemCurve = prresults + "chemCurve/"
#pathFolder.createFolder(prChemCurve)
#assayResults.ChemByCurve(cluc, prPNG, pathFolder.createFolder(prChemCurve + "luc/"))
#assayResults.ChemByCurve(chepg2, prPNG, pathFolder.createFolder(prChemCurve + "HepG2/"))
#assayResults.ChemByCurve(chek293, prPNG, pathFolder.createFolder(prChemCurve + "HEK293/"))

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


# IC50 hist -> by cell #
########################
#cluc.AC50Distribution()
#chek293.AC50Distribution()
#chepg2.AC50Distribution()

#dd
####### => no tested


#############################################=> first round

# MOLECULAR DESCRIPTOR #
########################

prDesc = pathFolder.createFolder(prMain + "Desc/")

# Create folder with smi and sdf files
#cluc.extractChemical(pSDFTox21, prDesc)
#chek293.extractChemical(pSDFTox21, prDesc)
#chepg2.extractChemical(pSDFTox21, prDesc)

prlogDesc = prMain + "log/"
pathFolder.createFolder(prlogDesc)

cDesc = analyseDB.Descriptors(prDesc, prDesc + "SMI_origin/", prresults, prlogDesc)
cDesc.computeAllDesc(compute=0)


# set FP for different parameters #
###################################

# FP type
#lMetrics = ['Tanimoto', "Dice", "Cosine", "Sokal", "Russel", "RogotGoldberg", "AllBit", "Kulczynski",
#                       "McConnaughey", "Asymmetric", "BraunBlanquet"]
#lAffFP = ["Mol", "MACCS", "pairs", "Torsion", "Morgan"]

#prFP = prMain + "FP/"
#pathFolder.createFolder(prFP)
#laggregtype = ["ward.D2", "complete", "single", "average"]
#distMeth = None
#corval = 0.9
#maxQuantile = 90
#clusterType = "hclust"
#optimalCluster = "silhouette"

# run FP have to add in function
#prcluster = prresults + "FPclusters/"
#pathFolder.createFolder(prcluster)


#for Metric in lMetrics:
#    for AffFP in lAffFP:
#        cDesc.computeFPMatrix(prFP, AffFP, Metric)

#        for aggregtype in laggregtype:
#            cclust = clusteringDB.clustering(cDesc, prcluster, corval, maxQuantile, distmeth=distMeth,
#                                             aggregtype=aggregtype, clusterType=clusterType,
#                                             optimalCluster=optimalCluster)
#            cclust.enrichmentIndex(pAC50All, FP=1)
#            cclust.optimalClusteringForEnrich(FP=1)
#            cclust.visualizeOptimalClustering(prresults, FP=1)

#dd
#cDesc.computeFPMatrix(prFP, "Mol", 'Sokal')
#cDesc.computeFPMatrix(prFP, "Mol", 'Tanimoto')
#cDesc.computeFPMatrix(prFP, "MACCS", 'Tanimoto')
#cDesc.computeFPMatrix(prFP, "MACCS", 'Dice')
#cDesc.computeFPMatrix(prFP, "pairs", 'Dice')
#cDesc.computeFPMatrix(prFP, "Torsion", 'Dice')
#cDesc.computeFPMatrix(prFP, "Morgan", 'Dice')

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
#distMeth = "euclidean"
#aggregtype = "ward.D2"#"ward.D2", "complete", "single", "average"
#clusterType = "hclust"#"hclust", "kmeans"
#optimalCluster = "silhouette"#"gap_stat", "silhouette", "wss", "gap_stat"
#doubleclustering = 0


#laggregtype = ["ward.D2", "complete", "single", "average"]
#ldistMeth = ["euclidean", "maximum", "manhattan", "canberra", "binary","minkowski"]
#prcluster = prresults + "Descclusters/"
#pathFolder.createFolder(prcluster)
#for aggregtype in laggregtype:
#    for distMeth in ldistMeth:
#        cclust = clusteringDB.clustering(cDesc, prcluster, corval, maxQuantile, distmeth=distMeth, aggregtype=aggregtype, clusterType=clusterType, optimalCluster=optimalCluster)
        #cclust.createMainClustering(doublecluster=doubleclustering)
        #cclust.enrichmentCluster(pAC50All)
#        cclust.enrichmentIndex(pAC50All, FP=0)
#        cclust.optimalClusteringForEnrich(FP=0)
#        cclust.visualizeOptimalClustering(prresults, FP=0)


#######
# clusterize only active chemicals
######
#-> by descriptors
#nbNA = 50
#prDescAct = pathFolder.createFolder(prresults + "DescActive/")
#cclustAct = clusteringDB.clustering(cDesc, prDescAct, corval, maxQuantile, nbNA=nbNA, distmeth=distMeth, aggregtype=aggregtype, clusterType=clusterType, optimalCluster=optimalCluster)
#cclustAct.clusterActive(pAC50All)
# using FP
# => to do but no used at the moment


#####################
# clustering for FP #
#####################

# for dissimilarity matrix frey, mcclain, cindex, silhouette and dunn
#prcluster = prresults + "FPclusters/"
#optimalCluster = "silhouette"
#laggregtype = ["ward.D2", "complete", "single", "average"]# take a look to be sure it is the best aggregation methods
#distMeth = None
#pathFolder.createFolder(prcluster)
#for aggregtype in laggregtype:
#    cclust = clusteringDB.clustering(cDesc, prcluster, corval, maxQuantile, distmeth=distMeth, aggregtype=aggregtype, clusterType=clusterType, optimalCluster=optimalCluster)
    ####cclust.createMainClustering(doublecluster=doubleclustering)
    ###cclust.enrichmentCluster(pAC50All)
#    cclust.enrichmentIndex(pAC50All, FP=1)
#    cclust.optimalClusteringForEnrich(FP=1)
#    cclust.visualizeOptimalClustering(prresults, FP=1)


#################
# Student Tests #
#################

#prTtest = pathFolder.createFolder(prresults + "Ttest/")
#nbNA = 500
#cDesc.setConstantPreproc("0", corval, maxQuantile, nbNA, prTtest)
#cDesc.Ttest(pAC50All, prTtest)

###############
# MAIN - SOM  #
###############
prSOM = prresults + "SOM/"
pathFolder.createFolder(prSOM)
nbNA = 100
cDesc.setConstantPreproc("0", corval, maxQuantile, nbNA, prSOM)
pmodelSOM = cDesc.MainSOM(15)


# SOM active  #
###############
# for luciferase
#prSOMAct = prresults + "SOMactiveluc/"
#pathFolder.createFolder(prSOMAct)

#nbNA = 100
#cDesc.prepareActiveMatrix(corval, maxQuantile, nbNA, pAC50All, prSOMAct, luciferase=1)
#cDesc.createActiveSOM(15, prSOMAct, pmodelSOM)
#cDesc.extractActivebySOM()

# for autofluorescence
#prSOMAct = prresults + "SOMactivefluo/"
#pathFolder.createFolder(prSOMAct)
#cDesc.prepareActiveMatrix(corval, maxQuantile, nbNA, pAC50All, prSOMAct, luciferase=0)
#cDesc.createActiveSOM(15, prSOMAct, pmodelSOM)
#cDesc.extractActivebySOM()


### for luc  ###
################
################
# prep
#nbNA = 1000
#pranalysis = cluc.proutSP + "Stat/"
#pathFolder.createFolder(pranalysis)
#cDesc.setConstantPreproc(cluc.pAC50, corval, maxQuantile, nbNA, pranalysis)
#cluc.summarize(pranalysis)
##ranking chemical based on AC50
#cDesc.rankingAC50()
#prRank = cluc.proutSP + "ranking/"
#pathFolder.createFolder(prRank)
#cluc.rankingTop(100, prPNG, prRank)
#prRank = cluc.proutSP + "rankinggood/"
#pathFolder.createFolder(prRank)
#cluc.rankingTop(100, prPNG, prRank, 1)

# for HEPG #
############
############
#nbNA = 500
#pranalysis = chepg2.proutSP + "Stat/"
#pathFolder.createFolder(pranalysis)
#cDesc.setConstantPreproc(chepg2.pAC50, corval, maxQuantile, nbNA, pranalysis)

#prVenn = pathFolder.createFolder(pranalysis + "Venn/")
#chepg2.drawVennPlot(prVenn, prPNG)

# cor with different AC50 available
###################################
#chepg2.corAC50()
# rank AC50
#prRank = chepg2.proutSP + "ranking/"
#pathFolder.createFolder(prRank)
#chepg2.rankingTop(100, prPNG, prRank)

#prRank = chepg2.proutSP + "rankinggood/"
#pathFolder.createFolder(prRank)
#chepg2.rankingTop(100, prPNG, prRank, 1)
#cDesc.rankingAC50()

# clustering
#disttype = "euc"
#aggregtype = "ward.D2"#"ward.D2", "complete", "single", "average"
#clusterType = "hclust"#"hclust", "kmeans"
#optimalCluster = "wss"#"silhouette", "wss", "gap_stat"
#cDesc.clustering(disttype, aggregtype, clusterType, optimalCluster)

#optimalCluster = "silhouette"
#cDesc.clustering(disttype, aggregtype, clusterType, optimalCluster)

#optimalCluster = "gap_stat"
#cDesc.clustering(disttype, aggregtype, clusterType, optimalCluster)




# for HEK #
############

# prep and preliminary analysis
#chek293.corAC50()
#pranalysis = chek293.proutSP + "Stat/"
#pathFolder.createFolder(pranalysis)
#cDesc.setConstantPreproc(chek293.pAC50, corval, maxQuantile, nbNA, pranalysis)
#chek293.summarize(pranalysis)
#prVeen = pathFolder.createFolder(pranalysis + "Venn/")
#chek293.drawVennPlot(prVeen, prPNG)

#rank by AC50
#cDesc.rankingAC50()
#prRank = chek293.proutSP + "ranking/"
#pathFolder.createFolder(prRank)
#chek293.rankingTop(100, prPNG, prRank)

#prRank = chek293.proutSP + "rankinggood/"
#pathFolder.createFolder(prRank)
#chek293.rankingTop(100, prPNG, prRank, 1)

# clustering
#disttype = "euc"
#aggregtype = "ward.D2"#"ward.D2", "complete", "single", "average"
#clusterType = "hclust"#"hclust", "kmeans"
#optimalCluster = "wss"#"silhouette", "wss", "gap_stat"
#cDesc.clustering(disttype, aggregtype, clusterType, optimalCluster)

#optimalCluster = "silhouette"
#cDesc.clustering(disttype, aggregtype, clusterType, optimalCluster)

#optimalCluster = "gap_stat"
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
#nbNA = 1000
#prSOM = cluc.proutSP + "SOM/"
#pathFolder.createFolder(prSOM)
#clusteringDB.createSOM(cDesc.pdesc1D2D, cluc.pAC50, corval, maxQuantile, pmodelSOM, nbNA, prSOM)

#### for chepg2  ####
#####################
#nbNA = 500
#prSOM = chepg2.proutSP + "SOM/"
#pathFolder.createFolder(prSOM)
#clusteringDB.createSOM(cDesc.pdesc1D2D, pAC50All, corval, maxQuantile, pmodelSOM, nbNA, prSOM)

#### for HEK293  ####
#####################
#nbNA = 500
#prSOM = chek293.proutSP + "SOM/"
#pathFolder.createFolder(prSOM)
#clusteringDB.createSOM(cDesc.pdesc1D2D, pAC50All, corval, maxQuantile, pmodelSOM, nbNA, prSOM)


####################
#   PCA analysis   #
####################

##### for luc #####
###################
#nbNA = 1000
#prPCA = cluc.proutSP + "PCA/"
#pathFolder.createFolder(prPCA)
#cluc.createPCA(cDesc.pdesc1D2D, cluc.pAC50, corval, maxQuantile, nbNA, prPCA)

#### for chepg2  ####
#####################
#nbNA = 500
#prPCA = chepg2.proutSP + "PCA/"
#pathFolder.createFolder(prPCA)
#chepg2.createPCA(cDesc.pdesc1D2D, chepg2.pAC50, corval, maxQuantile, nbNA, prPCA)

#### for chek293 ####
#####################
#prPCA = chek293.proutSP + "PCA/"
#pathFolder.createFolder(prPCA)
#chek293.createPCA(cDesc.pdesc1D2D, chek293.pAC50, corval, maxQuantile, nbNA, prPCA)



##################
# MDS analysis   #
##################

#### for chepg2  ####
#####################
#prMDS = chepg2.proutSP + "MDS/"
#pathFolder.createFolder(prMDS)
#chepg2.createMDS(cDesc.pdesc1D2D, chepg2.pAC50, corval, maxQuantile, prMDS)

#### for chek293 ####
#####################
#prMDS = chek293.proutSP + "MDS/"
#pathFolder.createFolder(prMDS)
#chek293.createMDS(cDesc.pdesc1D2D, chek293.pAC50, corval, maxQuantile, prMDS)

##### QSAR modeling ######
##########################
ratioAct = 0.3
nbRepeat = 10
nbNA = 500
ltypeCellChannel = ["cell_blue_n", "cell_green_n", "cell_red_n", "med_blue_n", "med_green_n", "med_red_n"]
typeQSAR = "class"

##### Classification  #####
###########################

# for each chanel and cell line #
#################################
typeData = "all"
#QSARModel.runQSARClass(cDesc, cluc, pAC50All, corval, maxQuantile, splitratio, nbCV, ratioAct, nbRepeat, nbNA, "Luc", ["IC50"], typeData, cluc.proutSP + "QSARclass/")
#QSARModel.runQSARClass(cDesc, chepg2, pAC50All, corval, maxQuantile, splitratio, nbCV, ratioAct, nbRepeat, nbNA,"hepg2", ltypeCellChannel, typeData, chepg2.proutSP + "QSARclass/")
QSARModel.runQSARClass(cDesc, chek293, pAC50All, corval, maxQuantile, splitratio, nbCV, ratioAct, nbRepeat, nbNA,"hek293", ltypeCellChannel, typeData, chek293.proutSP + "QSARclass/")
ss

# for each chanel and favorize active chemical #
################################################
#typeData = "active"
#QSARModel.runQSARClass(cDesc, chepg2, pAC50All, corval, maxQuantile, splitratio, nbCV, ratioAct, nbRepeat, nbNA,"hepg2", ltypeCellChannel, typeData, chepg2.proutSP + "QSARclassActive/")
#QSARModel.runQSARClass(cDesc, chek293, pAC50All, corval, maxQuantile, splitratio, nbCV, ratioAct, nbRepeat, nbNA,"hek293", ltypeCellChannel, typeData, chek293.proutSP + "QSARclassActive/")

# for each color #
##################
typeData = "color"
nbNA = 500
QSARModel.runQSARClass(cDesc, chepg2, pAC50All, corval, maxQuantile, splitratio, nbCV, ratioAct, nbRepeat, nbNA,"blue", ltypeCellChannel, typeData, prresults + "QSARclassColor/")
QSARModel.runQSARClass(cDesc, chepg2, pAC50All, corval, maxQuantile, splitratio, nbCV, ratioAct, nbRepeat,nbNA, "green", ltypeCellChannel, typeData, prresults + "QSARclassColor/")
QSARModel.runQSARClass(cDesc, chepg2, pAC50All, corval, maxQuantile, splitratio, nbCV, ratioAct, nbRepeat,nbNA, "red", ltypeCellChannel, typeData, prresults + "QSARclassColor/")
ddd
# for all color #
##################
typeData = "crosscolor"
nbNA = 500
QSARModel.runQSARClass(cDesc, chepg2, pAC50All, corval, maxQuantile, splitratio, nbCV, ratioAct, nbRepeat, nbNA, "all", ltypeCellChannel, typeData, prresults + "QSARclassCrossColor/")



# PUBCHEM classification model #
################################


# for clustering
#for i in range(1, 10):
    # for luc #
    ###########
    # --> classification
#    prQSAR = cluc.proutSP + "QSARclassAct/" + str(i) + "/"
#    pathFolder.createFolder(prQSAR)
#    if len(listdir(prQSAR)) == 0:
#        cluc.combineAC50()
#        cModelluc = QSARModel.Model(cDesc.pdesc1D2D, cluc.pAC50, pAC50All, typeQSAR, corval, maxQuantile, splitratio, nbCV, ratioAct, nbNA,"Luc", ["IC50"], prQSAR)
#        cModelluc.prepData()
#        cModelluc.buildQSARClass()

    # for HEK293 #
    ##############
    # --> classification
#    prQSARClass = chek293.proutSP + "QSARclassAct/" + str(i) + "/"
#    pathFolder.createFolder(prQSARClass)
#    if len(listdir(prQSARClass)) == 0:
#        cModelHEK293 = QSARModel.Model(cDesc.pdesc1D2D, chek293.pAC50, pAC50All, typeQSAR, corval, maxQuantile, splitratio, nbCV, ratioAct, nbNA,chek293.proutSP.split("-")[-2], ltypeCellChannel, prQSARClass)
#        cModelHEK293.prepData()
#        cModelHEK293.buildQSARClass()


    # for HEPG2 #
    #############
    # --> classification
#    prQSARClass = chepg2.proutSP + "QSARclassAct/" + str(i) + "/"
#    pathFolder.createFolder(prQSARClass)
#    if len(listdir(prQSARClass)) == 0:
#        cModelHEPG2 = QSARModel.Model(cDesc.pdesc1D2D, chepg2.pAC50, pAC50All, typeQSAR, corval, maxQuantile, splitratio, nbCV, ratioAct, nbNA,chepg2.proutSP.split("-")[-2], ltypeCellChannel, prQSARClass)
#        cModelHEPG2.prepData()
#        cModelHEPG2.buildQSARClass()





# for luc #
###########
# --> merge
#prQSARAV = pathFolder.createFolder(cluc.proutSP + "QSARclass/Average/")
#QSARModel.mergeResults(cluc.proutSP + "QSARclass/" , prQSARAV)
# for HEPG2#
############
# --> merge
#prQSARAV = pathFolder.createFolder(chepg2.proutSP + "QSARclass/Average/")
#QSARModel.mergeResults(chepg2.proutSP + "QSARclass/" , prQSARAV)

# for HEK293#
#############
# --> merge
#prQSARAV = pathFolder.createFolder(chek293.proutSP + "QSARclass/Average/")
#QSARModel.mergeResults(chek293.proutSP + "QSARclass/" , prQSARAV)


##### Regression  #####
#######################

# luciferase
# ----> Regression
#typeQSAR = "Reg"
#prQSAR = cluc.proutSP + "QSARReg/"
#pathFolder.createFolder(prQSAR)

#cluc.combineAC50()
#cModelluc = QSARModel.Model(cDesc.pdesc1D2D, cluc.pAC50, typeQSAR, corval, maxQuantile, splitratio, nbCV, prQSAR)
#cModelluc.prepData()
#cModelluc.buildQSARReg()

# hepg2
# ---> regression
#typeQSAR = "Reg"
#prQSARReg = pathFolder.createFolder(chepg2.proutSP + "QSARreg/")


#hek293
# ---> regression
#typeQSAR = "Reg"
#prQSARReg = pathFolder.createFolder(chepg2.proutSP + "QSARreg/")

#cModelHEPG2 = QSARModel.Model(cDesc.pdesc1D2D, chepg2.pAC50, typeQSAR, corval, maxQuantile, splitratio, nbCV, nbNA,prQSARClass)
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
#prCrossVenn = pathFolder.createFolder(prresults + "CrossVenn/")
#analyseDB.VennCross(cluc, chepg2, chek293, prPNG, prCrossVenn)



##################################
# Venn diagram with PUBCHEM data #
##################################
#prVennPubmed = pathFolder.createFolder(prMain + "VennTox21Pubmed/")
#analyseDB.VennTox21PubMed(pAC50All, prMain, prVennPubmed)



#################
#  cross PCA    #
#################
#nbNA = 500
#prCrossPCA = pathFolder.createFolder(prresults + "CrossPCA/")
#analyseDB.PCACross(cDesc.pdesc1D2D, chepg2.pAC50, chek293.pAC50, nbNA, corval, maxQuantile, prCrossPCA)
