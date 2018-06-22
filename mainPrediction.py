from os import listdir
from re import search

import analyseDB
import predictInterference


prMain = "/home/borrela2/interference/"
prresults = "/home/borrela2/interference/spDataAnalysis/"
prclusters = "/home/borrela2/interference/spDataAnalysis/FinalClustering/"
pAC50all = "/home/borrela2/interference/spDataAnalysis/AC50_all"

prpredict = "/home/borrela2/interference/spDataAnalysis/predictions/"

prSMI = prMain + "SMI/"
prDesc = prMain + "Desc/"
prlogDesc = prMain + "log/"
prPNG = prMain + "PNG/"

corval = 0.90
maxquantile = 90
distMeth = "euclidian"
distAgg = "ward.D2"

cDesc = analyseDB.Descriptors(prSMI, prDesc, prPNG, prresults, prlogDesc)
cDesc.loadClassForPred(corval, maxquantile)

lcasID = []
predictor = predictInterference.predictor(cDesc, prclusters, lcasID, prpredict, distMeth, distAgg)
#predictor.summarizePredictor()

predictor.validationPredictor(pAC50all)

#lpSMI = [prpredict + "smi/" + i for i in listdir(prpredict + "smi/")]
lpSMI = [prpredict + "luciferin.smi"]

#predictor.predictlpSMI(lpSMI)

