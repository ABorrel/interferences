import analyseDB
import predictInterference


prMain = "/home/borrela2/interference/"
prresults = "/home/borrela2/interference/spDataAnalysis/"
prclusters = "/home/borrela2/interference/spDataAnalysis/FinalClustering/"

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

lpSMI = ["/home/borrela2/interference/spDataAnalysis/predictions/test.smi"]

predictor.predictlpSMI(lpSMI)

