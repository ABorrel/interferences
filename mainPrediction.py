import analyseDB



prMain = "/home/borrela2/interference/"
prresults = "/home/borrela2/interference/spDataAnalysis/"



prSMI = prMain + "SMI/"
prDesc = prMain + "Desc/"
prlogDesc = prMain + "log/"
prPNG = prMain + "PNG/"

corval = 0.90
maxquantile = 90

cDesc = analyseDB.Descriptors(prSMI, prDesc, prPNG, prresults, prlogDesc)
cDesc.loadClassForPred(corval, maxquantile)