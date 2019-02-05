import QSARModel
import pathFolder





PRMAIN = "/home/borrela2/interference/"
PRTESTING = "/home/borrela2/interference/testing/"
PRMODELS = PRTESTING + "QSARmodel/"
PRPUBCHEM = PRMAIN + "PUBCHEM/"
PRPUBCHEMSDF = PRMAIN + "PUBCHEMSDF/"

PRDESC = PRTESTING + "DESC/"
PRSMI = PRTESTING + "SMI/"



passays = "/home/borrela2/interference/testing/data/AID_411_luc.csv"
corval = 0.85
maxQuantile = 90
splitratio = 0.15
nbCV = 10
nbNA = 10000
ratioAct = 1
nbrepeat = 1
nameCell = "IC50"
lchannels = ["Aff"]
typeData = "all"
prout = pathFolder.createFolder(PRPUBCHEM + "411/")


QSARModel.runQSARClassForPubChem(passays, PRPUBCHEMSDF, PRDESC, PRSMI, corval, maxQuantile, splitratio, nbCV, ratioAct, nbNA, nameCell, lchannels, typeData,  prout)

