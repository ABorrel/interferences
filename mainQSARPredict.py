import toolbox
import QSARpredictor
import pubmed
import pathFolder
import chemical
import runExternalSoft

from os import path, listdir, remove
from re import search
from shutil import move
from random import shuffle
from numpy import mean, std


PRMAIN = "/home/borrela2/interference/"
PRTESTING = "/home/borrela2/interference/testing/"
PRMODELS = PRTESTING + "QSARmodel/"
PRPUBCHEMSDF = PRMAIN + "PUBCHEMSDF/"

PRDESC = PRTESTING + "DESC/"
PRSMI = PRTESTING + "SMI/"



def formatPRModels(PRMODELS):

    lf1 = listdir(PRMODELS)
    lfile = []
    for f1 in lf1:
        for f2 in listdir(PRMODELS + f1):
            print f2, "f2"
            lfile.append(PRMODELS + f1 + "/" + f2)
            for f3 in listdir(PRMODELS + f1 + "/" + f2):
                lfile.append(PRMODELS + f1 + "/" + f2 + "/" + f3)
                print f3, "f3"
                if path.exists(PRMODELS + f1 + "/" + f2 + "/" + f3 + "/"):
                    for f4 in listdir(PRMODELS + f1 + "/" + f2 + "/" + f3):
                        print f4, "f4"
                        lfile.append(PRMODELS + f1 + "/" + f2 + "/" + f3 + "/" + f4)
                        if path.exists(PRMODELS + f1 + "/" + f2 + "/" + f3 + "/" + f4 + "/"):
                            for f5 in listdir(PRMODELS + f1 + "/" + f2 + "/" + f3 + "/" + f4):
                                print f5, "f5"
                                #lfile.append(PRMODELS + f1 + "/" + f2 + "/" + f3 + "/" + f4 + "/" + f5)
                                if not f5 == "model.RData":
                                    print "del", PRMODELS + f1 + "/" + f2 + "/" + f3 + "/" + f4 + "/" + f5
                                    remove(PRMODELS + f1 + "/" + f2 + "/" + f3 + "/" + f4 + "/" + f5)
                                else:
                                    print "rename", PRMODELS + f1 + "/" + f2 + "/" + f3 + "/" + f4 + "/model" + f2 + ".Rdata"
                                    move(PRMODELS + f1 + "/" + f2 + "/" + f3 + "/" + f4 + "/model.RData", PRMODELS + f1 + "/" + f2 + "/" + f3 + "/" + f4 + "/model" + f2 + ".Rdata")

    for f in lfile:
        if search("png", f) or search("csv", f):
            remove(f)





def applyModel(pdesc, prmodels, prout):

    lfolderModel = listdir(prmodels)
    dmodel = {}
    for folderModel in lfolderModel:
        dmodel[folderModel] = {}
        lmodel = listdir(prmodels + folderModel)
        for model in lmodel:
            dmodel[folderModel][model] = {}
            prmodel = prmodels + folderModel + "/" + model
            for ML in listdir(prmodel):
                prmodelML = prmodels + folderModel + "/" + model + "/" + ML + "/"
                lmodelR = listdir(prmodelML)
                dmodel[folderModel][model][ML] = [prmodels + folderModel + "/" + model + "/" + ML + "/" + modelR for modelR in lmodelR]



    for typeModel in dmodel.keys():
        pathFolder.createFolder(prout + typeModel + "/")
        for color in dmodel[typeModel].keys():
            pathFolder.createFolder(prout + typeModel + "/" + color + "/")
            dperf = {}
            for ML in dmodel[typeModel][color].keys():
                proutbyRmodel = pathFolder.createFolder(prout + typeModel + "/" + color + "/" + ML + "/")
                lpredict = []
                for modelR in dmodel[typeModel][color][ML]:
                    ppredict = runExternalSoft.predictModel(pdesc, modelR, ML, proutbyRmodel)
                    lpredict.append(ppredict)


                dprob = {}
                for ppredict in lpredict:
                    try: dpredict = toolbox.loadMatrix(ppredict, sep=",")
                    except: continue
                    for chem in dpredict.keys():
                        if not chem in dprob:
                            dprob[chem] = {}
                            dprob[chem]["Pred"] = []
                            dprob[chem]["Aff"] = dpredict[chem]["Aff"]
                        if dpredict[chem]["Pred"] != "NA":
                            dprob[chem]["Pred"].append(float(dpredict[chem]["Pred"]))

                pfsumML = proutbyRmodel + "sumProb"
                fsumML = open(pfsumML, "w")
                fsumML.write("ID,Mpred,SDpred,Real\n")
                for chem in dprob.keys():
                    if len(dprob[chem]["Pred"]) == 0:
                        fsumML.write("%s,NA,NA,%s\n"%(chem, dprob[chem]["Aff"]))
                    else:
                        fsumML.write("%s,%f,%f,%s\n"%(chem, mean(dprob[chem]["Pred"]), std(dprob[chem]["Pred"]), dprob[chem]["Aff"]))
                fsumML.close()
                pquality = runExternalSoft.qualityPred(pfsumML)
                runExternalSoft.plotAC50VSProb(pfsumML, proutbyRmodel)
                dperf[ML] = toolbox.loadMatrix(pquality, sep = ",")

            pfsum = prout + typeModel + "/" + color + "/sumPerf.csv"
            fsum = open(pfsum, "w")
            lh = ["TP", "TN", "FP", "FN", "acc", "se", "sp", "mcc", "MpbTP", "SDpbTP", "MpbTN", "SDpbTN", "MpbFP", "SDpbFP", "MpbFN", "SDpbFN"]
            fsum.write("Model," + ",".join(lh) + "\n")
            for ML in dperf.keys():
                lw = [dperf[ML][h]["x"]for h in lh]
                i = 0
                imax = len(lw)
                while i<imax:
                    try:
                        lw[i] = str(round(float(lw[i]),2))
                    except:
                        lw[i] = str(lw[i])
                    i += 1
                fsum.write("%s,%s\n"%(ML, ",".join(lw)))
            fsum.close()


def computePCACross(pdesc, pdescModel, paff, prout):

    runExternalSoft.drawPCACross2Set(pdescModel, paff, pdesc, prout)





for assay in listdir("/home/borrela2/interference/testing/data/"):
    #print assay
    #continue
    #assay = "AID_589_4MU.csv"

    dPubChemBioassays = "/home/borrela2/interference/testing/data/" + assay
    prout = pathFolder.createFolder(PRTESTING + assay[4:-4] + "/")
    passay = pubmed.formatPubChemTable(dPubChemBioassays, PRPUBCHEMSDF, prout, update=0)
    pdesc = pubmed.computeDesc(passay, PRDESC, PRSMI, prout, update=0)

    pdescModel = "/home/borrela2/interference/spDataAnalysis/CrossPCA/descClean.csv"
    paff = "/home/borrela2/interference/spDataAnalysis/AC50_all"
    prPCA = pathFolder.createFolder(prout + "PCA/")
    computePCACross(pdesc, pdescModel, paff, prPCA)

    applyModel(pdesc, PRMODELS, prout)






