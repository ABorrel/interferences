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
                                    pathFolder.createFolder(PRMODELS + f1 + "/" + f3 + "/" + f4 + "/")
                                    move(PRMODELS + f1 + "/" + f2 + "/" + f3 + "/" + f4 + "/model.RData", PRMODELS + f1 + "/" + f3 + "/" + f4 + "/model" + f2 + ".Rdata")

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



def combineResult(cIn, prout):

    lfolders = listdir(prout)
    lmodel = ["HEK293", "HepG2", "QSARclassColor", "QSARclassCrossColor"]
    dpred = {}
    for folder in lfolders:
        if folder in lmodel:
            lsubfolder = listdir(prout + folder + "/")
            for subfolder in lsubfolder:
                lsubsubfolder = listdir(prout + folder + "/" + subfolder + "/")
                for subsubfolder in lsubsubfolder:
                    if subsubfolder == "RFclass":
                        if not folder in dpred.keys():
                            dpred[folder] = {}
                        if not subfolder in dpred[folder].keys():
                            dpred[folder][subfolder] = {}
                        dpred[folder][subfolder][subsubfolder] = toolbox.loadMatrixToDict(prout + folder + "/" + subfolder + "/" + subsubfolder + "/sumProb", sep = ",")


    dw = {}
    for k1 in dpred.keys():
        dw[k1] = {}
        for k2 in dpred[k1].keys():
            for k3 in dpred[k1][k2].keys():
                if not k3 in dw[k1].keys():
                    dw[k1][k3] = {}
                if not k2 in dw[k1][k3].keys():
                    dw[k1][k3][k2] = {}
                    dw[k1][k3][k2]["Prob"] = []
                    dw[k1][k3][k2]["ID"] = []
                    dw[k1][k3][k2]["Abs"] = []


                labs = cIn.values()
                labs.sort()

                ltemp = []
                for Abs in labs:
                    for ID in cIn.keys():
                        if Abs == cIn[ID] and ID not in ltemp:
                            dw[k1][k3][k2]["Prob"].append(str(dpred[k1][k2][k3][ID]["Mpred"]))
                            dw[k1][k3][k2]["ID"].append(ID)
                            dw[k1][k3][k2]["Abs"].append(str(cIn[ID]))
                            ltemp.append(ID)


    for k1 in dw.keys():
        for k2 in dw[k1].keys():
            pfiloutProb = prout + k1 + "_" + k2 + "predict_prob.csv"
            pfiloutPred = prout + k1 + "_" + k2 + "predict.csv"

            filoutProb = open(pfiloutProb, "w")
            filoutPred = open(pfiloutPred, "w")

            lk3 = dw[k1][k2].keys()
            filoutPred.write("ID\tAbs\t" + "\t".join(lk3) + "\n")
            filoutProb.write("ID\tAbs\t" + "\t".join(lk3) + "\n")
            i = 0
            imax = len(dw[k1][k2][dw[k1][k2].keys()[0]]["ID"])
            while i < imax:
                filoutPred.write(dw[k1][k2][lk3[0]]["ID"][i] + "\t" + dw[k1][k2][lk3[0]]["Abs"][i])
                filoutProb.write(dw[k1][k2][lk3[0]]["ID"][i] + "\t" + dw[k1][k2][lk3[0]]["Abs"][i])
                for k3 in lk3:
                    filoutProb.write("\t" + dw[k1][k2][k3]["Prob"][i])

                    if dw[k1][k2][k3]["Prob"][i] == "NA":
                        filoutPred.write("\tNA")
                    elif float(dw[k1][k2][k3]["Prob"][i]) < 0.5:
                        filoutPred.write("\t0")
                    else:
                        filoutPred.write("\t1")

                i += 1

                filoutProb.write("\n")
                filoutPred.write("\n")

            filoutPred.close()
            filoutProb.close()






# prepare model R
#formatPRModels(PRMODELS)


for assay in listdir("/home/borrela2/interference/testing/data/"):
    #print assay
    #continue
    #assay = "AID_589_4MU.csv"
    # for PCA
    pdescModel = "/home/borrela2/interference/spDataAnalysis/CrossPCA/descClean.csv"
    paff = "/home/borrela2/interference/spDataAnalysis/AC50_all"

    # case of PUBMED
    if search("AID_", assay):
        continue
        dPubChemBioassays = "/home/borrela2/interference/testing/data/" + assay
        prout = pathFolder.createFolder(PRTESTING + assay[4:-4] + "/")
        passay = pubmed.formatPubChemTable(dPubChemBioassays, PRPUBCHEMSDF, prout, update=0)
        pdesc = pubmed.computeDesc(passay, PRDESC, PRSMI, prout, update=0)


        prPCA = pathFolder.createFolder(prout + "PCA/")
        computePCACross(pdesc, pdescModel, paff, prPCA)

        applyModel(pdesc, PRMODELS, prout)
    elif search("Spectre", assay):
        pBioassays = "/home/borrela2/interference/testing/data/" + assay
        prout = pathFolder.createFolder(PRTESTING + assay.split(".") [0] + "/")

        import photoChemDB
        cSpectrum = photoChemDB.photoChem(pBioassays, [], prout)
        pdesc = cSpectrum.importDescriptors()

        prPCA = pathFolder.createFolder(prout + "PCA/")
        #computePCACross(pdesc, pdescModel, paff, prPCA)
        #applyModel(pdesc, PRMODELS, prout)

        combineSpectreResult(cSpectrum.dwave, prout)

    elif search("Dye", assay):
        pBioassays = "/home/borrela2/interference/testing/data/" + assay
        prout = pathFolder.createFolder(PRTESTING + assay.split(".")[0] + "/")

        import dyeAnalysis
        cDYE = dyeAnalysis.DYE(pBioassays, [], prout)
        pdesc = cDYE.importDescriptors()

        #prPCA = pathFolder.createFolder(prout + "PCA/")
        #computePCACross(pdesc, pdescModel, paff, prPCA)
        #applyModel(pdesc, PRMODELS, prout)

        combineResult(cDYE.dcolor, prout)



