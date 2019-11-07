import runExternalSoft
import pathFolder
import toolbox
import pubmed

from os import path, listdir
from shutil import rmtree
from numpy import mean, std
from copy import deepcopy
from re import search
from random import shuffle


class Model:

    def __init__(self, pdesc, pAC50, pAC50All, typeQSAR, corval, maxQuantile, splitRatio, nbCV, ratioAct, nbNA, cell, lchannels, prresult):

        self.corval = corval
        self.maxQauntile = maxQuantile
        self.prresult = prresult
        self.pdesc = pdesc
        self.pAC50 = pAC50
        self.pAC50All = pAC50All
        self.splitRatio = splitRatio
        self.ratioAct = ratioAct
        self.nbCV = nbCV
        self.typeQSAR = typeQSAR
        self.lchannel = lchannels
        self.cell = cell
        self.nbNA = nbNA


    def prepDataColor(self):

        # format by type of AC50
        # change self with one folder by type of AC50

        presult = pathFolder.createFolder(self.prresult + self.cell + "/")
        pClass = presult + "AC50_" + str(self.cell)

        # by pass
        if path.exists(presult + "trainSet.csv") and path.exists(presult + "testSet.csv") and path.exists(pClass):
            dtrain = {}
            dtrain[self.cell] = presult + "trainSet.csv"

            dtest = {}
            dtest[self.cell] = presult + "testSet.csv"

            self.dptrain = dtrain
            self.dptest = dtest

            self.dpAC50 = {}
            self.dpAC50[self.cell] = pClass

            self.dpresult = {}
            self.dpresult[self.cell] = presult
            return 0


        

        color = self.cell + "_n"
        dAC50 = toolbox.loadMatrix(self.pAC50All, sep = "\t")

        fclass = open(pClass, "w")
        fclass.write("CAS\tAff\n")

        lCASID = list(dAC50.keys())[1:]# remove ""
        shuffle(lCASID)

        lact = []
        linact = []
        for CASID in lCASID:
            flagAct = 0
            for channel in list(dAC50[CASID].keys()):
                if search(color, channel):
                    #print dAC50[CASID][channel]
                    if dAC50[CASID][channel] != "NA":
                        lact.append(str(CASID) + "\t1")
                        flagAct = 1
                        break
            if flagAct == 0:
                linact.append(str(CASID) + "\t0")

        nbinact = int(100 * len(lact) / (100 * self.ratioAct)) - len(lact)

        lw = lact + linact[:nbinact]
        shuffle(lw)

        fclass.write("\n".join(lw))
        fclass.close()

        runExternalSoft.prepDataQSAR(self.pdesc, pClass, presult, self.corval, self.maxQauntile, self.splitRatio, self.nbNA, "0")

        dtrain = {}
        dtrain[self.cell] = presult + "trainSet.csv"

        dtest = {}
        dtest[self.cell] = presult + "testSet.csv"

        self.dptrain = dtrain
        self.dptest = dtest

        self.dpAC50 = {}
        self.dpAC50[self.cell] = pClass

        self.dpresult = {}
        self.dpresult[self.cell] = presult




    def prepDataCrossColor(self):

        # format by type of AC50
        # change self with one folder by type of AC50

        presult = pathFolder.createFolder(self.prresult + "crossColor/")
        pClass = presult + "AC50_" + str(self.cell)
        # by pass
        if path.exists(presult + "trainSet.csv") and path.exists(presult + "testSet.csv") and path.exists(pClass):
            dtrain = {}
            dtrain[self.cell] = presult + "trainSet.csv"

            dtest = {}
            dtest[self.cell] = presult + "testSet.csv"

            self.dptrain = dtrain
            self.dptest = dtest

            self.dpAC50 = {}
            self.dpAC50[self.cell] = pClass

            self.dpresult = {}
            self.dpresult[self.cell] = presult

            return 0


        lcolors = ["blue_n", "green_n", "red_n"]
        dAC50 = toolbox.loadMatrix(self.pAC50All, sep="\t")



        fclass = open(pClass, "w")
        fclass.write("CAS\tAff\n")

        lCASID = list(dAC50.keys())[1:]  # remove ""
        shuffle(lCASID)

        lact = []
        linact = []
        for CASID in lCASID:
            flag = 0
            for color in lcolors:
                #print flag
                if flag == 4:
                    break
                else:
                    flag = 0
                for channel in list(dAC50[CASID].keys()):
                    if search(color, channel):
                        #print color, channel
                        if dAC50[CASID][channel] != "NA":
                            flag = flag + 1

            if flag == 4:
                lact.append(str(CASID) + "\t1")
            else:
                linact.append(str(CASID) + "\t0")

        nbinact = int(100 * len(lact) / (100 * self.ratioAct)) - len(lact)

        lw = lact + linact[:nbinact]
        shuffle(lw)

        fclass.write("\n".join(lw))
        fclass.close()

        runExternalSoft.prepDataQSAR(self.pdesc, pClass, presult, self.corval, self.maxQauntile, self.splitRatio, self.nbNA)

        dtrain = {}
        dtrain[self.cell] = presult + "trainSet.csv"

        dtest = {}
        dtest[self.cell] = presult + "testSet.csv"

        self.dptrain = dtrain
        self.dptest = dtest

        self.dpAC50 = {}
        self.dpAC50[self.cell] = pClass

        self.dpresult = {}
        self.dpresult[self.cell] = presult



    def prepData(self, typeData):
        # format by type of AC50
        # change self with one folder by type of AC50

        dAC50 = toolbox.loadMatrix(self.pAC50, sep = "\t")

        dfileAC50 = {}
        dprresult = {}


        imax = len(self.lchannel)
        i = 0
        while i < imax:
            AC50type = self.lchannel[i]
            presult = pathFolder.createFolder(self.prresult + AC50type + "/")
            dprresult[AC50type] = presult

            dfileAC50[AC50type] = open(presult + "AC50_" + str(AC50type), "w")
            dfileAC50[AC50type].write("CAS\tAff\n")

            i += 1

        for CAS in list(dAC50.keys()):
            for channel in self.lchannel:
                dfileAC50[channel].write(str(CAS) + "\t" + str(dAC50[CAS][channel]) + "\n")


        for typeAC50 in self.lchannel:
            dfileAC50[typeAC50].close()
            dfileAC50[typeAC50] = dfileAC50[typeAC50].name

        self.dpAC50 = dfileAC50
        self.dpresult = dprresult

        dtrain = {}
        dtest = {}
        for typeAC50 in list(self.dpAC50.keys()):
            if self.typeQSAR == "Reg":
                runExternalSoft.prepDataQSAR(self.pdesc, self.dpAC50[typeAC50], self.dpresult[typeAC50], self.corval, self.maxQauntile, self.splitRatio, self.nbNA)

            else:
                if typeData == "all":
                    self.writeClass()
                elif typeData == "active":
                    self.writeClassActive()

                ptrain = self.dpresult[typeAC50] + "trainSet.csv"
                ptest = self.dpresult[typeAC50] + "testSet.csv"

                print(ptrain)
                print(ptest)

                if not path.exists(ptrain) and not path.exists(ptest):
                    runExternalSoft.prepDataQSAR(self.pdesc, self.dpAC50[typeAC50], self.dpresult[typeAC50], self.corval, self.maxQauntile,
                                             self.splitRatio, self.nbNA)

            dtrain[typeAC50] = ptrain
            dtest[typeAC50] = ptest


        self.dptrain = dtrain
        self.dptest = dtest


    def buildQSARReg(self):

        for typeAC50 in self.dpAC50:
            runExternalSoft.QSARReg(self.dptrain[typeAC50], self.dptest[typeAC50], "0", self.dpresult[typeAC50], self.nbCV)



    def buildQSARClass(self):

        for typeAC50 in self.dpAC50:
            #if not path.exists(self.dpresult[typeAC50] + "perf.txt") and not path.exists(self.dpresult[typeAC50] + "perfCV.txt"):
            runExternalSoft.QSARClass(self.dptrain[typeAC50], self.dptest[typeAC50], self.dpresult[typeAC50], self.nbCV)




    def writeClassActive(self):

        from random import shuffle

        print(self.dpresult)
        print(self.pAC50All)
        dAC50All = toolbox.loadMatrix(self.pAC50All, sep="\t")

        for typeAC50 in self.dpresult:
            pclass = self.dpresult[typeAC50] + "actClass.txt"

            if path.exists(pclass) and path.getsize(pclass) > 10:
                self.dpAC50[typeAC50] = pclass

            else:
                filin = open(self.dpAC50[typeAC50], "r")
                llines = filin.readlines()
                filin.close()

                filout = open(pclass, "w")
                filout.write(llines[0])

                # shuffle lines
                llines = llines[1:]
                shuffle(llines)

                nbact = 0
                for lineChem in llines:
                    AC50 = lineChem.strip().split("\t")[-1]
                    if AC50 != "NA":
                        nbact = nbact + 1

                nbinact = int(100 * nbact / (100 * self.ratioAct)) - nbact

                # select active chemical
                llineAct = []
                for lineChem in llines[1:]:
                   lAC50 = lineChem.strip().split("\t")
                   lnew = [lAC50[0]]
                   for AC50 in lAC50[1:]:
                       if AC50 != "NA":
                           lnew.append("1")
                           llineAct.append("\t".join(lnew))

                # select inactive but select active for other channel
                if typeAC50 != "Luc_IC50":
                # add channel active in the set
                   llineInact = []
                   for CASID in list(dAC50All.keys()):
                       if dAC50All[CASID][self.cell + "_" + typeAC50] != "NA":
                           continue
                       else:
                           for channel in list(dAC50All[CASID].keys()):
                               if not search("Luc_IC50", channel):
                                   if channel != "CASID":
                                       if dAC50All[CASID][channel] != "NA":
                                           lnew = [CASID, "0"]
                                           llineInact.append("\t".join(lnew))
                                           break

                # random active
                nbinactselected = len(llineInact)
                #print nbinact, nbinactselected

                if nbinactselected >= nbinact:
                    shuffle(llineInact)
                    llineInact = llineInact[:nbinact]
                    lw = llineAct + llineInact
                    shuffle(lw)
                else:
                    nwinact = nbinactselected

                    # first loop to take inactive
                    for lineChem in llines[1:]:
                        lAC50 = lineChem.strip().split("\t")
                        lnew = [lAC50[0]]
                        for AC50 in lAC50[1:]:
                            if AC50 == "NA":
                                lnew.append("0")
                                lneww = "\t".join(lnew)
                                if not lneww in llineInact:
                                    llineInact.append(lneww)
                                    nwinact += 1
                                    break

                        if nwinact >= nbinact:
                            break
                lw = llineAct + llineInact
                shuffle(lw)

                filout.write("\n".join(lw))
                filout.close()
                self.dpAC50[typeAC50] = pclass


    def writeClass(self):
        from random import shuffle

        for typeAC50 in self.dpresult:
            pclass = self.dpresult[typeAC50] + "actClass.txt"
            #print pclass

            if path.exists(pclass) and path.getsize(pclass) > 10:
                self.dpAC50[typeAC50] = pclass

            else:
                filin = open(self.dpAC50[typeAC50], "r")
                llines = filin.readlines()
                filin.close()

                filout = open(pclass, "w")
                filout.write(llines[0])

                # shuffle lines
                llines = llines[1:]
                shuffle(llines)

                nbact = 0
                for lineChem in llines:
                    AC50 = lineChem.strip().split("\t")[-1]
                    if AC50 != "NA" and AC50 != "0":
                        nbact = nbact + 1

                if self.ratioAct == 1:
                    nbinact = len(llines) - nbact
                else:
                    nbinact = int(100*nbact/(100*self.ratioAct)) - nbact

                nwinact = 0
                flaginact = 0
                # first loop to take active
                for lineChem in llines[1:]:
                    lAC50 = lineChem.strip().split("\t")
                    lnew = [lAC50[0]]
                    for AC50 in lAC50[1:]:
                        if AC50 == "NA" or AC50 == "0":
                            lnew.append("0")
                            nwinact += 1
                            flaginact = 1
                        else:
                            lnew.append("1")
                            flaginact = 0

                    if flaginact == 1 and nwinact <= nbinact:
                        filout.write("\t".join(lnew) + "\n")
                    elif flaginact == 0:
                        filout.write("\t".join(lnew) + "\n")
                filout.close()
                self.dpAC50[typeAC50] = pclass


def runQSARClassForPubChem(passay, PRPUBCHEM, PRDESC, PRSMI, corval, maxQuantile, splitratio, nbCV, ratioAct, nbNA, nameCell, lchannels, typeData,  prout):


    passay = pubmed.formatPubChemTable(passay, PRPUBCHEM, prout)
    print(passay)

    lpdesc = pubmed.computeDesc(passay, PRDESC, PRSMI, prout, nbfile=2, update=0)

    prQSAR = pathFolder.createFolder(prout + "QSAR/")
    cModel = Model(lpdesc[0], lpdesc[1], "", "class", corval, maxQuantile, splitratio,
                   nbCV, ratioAct, nbNA, nameCell, lchannels, prQSAR)
    cModel.prepData(typeData = typeData)
    cModel.buildQSARClass()





def runQSARClass(cDesc, cAssay, pAC50All, corval, maxQuantile, splitratio, nbCV, ratioAct, nbrepeat, nbNA, nameCell, lchannels, typeData,  prout):

    for i in range(1, nbrepeat + 1):
        prQSAR = prout + str(i) + "/"
        #rmtree(prQSAR)############################################################################### to remove
        pathFolder.createFolder(prQSAR)

        cModel = Model(cDesc.pdesc1D2D, cAssay.pAC50, pAC50All, "class", corval, maxQuantile, splitratio,
                                        nbCV, ratioAct, nbNA, nameCell, lchannels, prQSAR)
        if typeData == "color":
            cModel.prepDataColor()
        elif typeData == "crosscolor":
            cModel.prepDataCrossColor()
        else:
            cModel.prepData(typeData)
        cModel.buildQSARClass()


    prQSARAV = pathFolder.createFolder(prout + "Average/")
    prQSARProb = pathFolder.createFolder(prout + "Prob/")
    mergeProba(prout, "RF", prQSARProb)
    mergeResults(prout, prQSARAV)
    prDescAV = pathFolder.createFolder(prout + "descImportance/")
    mergeDescInvolve(prout, "LDA", 10, prDescAV)
    mergeDescInvolve(prout, "RF", 10,  prDescAV)

def mergeDescInvolve(prin, ML, nbdesc, prout):

    dimportance = {}

    lprrun = listdir(prin)
    for prrun in lprrun:
        if prrun == "Average" or prrun == "descImportance" or prrun == "Prob":
            continue

        lprcell = listdir(prin + "/" + prrun + "/")
        for prcell in lprcell:
            if not prcell in list(dimportance.keys()):
                dimportance[prcell] = {}
            dimportance[prcell][prrun] = {}
            pimportance = prin + prrun + "/" + prcell + "/" + str(ML) + "class/ImportanceDesc"

            if path.exists(pimportance):
                ddescimportance = toolbox.loadMatrix(pimportance , sep="\t")
                dimportance[prcell][prrun] = ddescimportance
            else:
                pmodel = prin + prrun + "/" + prcell + "/" + str(ML) + "class/model.RData"
                ptrain = prin + prrun + "/" + prcell + "/trainSet.csv"
                if path.exists(pmodel):
                    runExternalSoft.createImportanceTable(pmodel, ML, ptrain, prin + prrun + "/" + prcell + "/" + str(ML) + "class/")
                    ddescimportance = toolbox.loadMatrix(pimportance, sep="\t")
                    dimportance[prcell][prrun] = ddescimportance


    # write global table
    for typeAssay in list(dimportance.keys()):
        pdesc = prout + "Importance" + str(ML) + "_" + typeAssay
        fdesc = open(pdesc, "w")
        lrun = list(dimportance[typeAssay].keys())

        ldesc = list(dimportance[typeAssay][lrun[0]].keys())
        fdesc.write("Desc\tRun\tval\n")
        for desc in ldesc:
            for run in lrun:
                try: fdesc.write(desc + "\t" + str(run) + "\t" + str(dimportance[typeAssay][run][desc]["x"]) + "\n")
                except: fdesc.write(desc + "\t" + str(run) + "\t0.0\n")
        fdesc.close()

        runExternalSoft.runImportanceDesc(pdesc, nbdesc)

    return 0



def mergeProba(prin, ML, prout):


    dprob = {}
    dreal = {}

    lprrun = listdir(prin)
    for prrun in lprrun:
        if prrun == "Average" or prrun == "descImportance" or prrun == "Prob":
            continue

        lprcell = listdir(prin + "/" + prrun + "/")
        for prcell in lprcell:
            if not prcell in list(dreal.keys()):
                dreal[prcell] = {}
            flag = 0
            for filin in listdir(prin + prrun + "/" + prcell + "/"):
                if search("AC50_", filin):
                    paff = prin + prrun + "/" + prcell + "/" + filin
                    flag = 1
                    break
            daff = toolbox.loadMatrix(paff, sep ="\t")
            dreal[prcell].update(deepcopy(daff))

            if not prcell in list(dprob.keys()):
                dprob[prcell] = {}
            dprob[prcell][prrun] = {}
            pCV = prin + prrun + "/" + prcell + "/" + str(ML) + "class/PerfRFClassCV10.txt"
            dCV = toolbox.loadMatrix(pCV)
            dprob[prcell][prrun]["CV"] = dCV

            ptrain = prin + prrun + "/" + prcell + "/" + str(ML) + "class/classTrain.csv"
            dtrain = toolbox.loadMatrix(ptrain, sep = ",")
            dprob[prcell][prrun]["train"] = dtrain

            ptest = prin + prrun + "/" + prcell + "/" + str(ML) + "class/classTest.csv"
            dtest = toolbox.loadMatrix(ptest, sep = ",")
            dprob[prcell][prrun]["test"] = dtest


    # write table for probability
    dw = {}
    for prcell in list(dprob.keys()):
        dw[prcell] = {}
        dw[prcell] = {}
        dw[prcell]["train"] = {}
        dw[prcell]["test"] = {}
        dw[prcell]["CV"] = {}

        for run in list(dprob[prcell].keys()):
            for IDtrain in list(dprob[prcell][run]["train"].keys()):
                if not IDtrain in list(dw[prcell]["train"].keys()):
                    dw[prcell]["train"][IDtrain] = []
                dw[prcell]["train"][IDtrain].append(float(dprob[prcell][run]["train"][IDtrain]["x"]))


            for IDtest in dprob[prcell][run]["test"]:
                if not IDtest in list(dw[prcell]["test"].keys()):
                    dw[prcell]["test"][IDtest] = []
                dw[prcell]["test"][IDtest].append(float(dprob[prcell][run]["test"][IDtest]["x"]))

            for IDCV in dprob[prcell][run]["CV"]:
                if not IDCV in list(dw[prcell]["CV"].keys()):
                    dw[prcell]["CV"][IDCV] = []
                dw[prcell]["CV"][IDCV].append(float(dprob[prcell][run]["CV"][IDCV]["Predict"]))


    for prcell in list(dw.keys()):

        # train
        pfiloutTrain = prout + prcell + "_train"
        filoutTrain = open(pfiloutTrain, "w")
        filoutTrain.write("ID\tMpred\tSDpred\tReal\n")
        for IDtrain in list(dw[prcell]["train"].keys()):
            try:filoutTrain.write("%s\t%.3f\t%.3f\t%s\n"%(IDtrain, mean(dw[prcell]["train"][IDtrain]), std(dw[prcell]["train"][IDtrain]), dreal[prcell][IDtrain]["Aff"]))
            except:
                print(dw[prcell]["train"][IDtrain])
                print(dreal[prcell][IDtrain]["Aff"])
                ddd
        filoutTrain.close()

        runExternalSoft.plotAC50VSProb(pfiloutTrain, prout)

        #test
        pfiloutTest = prout + prcell + "_test"
        filoutTest = open(pfiloutTest, "w")
        filoutTest.write("ID\tMpred\tSDpred\tReal\n")
        for IDtest in list(dw[prcell]["test"].keys()):
            filoutTest.write("%s\t%.3f\t%.3f\t%s\n"%(IDtest, mean(dw[prcell]["test"][IDtest]), std(dw[prcell]["test"][IDtest]), dreal[prcell][IDtest]["Aff"]))
        filoutTest.close()

        runExternalSoft.plotAC50VSProb(pfiloutTest, prout)

        #CV
        pfiloutCV = prout + prcell + "_CV"
        filoutCV = open(pfiloutCV, "w")
        filoutCV.write("ID\tMpred\tSDpred\tReal\n")
        for IDCV in list(dw[prcell]["CV"].keys()):
            filoutCV.write("%s\t%.3f\t%.3f\t%s\n" % (
            IDCV, mean(dw[prcell]["CV"][IDCV]), std(dw[prcell]["CV"][IDCV]), dreal[prcell][IDCV]["Aff"]))
        filoutCV.close()

        runExternalSoft.plotAC50VSProb(pfiloutCV, prout)

    return 0




def mergeResults(prin, prout):

    dresult = {}
    dperf = {}
    dperf["Acc"] = []
    dperf["Sp"] = []
    dperf["Se"] = []
    dperf["MCC"] = []

    lprrun = listdir(prin)
    for prrun in lprrun:
        if prrun == "Average" or prrun == "descImportance":
            continue
        lprcell = listdir(prin + "/" + prrun + "/")
        for prcell in lprcell:
            pperfCV = prin + "/" + prrun + "/" + prcell + "/perfCV.csv"
            pperftrain = prin + "/" + prrun + "/" + prcell + "/perfTrain.csv"
            pperftest = prin + "/" + prrun + "/" + prcell + "/perfTest.csv"

            try:
                MCV = toolbox.loadMatrix(pperfCV, sep=",")
                Mtrain = toolbox.loadMatrix(pperftrain, sep=",")
                Mtest = toolbox.loadMatrix(pperftest, sep=",")
            except:
                continue


            lML = list(MCV.keys())
            lcriteria = list(dperf.keys())
            lset = ["CV", "train", "test"]
            # create the structures
            if not prcell in list(dresult.keys()):
                dresult[prcell] = {}
                dresult[prcell]["CV"] = {}
                dresult[prcell]["train"] = {}
                dresult[prcell]["test"] = {}
                for ML in lML:
                    dresult[prcell]["CV"][ML] = deepcopy(dperf)
                    dresult[prcell]["train"][ML] = deepcopy(dperf)
                    dresult[prcell]["test"][ML] = deepcopy(dperf)

            for ML in lML:
                for criteria in lcriteria:
                    dresult[prcell]["CV"][ML][criteria].append(float(MCV[ML][criteria]))
                    dresult[prcell]["train"][ML][criteria].append(float(Mtrain[ML][criteria]))
                    dresult[prcell]["test"][ML][criteria].append(float(Mtest[ML][criteria]))


    dout = deepcopy(dresult)
    for celltype in list(dresult.keys()):
        for set in list(dresult[celltype].keys()):
            for ML in list(dresult[celltype][set].keys()):
                for criteria in  list(dresult[celltype][set][ML].keys()):

                    AV = round(mean(dresult[celltype][set][ML][criteria]),3)
                    SD = round(std(dresult[celltype][set][ML][criteria]),3)

                    dout[celltype][set][ML][criteria] = [AV, SD]


    # write result
    lperfcriteria = ["Acc", "Sp", "Se", "MCC"]
    for celltype in list(dout.keys()):
        pfilout = prout + celltype + ".csv"
        filout = open(pfilout, "w")
        for set in list(dout[celltype].keys()):
            filout.write(str(set) + "\n")
            filout.write("\t" + "\t".join(["M-" + str(c) + "\t" + "SD-" + str(c) for c in lperfcriteria]) + "\n")
            for ML in list(dout[celltype][set].keys()):
                filout.write(ML)
                for criteria in lperfcriteria:
                    filout.write("\t" + str(dout[celltype][set][ML][criteria][0]) + "\t" + str(dout[celltype][set][ML][criteria][1]))
                filout.write("\n")
        filout.close()


    return 0
