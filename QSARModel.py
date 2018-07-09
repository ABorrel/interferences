import runExternalSoft
import pathFolder
import toolbox

from os import path, listdir
from numpy import mean, std
from copy import deepcopy
from re import search


class Model:

    def __init__(self, pdesc, pAC50, pAC50All, typeQSAR, corval, maxQuantile, splitRatio, nbCV, ratioAct, cell, lchannels, prresult):

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


    def prepData(self):

        # format by type of AC50
        # change self with one folder by type of AC50

        fAC50 = open(self.pAC50, "r")
        llinesAC50 = fAC50.readlines()
        fAC50.close()

        dfileAC50 = {}
        dprresult = {}
        if self.lchannel == []:
            lAC50type = llinesAC50[0].strip().split("\t")[1:]
        else:
            lAC50type = self.lchannel
        imax = len(lAC50type)
        i = 0
        while i < imax:
            AC50type = lAC50type[i]
            presult = pathFolder.createFolder(self.prresult + AC50type + "/")
            dprresult[AC50type] = presult

            dfileAC50[AC50type] = open(presult + "AC50_" + str(AC50type), "w")
            dfileAC50[AC50type].write("CAS\tAff\n")

            i += 1

        for lineAC50 in llinesAC50[1:]:
            lem = lineAC50.strip().split("\t")
            CAS = lem[0]
            print CAS

            i = 0
            while i < len(lAC50type):
                dfileAC50[lAC50type[i]].write(str(CAS) + "\t" + str(lem[i+1]) + "\n")
                i += 1

        for typeAC50 in self.lchannel:
            dfileAC50[typeAC50].close()
            dfileAC50[typeAC50] = dfileAC50[typeAC50].name

        self.dpAC50 = dfileAC50
        self.dpresult = dprresult

        dtrain = {}
        dtest = {}
        for typeAC50 in self.dpAC50.keys():
            if self.typeQSAR == "Reg":
                runExternalSoft.prepDataQSAR(self.pdesc, self.dpAC50[typeAC50], self.dpresult[typeAC50], self.corval, self.maxQauntile, self.splitRatio)

            else:
                self.writeClassActive()
                runExternalSoft.prepDataQSAR(self.pdesc, self.dpAC50[typeAC50], self.dpresult[typeAC50], self.corval, self.maxQauntile,
                                             self.splitRatio, "0")

            dtrain[typeAC50] = self.dpresult[typeAC50] + "trainSet.csv"
            dtest[typeAC50] = self.dpresult[typeAC50] + "testSet.csv"


        self.dptrain = dtrain
        self.dptest = dtest


    def buildQSARReg(self):

        for typeAC50 in self.dpAC50:
            runExternalSoft.QSARReg(self.dptrain[typeAC50], self.dptest[typeAC50], "0", self.dpresult[typeAC50], self.nbCV)



    def buildQSARClass(self):

        for typeAC50 in self.dpAC50:
            runExternalSoft.QSARClass(self.dptrain[typeAC50], self.dptest[typeAC50], self.dpresult[typeAC50], self.nbCV)




    def writeClassActive(self):

        from random import shuffle

        print self.dpresult
        print self.pAC50All
        dAC50All = toolbox.loadMatrix(self.pAC50All, sep="\t")

        for typeAC50 in self.dpresult:
            pclass = self.dpresult[typeAC50] + "actClass.txt"
            print pclass

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
                   for CASID in dAC50All.keys():
                       if dAC50All[CASID][self.cell + "_" + typeAC50] != "NA":
                           continue
                       else:
                           for channel in dAC50All[CASID].keys():
                               if not search("Luc_IC50", channel):
                                   if channel != "CASID":
                                       if dAC50All[CASID][channel] != "NA":
                                           lnew = [CASID, "0"]
                                           llineInact.append("\t".join(lnew))
                                           break

                # random active
                nbinactselected = len(llineInact)
                print nbinact, nbinactselected

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


    def writeClass(self):
        from random import shuffle

        for typeAC50 in self.dpresult:
            pclass = self.dpresult[typeAC50] + "actClass.txt"
            print pclass

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

                nbinact = int(100*nbact/(100*self.ratioAct)) - nbact

                nwinact = 0
                flaginact = 0
                # first loop to take active
                for lineChem in llines[1:]:
                    lAC50 = lineChem.strip().split("\t")
                    lnew = [lAC50[0]]
                    for AC50 in lAC50[1:]:
                        if AC50 == "NA":
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


def mergeResults(prin, prout):

    dresult = {}
    dperf = {}
    dperf["Acc"] = []
    dperf["Sp"] = []
    dperf["Se"] = []
    dperf["MCC"] = []

    lprrun = listdir(prin)
    for prrun in lprrun:
        if prrun == "Average":
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


            lML = MCV.keys()
            lcriteria = dperf.keys()
            lset = ["CV", "train", "test"]
            # create the structures
            if not prcell in dresult.keys():
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
    for celltype in dresult.keys():
        for set in dresult[celltype].keys():
            for ML in dresult[celltype][set].keys():
                for criteria in  dresult[celltype][set][ML].keys():

                    AV = round(mean(dresult[celltype][set][ML][criteria]),3)
                    SD = round(std(dresult[celltype][set][ML][criteria]),3)

                    dout[celltype][set][ML][criteria] = [AV, SD]


    # write result
    for celltype in dout.keys():
        pfilout = prout + celltype + ".csv"
        filout = open(pfilout, "w")
        for set in dout[celltype].keys():
            filout.write(str(set) + "\n")
            filout.write("\t" + "\t".join(["M-" + str(c) + "\t" + "SD-" + str(c) for c in dperf.keys()]) + "\n")
            for ML in dout[celltype][set].keys():
                filout.write(ML)
                for criteria in dout[celltype][set][ML].keys():
                    filout.write("\t" + str(dout[celltype][set][ML][criteria][0]) + "\t" + str(dout[celltype][set][ML][criteria][1]))
                filout.write("\n")
        filout.close()


    return 0
