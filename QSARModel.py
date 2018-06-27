import runExternalSoft
import pathFolder

from os import path




class Model:

    def __init__(self, pdesc, pAC50, typeQSAR, corval, maxQuantile, splitRatio, nbCV, ratioAct, lchannels, prresult):

        self.corval = corval
        self.maxQauntile = maxQuantile
        self.prresult = prresult
        self.pdesc = pdesc
        self.pAC50 = pAC50
        self.splitRatio = splitRatio
        self.ratioAct = ratioAct
        self.nbCV = nbCV
        self.typeQSAR = typeQSAR
        self.lchannel = lchannels


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
                self.writeClass()
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

