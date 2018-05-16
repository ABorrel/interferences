import runExternalSoft
from os import path




class Model:

    def __init__(self, pdesc, pAC50, typeQSAR, corval, maxQuantile, splitRatio, nbCV, prresult):

        self.corval = corval
        self.maxQauntile = maxQuantile
        self.prresult = prresult
        self.pdesc = pdesc
        self.pAC50 = pAC50
        self.splitRatio = splitRatio
        self.nbCV = nbCV
        self.typeQSAR = typeQSAR


    def prepData(self):

        if self.typeQSAR == "Reg":
            runExternalSoft.prepDataQSAR(self.pdesc, self.pAC50, self.prresult, self.corval, self.maxQauntile, self.splitRatio)

        else:
            self.writeClass()
            runExternalSoft.prepDataQSAR(self.pdesc, self.pAC50, self.prresult, self.corval, self.maxQauntile,
                                         self.splitRatio)

        self.ptrain = self.prresult + "trainSet.csv"
        self.ptest = self.prresult + "testSet.csv"

        return

    def buildQSARReg(self):

        runExternalSoft.QSARReg(self.ptrain, self.ptest, "0", self.prresult, self.nbCV)



    def buildQSARClass(self):

        runExternalSoft.QSARClass(self.ptrain, self.ptest, self.prresult, self.nbCV)


    def writeClass(self):

        pclass = self.prresult + "actClass.txt"
        print pclass

        if path.exists(pclass) and path.getsize(pclass) > 10:
            self.pAC50 = pclass

        else:
            filin = open(self.pAC50, "r")
            llines = filin.readlines()
            filin.close()

            filout = open(pclass, "w")
            filout.write(llines[0])

            for lineChem in llines[1:]:
                lAC50 = lineChem.strip().split("\t")
                lnew = [lAC50[0]]
                for AC50 in lAC50[1:]:
                    if AC50 == "NA":
                        lnew.append("0")
                    else:
                        lnew.append("1")
                filout.write("\t".join(lnew) + "\n")
            filout.close()
            self.pAC50 = pclass

