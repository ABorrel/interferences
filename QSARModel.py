import runExternalSoft





class Model:

    def __init__(self, pdesc, pAC50, corval, maxQuantile, splitRatio, nvCV, prresult):

        self.corval = corval
        self.maxQauntile = maxQuantile
        self.prresult = prresult
        self.pdesc = pdesc
        self.pAC50 = pAC50
        self.splitRatio = splitRatio
        self.nbCV = nvCV


    def prepData(self):

        runExternalSoft.prepDataQSAR(self.pdesc, self.pAC50, self.prresult, self.corval, self.maxQauntile, self.splitRatio)

        self.ptrain = self.prresult + "trainSet.csv"
        self.ptest = self.prresult + "testSet.csv"

        return

    def buildQSARReg(self):

        runExternalSoft.QSARReg(self.ptrain, self.ptest, "0", self.prresult, self.nbCV)


