from os import path, listdir

from mxnet.symbol import mean

import pathFolder
from scipy import stats
from shutil import copyfile

import chemical
import runExternalSoft
import toolbox
from re import search

class Descriptors:

    def __init__(self, prSMI, prDesc, prPNG, prout, prlog):
        self.prSMI = prSMI
        self.prDesc = prDesc
        self.prlog = prlog
        self.prout = prout
        self.prPNG = prPNG



    def loadClassForPred(self, corval, maxQuantile, lCASID = []):

        #Desc 1D and 2D
        self.computeDesc()
        pdescClean = runExternalSoft.dataManager(self.pdesc1D2D, "0", corval, maxQuantile, self.prout)
        self.pdesc1D2Dclean = pdescClean

        if lCASID != []:
            self.reduceMatrixDesc(self.pdesc1D2Dclean, lCASID)


        self.computeFP("All")

        if lCASID != []:
            self.reduceMatrixFP(self.dFP, lCASID)


    def prepareActiveMatrix(self, corval, maxQuantile, pAC50All, prout):

        self.corval = corval
        self.maxQuantile = maxQuantile

        pdescAct = prout + "descActive"
        pAC50Act = prout + "AC50ACactive"

        if path.exists(pdescAct) and path.getsize(pdescAct) > 10 and path.exists(pAC50Act) and path.getsize(pAC50Act) > 10:
            pdescActClean = runExternalSoft.dataManager(pdescAct, pAC50Act, corval, maxQuantile, prout)
            self.pdescCleanActive = pdescActClean
            self.pAC50AllActive = pAC50Act
            return 0

        ddesc = toolbox.loadMatrix(self.pdesc1D2D)
        dAC50All = toolbox.loadMatrix(pAC50All)


        i = 0
        imax = len(ddesc.keys())

        while i < imax:
            casID = ddesc.keys()[i]
            nbNA = 0
            for kAC50 in dAC50All[casID].keys():
                if kAC50 == "CASID" or kAC50 == "Luc_IC50":
                    continue
                else:
                    if dAC50All[casID][kAC50] == "NA":
                        nbNA += 1
            print nbNA, len(dAC50All[casID].keys())
            if nbNA == (len(dAC50All[casID].keys()) -2):
                del dAC50All[casID]
                del ddesc[casID]
                imax = imax -1
            else:
                i += 1

        toolbox.writeMatrix(ddesc, pdescAct)
        toolbox.writeMatrix(dAC50All, pAC50Act)

        pdescActClean = runExternalSoft.dataManager(pdescAct, pAC50Act, corval, maxQuantile, prout)

        self.pdescCleanActive = pdescActClean
        self.pAC50AllActive = pAC50Act

        return 0


    def createActiveSOM(self, sizeMap, prout, pmodelAll):

        import clusteringDB

        # create model
        prall = pathFolder.createFolder(prout + "all/")
        runExternalSoft.drawEnrichSOM(self.pdescCleanActive, self.pAC50AllActive, pmodelAll, prall)

        pModel = prout + "SOMmodel.Rdata"
        if not path.exists(pModel):
            runExternalSoft.generateMainSOM(self.pdescCleanActive, prout, sizeMap, "0")
        else:
            runExternalSoft.generateMainSOM(self.pdescCleanActive, prout, sizeMap, pModel)
        clusteringDB.createSOM(self.pdescCleanActive, self.pAC50AllActive, self.corval, self.maxQuantile, pModel, prout)

        self.prSOMactive = prout


    def extractActivebySOM(self):

        lfolders = listdir(self.prSOMactive)
        dAC50all = toolbox.loadMatrix(self.pAC50AllActive)
        for assay in lfolders:
            pclust = self.prSOMactive + assay + "/SOMClust"

            if not path.exists(pclust):
                continue
            fclust = open(pclust, "r")
            lchemicals = fclust.readlines()
            fclust.close()

            for lineChem in lchemicals[1:]:
                lchemClust = lineChem.strip().replace("\"", "").split(",")
                CAS = lchemClust[0]
                clust = lchemClust[-1]

                if not assay in dAC50all[CAS].keys():
                    continue
                if dAC50all[CAS][assay] != "NA":
                    pclust = pathFolder.createFolder(self.prSOMactive + assay + "/" + str(clust) + "/")
                    copyfile(self.prPNG + CAS + ".png", pclust + CAS + ".png")




    def computeDesc(self):

        pdesc1D2D = self.prDesc + "tableDesc1D2D"
        self.pdesc1D2D = pdesc1D2D

        prSMIclean = self.prDesc + "SMIclean/"
        pathFolder.createFolder(prSMIclean)
        self.prSMIclean = prSMIclean

        prDescbyCAS = self.prDesc + "DESCbyCAS/"
        pathFolder.createFolder(prDescbyCAS)
        self.prDescByCAS = prDescbyCAS


        if path.exists(pdesc1D2D) and path.getsize(pdesc1D2D) > 100:
            return pdesc1D2D
        else:
            ddd
            fdesc1D2D = open(pdesc1D2D, "w")
            ldesc = chemical.getLdesc("1D2D")
            fdesc1D2D.write("CAS\t" + "\t".join(ldesc) + "\n")



        for pSMI in listdir(self.prSMI):
        #for pSMI in ["/home/borrela2/interference/spDataAnalysis/Desc/SMIclean/1212-72-2.smi"]: # to verify for one chem
            cas = pSMI.split("/")[-1].split(".")[0]
            #print cas

            psmiles = self.prSMI + cas + ".smi"
            if path.exists(self.prSMI + cas + ".smi"):
                fsmiles = open(psmiles, "r")
                smiles = fsmiles.readlines()[0].strip()
                fsmiles.close()

                # chemical
                chem = chemical.chemical(cas, smiles)
                chem.prepareChem(prSMIclean)
                chem.compute1D2DDesc(prDescbyCAS)
                err = chem.writeTablesDesc(prDescbyCAS)#
                if err == 1: chem.writelog(self.prlog)
                #Write in the table
                chem.writeDesc(ldesc, fdesc1D2D)
        fdesc1D2D.close()



    def computeFP(self, FPtype):

        # set SMI after cleanning
        prSMIclean = self.prDesc + "SMIclean/"
        pathFolder.createFolder(prSMIclean)
        self.prSMIclean = prSMIclean


        dFP = {}
        i = 1
        for pSMI in listdir(self.prSMI): # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            # for pSMI in ["/home/borrela2/interference/spDataAnalysis/Desc/SMIclean/1212-72-2.smi"]: # to verify for one chem
            cas = pSMI.split("/")[-1].split(".")[0]
            print cas, i, len(listdir(self.prSMI))
            i += 1

            psmiles = self.prSMI + cas + ".smi"
            if path.exists(self.prSMI + cas + ".smi"):
                fsmiles = open(psmiles, "r")
                smiles = fsmiles.readlines()[0].strip()
                fsmiles.close()

                # chemical
                chem = chemical.chemical(cas, smiles)
                chem.prepareChem(prSMIclean)
                error = chem.computeFP(FPtype)

                if error == 1:
                    print "ERROR FP"
                    continue
                else:
                    dFP[cas] = chem.FP

        self.dFP = dFP


    def computeFPMatrix(self, prFP, FPtype, typeMetric):

        from rdkit import DataStructs

        # set FP
        self.prFP = prFP

        # to just load the file
        pfilout = prFP + str(FPtype) + "-" + str(typeMetric)
        if path.exists(pfilout):
            self.pFP = pfilout
            return 0

        self.computeFP(FPtype)

        if FPtype == "pairs" or FPtype == "Torsion" or FPtype == "Morgan":
            if typeMetric != "Dice":
                print "Similarity metric incompatible for ", FPtype, typeMetric
                return 1

        lcas = self.dFP.keys()
        dmetric = {}
        i = 0
        imax = len(self.dFP.keys())
        while i < imax:
            if not lcas[i] in dmetric.keys():
                dmetric[lcas[i]] = {}
            j = i
            while j < imax:
                if not lcas[j] in dmetric[lcas[i]].keys():
                    dmetric[lcas[i]][lcas[j]] = {}

                if typeMetric == 'Tanimoto':
                    dmetric[lcas[i]][lcas[j]][typeMetric] = DataStructs.FingerprintSimilarity(self.dFP[lcas[i]][FPtype],
                                                                                              self.dFP[lcas[j]][FPtype],
                                                                                                  metric=DataStructs.TanimotoSimilarity)
                elif typeMetric == "Dice":
                    dmetric[lcas[i]][lcas[j]][typeMetric] = DataStructs.DiceSimilarity(self.dFP[lcas[i]][FPtype],self.dFP[lcas[j]][FPtype])

                elif typeMetric == "Cosine":
                    dmetric[lcas[i]][lcas[j]][typeMetric] = DataStructs.FingerprintSimilarity(self.dFP[lcas[i]][FPtype],
                                                                                              self.dFP[lcas[j]][FPtype],
                                                                                                  metric=DataStructs.CosineSimilarity)
                elif typeMetric == "Sokal":
                    dmetric[lcas[i]][lcas[j]][typeMetric] = DataStructs.FingerprintSimilarity(self.dFP[lcas[i]][FPtype],
                                                                                              self.dFP[lcas[j]][FPtype],
                                                                                                  metric=DataStructs.SokalSimilarity)
                elif typeMetric == "Russel":
                    dmetric[lcas[i]][lcas[j]][typeMetric] = DataStructs.FingerprintSimilarity(self.dFP[lcas[i]][FPtype],
                                                                                              self.dFP[lcas[j]][FPtype],
                                                                                                  metric=DataStructs.RusselSimilarity)
                elif typeMetric == "RogotGoldberg":
                    dmetric[lcas[i]][lcas[j]][typeMetric] = DataStructs.FingerprintSimilarity(self.dFP[lcas[i]][FPtype],
                                                                                              self.dFP[lcas[j]][FPtype],
                                                                                                  metric=DataStructs.RogotGoldbergSimilarity)
                elif typeMetric == "AllBit":
                    dmetric[lcas[i]][lcas[j]][typeMetric] = DataStructs.FingerprintSimilarity(self.dFP[lcas[i]][FPtype],
                                                                                              self.dFP[lcas[j]][FPtype],
                                                                                                  metric=DataStructs.AllBitSimilarity)
                elif typeMetric == "Kulczynski":
                    dmetric[lcas[i]][lcas[j]][typeMetric] = DataStructs.FingerprintSimilarity(self.dFP[lcas[i]][FPtype],
                                                                                              self.dFP[lcas[j]][FPtype],
                                                                                                  metric=DataStructs.KulczynskiSimilarity)
                elif typeMetric == "McConnaughey":
                    dmetric[lcas[i]][lcas[j]][typeMetric] = DataStructs.FingerprintSimilarity(self.dFP[lcas[i]][FPtype],
                                                                                              self.dFP[lcas[j]][FPtype],
                                                                                                  metric=DataStructs.McConnaugheySimilarity)
                elif typeMetric == "Asymmetric":
                    dmetric[lcas[i]][lcas[j]][typeMetric] = DataStructs.FingerprintSimilarity(self.dFP[lcas[i]][FPtype],
                                                                                              self.dFP[lcas[j]][FPtype],
                                                                                                  metric=DataStructs.AsymmetricSimilarity)
                elif typeMetric == "BraunBlanquet":
                    dmetric[lcas[i]][lcas[j]][typeMetric] = DataStructs.FingerprintSimilarity(self.dFP[lcas[i]][FPtype],
                                                                                              self.dFP[lcas[j]][FPtype],
                                                                                                  metric=DataStructs.BraunBlanquetSimilarity)
                j += 1
            i += 1
        #write matrix of similarity
        filout = open(pfilout, "w")
        filout.write("\t".join(lcas) + "\n")

        i = 0
        imax = len(lcas)
        while i < imax:
            j = 0
            lw = []
            while j < imax:
                try:
                    score = dmetric[lcas[i]][lcas[j]][typeMetric]
                except:
                    score = dmetric[lcas[j]][lcas[i]][typeMetric]
                lw.append(str(score))
                j += 1
            filout.write(lcas[i] + "\t" + "\t".join(lw) + "\n")
            i += 1
        filout.close()
        self.pFP = pfilout



    def generatePNG(self):

        pathFolder.createFolder(self.prPNG)
        lnSMIs = listdir(self.prSMIclean)

        for nSMI in lnSMIs:
            runExternalSoft.molconvert(self.prSMIclean + nSMI, self.prPNG + nSMI[:-3] + "png")




    def setConstantPreproc(self, pAC50, corval, maxQuantile, prAnalysis):

        self.corval = corval
        self.maxQauntile = maxQuantile
        self.pAC50 = pAC50
        self.prAnalysis = prAnalysis


        # output
        paffclean = self.prAnalysis + "IC50Clean.csv"
        pdesc1D2Dclean = self.prAnalysis + "descClean.csv"

        if path.exists(paffclean):
            self.pAC50clean = paffclean
        if path.exists(pdesc1D2Dclean) and pAC50 == "0":
            self.pdesc1D2Dclean = pdesc1D2Dclean
            return 0

        elif path.exists(self.pdesc1D2D) and path.getsize(self.pdesc1D2D) > 10:
            # preproc
            runExternalSoft.dataManager(self.pdesc1D2D, self.pAC50, self.corval, self.maxQauntile, self.prAnalysis)

            if path.exists(paffclean) and path.exists(pdesc1D2Dclean):
                self.pAC50clean = paffclean
                self.pdesc1D2Dclean = pdesc1D2Dclean
                return 0
            else:
                return 1
        else:
            return 1



    def reduceMatrixFP(self, lCASID, pout):

        if path.exists(pout) and path.getsize(pout) > 100:
            return 0

        if not "FP" in self.__dict__:
            self.loadFPmatrix()


        filout = open(pout, "w")
        filout.write("\t".join(lCASID) + "\n")
        for casID in lCASID:
            filout.write(casID)
            i = 0
            imax = len(lCASID)
            while i < imax:
                try:filout.write("\t" + str(self.FP[casID][lCASID[i]]))
                except : filout.write("\t" + str(self.FP[lCASID[i]][casID]))
                i += 1
            filout.write("\n")
        filout.close()
        return 0


    def loadFPmatrix(self):

        dout = {}
        if not "pFP" in self.__dict__:
            print "Error: no FP computed"
            return 1
        else:
            filin = open(self.pFP, "r")
            llinesFP = filin.readlines()
            filin.close()

            lcas = llinesFP[0].strip().split("\t")
            j = 1
            for linesFP in llinesFP[1:]:
                lscore = linesFP.strip().split("\t")
                imax = len(lscore)
                cas = lscore[0]
                i = 1
                dout[cas] = {}
                while i < imax:
                    dout[cas][lcas[i-1]] = lscore[i]
                    i += 1
                print j
                j += 1
            self.FP = dout
        return 0


    def clustering (self, disttype, clustype, agregmethod, optnbclustmethod):

        self.dist = disttype
        self.clusteringMeth = clustype
        self.aggreg = agregmethod
        self.optCluster = optnbclustmethod

        prcluster = self.prAnalysis + "clustering/" + str(disttype) + "_" + str(clustype) + "_" + str(agregmethod) + \
                    "_" + str(optnbclustmethod) + "/"
        pathFolder.createFolder(prcluster)

        pcluster = runExternalSoft.clustering(self.pdesc1D2Dclean, self.pAC50clean, prcluster, self.dist, self.aggreg, self.clusteringMeth, self.optCluster)
        if pcluster == 0:
            return 1
        else:
            self.pcluster = pcluster
            return 0


    def MainSOM(self, sizeMap):


        pModel = self.prAnalysis + "SOMmodel.Rdata"
        runExternalSoft.generateMainSOM(self.pdesc1D2Dclean, self.prAnalysis, sizeMap, pModel)
        if path.exists(pModel):
            return pModel
        return "0"


    def preliminaryAnalysis(self):

        return


    def rankingAC50(self):

        prRank = self.prAnalysis + "ranking/"
        pathFolder.createFolder(prRank)
        fchemAC50 = open(self.pAC50clean, "r")
        lchemac50 = fchemAC50.readlines()
        fchemAC50.close()

        dstock = {}
        lheader = lchemac50[0].strip().split(",")
        i = 1
        dstock["name"] = []
        while i < len(lheader):
            dstock[i] = []
            i += 1


        for chemAC50 in lchemac50[1:]:
            lelem = chemAC50.strip().split(",")
            name = lelem[0]
            dstock["name"].append(name)

            i = 1
            while i < len(lelem):
                dstock[i].append(lelem[i])
                i += 1

        for i in range(1, len(lheader)):
            dstock[i] = toolbox.rankList(dstock[i])


        # generate 3 files ordered differently
        for i in range(1, len(lheader)):
            prank = prRank + str(lheader[i].replace("\"", "")) + "_rank.txt"
            frank = open(prank, "w")
            frank.write(lchemac50[0])


            lisorted = sorted(range(len(dstock[i])), key=lambda k: dstock[i][k])

            for isorted in lisorted:
                j = 1
                frank.write(str(dstock["name"][isorted]))
                while j < len(lheader):
                    frank.write("," + str(dstock[j][isorted]))
                    j += 1
                frank.write("\n")
            frank.close()



def VennCross(cluc, chepg2, chek293, prPNG, prout, verbose = 0):


    if not "dresponse" in chepg2.__dict__:
        chepg2.responseCurves(drawn=0)

    if not "dresponse" in chek293.__dict__:
        chek293.responseCurves(drawn=0)

    lsample = chepg2.dresponse[chepg2.dresponse.keys()[0]].keys()
    lcolor = ["blue", "blue_n", "red", "red_n", "green", "green_n"]

    for sample in lsample:
        prsub = pathFolder.createFolder(prout + str(sample) + "/")

        for CASID in chepg2.dresponse.keys():
            if chepg2.dresponse[CASID][sample]["AC50"] == "NA" or chek293.dresponse[CASID][sample]["AC50"] == "NA":
                continue
            if float(chepg2.dresponse[CASID][sample]["CURVE_CLASS2"]) >= 4 or float(chek293.dresponse[CASID][sample]["CURVE_CLASS2"]) >= 4:
                continue

            if path.exists(prPNG + CASID + ".png"):
                copyfile(prPNG + CASID + ".png", prsub + CASID + ".png")

    #for color
    for color in lcolor:
        prsub = pathFolder.createFolder(prout + str(color) + "/")

        if verbose == 1: print len(chepg2.dresponse.keys())

        for CASID in chepg2.dresponse.keys():
            if chepg2.dresponse[CASID]["cell_" + color]["AC50"] == "NA" or chek293.dresponse[CASID]["cell_" + color]["AC50"] == "NA"\
                    or chepg2.dresponse[CASID]["med_" + color]["AC50"] == "NA" or chek293.dresponse[CASID]["med_" + color]["AC50"] == "NA":
                continue
            if float(chepg2.dresponse[CASID]["cell_" + color]["CURVE_CLASS2"]) >= 4 or float(chek293.dresponse[CASID]["cell_" + color]["CURVE_CLASS2"]) >= 4 \
                    or float(chepg2.dresponse[CASID]["med_" + color]["CURVE_CLASS2"]) >= 4 or float(chek293.dresponse[CASID]["med_" + color]["CURVE_CLASS2"]) >= 4:
                continue

            if path.exists(prPNG + CASID + ".png"):
                copyfile(prPNG + CASID + ".png", prsub + CASID + ".png")

    runExternalSoft.crossVenn(cluc.pAC50, chepg2.pAC50, chek293.pAC50, prout)



def PCACross(pdesc, pAC50_hepg2, pAC50_hek293, corval, maxQuantile, prCrossPCA):


    # output
    pdesc1D2Dclean = prCrossPCA + "descClean.csv"

    if not path.exists(pdesc1D2Dclean):

        if path.exists(pdesc) and path.getsize(pdesc) > 10:
            # preproc
            runExternalSoft.dataManager(pdesc, 0, corval, maxQuantile, prCrossPCA)
        else:
            print "Error ->", pdesc

    runExternalSoft.drawPCACross(pdesc1D2Dclean, pAC50_hepg2, pAC50_hek293, prCrossPCA)


