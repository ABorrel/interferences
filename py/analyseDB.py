from os import path, listdir

import pathFolder
from scipy import stats
from shutil import copyfile
from copy import deepcopy

#import chemical
import runExternalSoft
import toolbox
from re import search



def loadAllOperaDesc(pOperaDesc):
    dDTX = toolbox.loadMatrix(pOperaDesc, ',')
    dCAS = {}
    for DTXID in list(dDTX.keys()):
        CASID = dDTX[DTXID]["CASRN"]
        dCAS[CASID] = {}
        for desc in chemical.LOPERA:
            if dDTX[DTXID][desc] == "NaN":
                dDTX[DTXID][desc] = "NA"
            dCAS[CASID][desc] = dDTX[DTXID][desc]
    return dCAS




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
        lpdescClean = runExternalSoft.dataManager(self.pdesc1D2D, "0", corval, maxQuantile, self.prout)

        self.pdesc1D2Dclean = lpdescClean[0]

        if lCASID != []:
            self.reduceMatrixDesc(self.pdesc1D2Dclean, lCASID)


        self.computeFP("All")

        if lCASID != []:
            self.reduceMatrixFP(self.dFP, lCASID)


    def prepareActiveMatrix(self, corval, maxQuantile, NBNA, pAC50All, prout, luciferase=0):

        self.corval = corval
        self.maxQuantile = maxQuantile

        pdescAct = prout + "descActive"
        pAC50Act = prout + "AC50Active"

        if path.exists(pdescAct) and path.getsize(pdescAct) > 10 and path.exists(pAC50Act) and path.getsize(pAC50Act) > 10:
            lpdescActClean = runExternalSoft.dataManager(pdescAct, pAC50Act, corval, maxQuantile, NBNA, prout)
            self.pdescCleanActive = lpdescActClean[0]
            self.pAC50AllActive = lpdescActClean[1]
            return [self.pdescCleanActive, self.pAC50AllActive]

        ddesc = toolbox.loadMatrix(self.pdesc1D2D)
        dAC50All = toolbox.loadMatrix(pAC50All)


        if luciferase == 0:
            i = 0
            imax = len(list(ddesc.keys()))

            while i < imax:
                casID = list(dAC50All.keys())[i]
                nbNA = 0
                for kAC50 in list(dAC50All[casID].keys()):
                    if kAC50 == "CASID" or kAC50 == "Luc_IC50":# not considered luciferase
                        continue
                    else:
                        if dAC50All[casID][kAC50] == "NA":
                            nbNA += 1
                if nbNA == (len(list(dAC50All[casID].keys())) -2):
                    del dAC50All[casID]
                    try:
                        del ddesc[casID]
                    except:
                        pass
                    imax = imax -1
                else:
                    i += 1

            toolbox.writeMatrix(ddesc, pdescAct)
            toolbox.writeMatrix(dAC50All, pAC50Act)

            lpdescActClean = runExternalSoft.dataManager(pdescAct, pAC50Act, corval, maxQuantile, NBNA, prout)

            self.pdescCleanActive = lpdescActClean[0]
            self.pAC50AllActive = lpdescActClean[1]

            return [self.pdescCleanActive, self.pAC50AllActive]


        else:
            i = 0
            imax = len(list(dAC50All.keys()))

            while i < imax:
                casID = list(dAC50All.keys())[i]
                if not casID in list(ddesc.keys()):
                    del dAC50All[casID]
                    imax = imax - 1
                    i = i - 1
                    continue
                for kAC50 in list(dAC50All[casID].keys()):
                    if kAC50 != "Luc_IC50" and kAC50 != "CASID":  # not considered luciferase
                        del dAC50All[casID][kAC50]
                    else:
                        if dAC50All[casID][kAC50] == "NA":
                            del dAC50All[casID]
                            try:del ddesc[casID]
                            except: pass
                            imax = imax - 1
                            i = i -1
                            break
                i += 1

            toolbox.writeMatrix(ddesc, pdescAct)
            toolbox.writeMatrix(dAC50All, pAC50Act)

            lpdescActClean = runExternalSoft.dataManager(pdescAct, pAC50Act, corval, maxQuantile, NBNA, prout)

            self.pdescCleanActive = lpdescActClean[0]
            self.pAC50AllActive = lpdescActClean[1]

            return [self.pdescCleanActive, self.pAC50AllActive]




    def createActiveSOM(self, sizeMap, prout, pmodelAll):

        import clusteringDB

        # create model
        prall = pathFolder.createFolder(prout + "global/")
        runExternalSoft.drawEnrichSOM(self.pdescCleanActive, self.pAC50AllActive, pmodelAll, prall)
        self.extractActivebySOM(prall)
        # to better calibrate => change the max in the R script

        pModel = prout + "SOMmodel.Rdata"
        if not path.exists(pModel):
            runExternalSoft.drawEnrichSOM(self.pdescCleanActive, self.pAC50AllActive, "0", prout)
        else:
            runExternalSoft.drawEnrichSOM(self.pdescCleanActive, self.pAC50AllActive, pModel, prout)

        self.prSOMactive = prout


    def extractActivebySOM(self, prin = ""):

        if prin == "":
            prin = self.prSOMactive

        lfolders = listdir(prin)

        dAC50all = toolbox.loadMatrix(self.pAC50AllActive, sep=",")
        for assay in lfolders:
            pclust = prin + assay + "/SOMClust"

            if not path.exists(pclust):
                continue
            fclust = open(pclust, "r")
            lchemicals = fclust.readlines()
            fclust.close()

            for lineChem in lchemicals[1:]:
                lchemClust = lineChem.strip().replace("\"", "").split(",")
                CAS = lchemClust[0]
                clust = lchemClust[-1]
                if CAS == "NA":
                    continue

                if assay in list(dAC50all[CAS].keys()):
                    if dAC50all[CAS][assay] != "NA":
                        pclust = pathFolder.createFolder(prin + assay + "/" + str(clust) + "/")
                        copyfile(self.prPNG + CAS + ".png", pclust + CAS + ".png")
                    continue
                elif assay == "red" or assay == "green" or assay == "blue" or assay == "allcolors":
                    lassays = ["hepg2_cell_X_n", "hepg2_med_X_n", "hek293_med_X_n", "hek293_cell_X_n"]
                    if assay == "allcolors":
                        lassayout = []
                        lassayout = lassayout + [i.replace("X", "blue") for i in lassays] + [i.replace("X", "green") for i in lassays] + [i.replace("X", "red") for i in lassays]
                        lassays = lassayout
                    else:
                        lassays = [i.replace("X", assay) for i in lassays]
                    for ass in lassays:
                        if dAC50all[CAS][ass] != "NA":
                            pclust = pathFolder.createFolder(prin + assay + "/" + str(clust) + "/")
                            copyfile(self.prPNG + CAS + ".png", pclust + CAS + ".png")
                            break
                elif search("hepg2", assay) or search("hek293", assay):
                    lassays = ["X_cell_Y_n", "X_med_Y_n"]
                    lass = assay.split("_")
                    lassays = [i.replace("X", lass[0]).replace("Y", lass[1]) for i in lassays]

                    for ass in lassays:
                        if dAC50all[CAS][ass] != "NA":
                            pclust = pathFolder.createFolder(prin + assay + "/" + str(clust) + "/")
                            copyfile(self.prPNG + CAS + ".png", pclust + CAS + ".png")
                            break




    def computeDesc(self, opera=0, RDkitPhysico=1, pOperaDesc=""):

        pdesc1D2D = self.prDesc + "tableDesc1D2D"
        self.pdesc1D2D = pdesc1D2D

        if opera == 1:
            pdesc1D2D = pdesc1D2D + "Opera"
            self.pdesc1D2D = pdesc1D2D

        if RDkitPhysico == 0:
            pdesc1D2D = pdesc1D2D + "NoRDKITPhyChem"
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
            plog = self.prDesc + "log.log"
            flog = open(plog, "w")
            print((pdesc1D2D, "No found"))
            fdesc1D2D = open(pdesc1D2D, "w")
            if opera == 0:
                ldesc = chemical.getLdesc("1D2D", RDkitPhysico)
                fdesc1D2D.write("CAS\t" + "\t".join(ldesc) + "\n")
            else:
                ldesc = chemical.getLdesc("1D2D", RDkitPhysico)
                ldesc = ldesc + chemical.getLdesc("Opera", RDkitPhysico)
                doperaDesc = loadAllOperaDesc(pOperaDesc)
                fdesc1D2D.write("CAS\t" + "\t".join(ldesc) + "\n")


        for pSMI in listdir(self.prSMI):
        #for pSMI in ["/home/borrela2/interference/spDataAnalysis/Desc/SMIclean/1212-72-2.smi"]: # to verify for one chem
            cas = pSMI.split("/")[-1].split(".")[0]

            psmiles = self.prSMI + cas + ".smi"
            if path.exists(self.prSMI + cas + ".smi"):
                fsmiles = open(psmiles, "r")
                smiles = fsmiles.readlines()[0].strip()
                fsmiles.close()

                # chemical
                chem = chemical.chemical(cas, smiles)
                chem.prepareChem(prSMIclean)
                chem.compute1D2DDesc(prDescbyCAS)
                if opera == 1:
                    chem.loadOperaDesc(doperaDesc, flog)
                err = chem.writeTablesDescCAS(prDescbyCAS)#
                if err == 1: chem.writelog(self.prlog)
                #Write in the table
                chem.writeDesc(ldesc, fdesc1D2D)
        fdesc1D2D.close()
        flog.close()



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
            print((cas, i, len(listdir(self.prSMI))))
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
                    print ("ERROR FP")
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
                print(("Similarity metric incompatible for ", FPtype, typeMetric))
                return 1

        lcas = list(self.dFP.keys())
        dmetric = {}
        i = 0
        imax = len(list(self.dFP.keys()))
        while i < imax:
            if not lcas[i] in list(dmetric.keys()):
                dmetric[lcas[i]] = {}
            j = i
            while j < imax:
                if not lcas[j] in list(dmetric[lcas[i]].keys()):
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




    def setConstantPreproc(self, pAC50, corval, maxQuantile, nbNA, prAnalysis):

        self.corval = corval
        self.maxQauntile = maxQuantile
        self.pAC50 = pAC50
        self.prAnalysis = prAnalysis


        # output
        paffclean = self.prAnalysis + "AC50Clean.csv"
        pdesc1D2Dclean = self.prAnalysis + "descClean.csv"

        if path.exists(paffclean):
            self.pAC50clean = paffclean
        if path.exists(pdesc1D2Dclean):
            self.pdesc1D2Dclean = pdesc1D2Dclean
            return 0

        elif path.exists(self.pdesc1D2D) and path.getsize(self.pdesc1D2D) > 10:
            # preproc
            runExternalSoft.dataManager(self.pdesc1D2D, self.pAC50, self.corval, self.maxQauntile, nbNA, self.prAnalysis)

            if path.exists(paffclean):
                self.pAC50clean = paffclean
            if path.exists(pdesc1D2Dclean):
                self.pdesc1D2Dclean = pdesc1D2Dclean
        return 0



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
            print ("Error: no FP computed")
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
                print (j)
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
        if path.exists(pModel):
            return pModel
        else:
            runExternalSoft.generateMainSOM(self.pdesc1D2Dclean, self.prAnalysis, sizeMap, "0")
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


            lisorted = sorted(list(range(len(dstock[i]))), key=lambda k: dstock[i][k])

            for isorted in lisorted:
                j = 1
                frank.write(str(dstock["name"][isorted]))
                while j < len(lheader):
                    frank.write("," + str(dstock[j][isorted]))
                    j += 1
                frank.write("\n")
            frank.close()



    def Ttest(self, pAC50All, presult):

        dAC50 = toolbox.loadMatrix(pAC50All, sep="\t")
        ddesc = toolbox.loadMatrix(self.pdesc1D2Dclean, sep=",")

        print((list(ddesc.keys())[:20]))
        print((list(dAC50.keys())[:20]))

        runExternalSoft.TtestDesc(self.pdesc1D2Dclean, pAC50All, presult)




def VennCross(cluc, chepg2, chek293, prPNG, prout, verbose = 0):


    lsample = list(chepg2.dAC50[list(chepg2.dAC50.keys())[0]].keys())

    if "CAS" in lsample:
        del lsample[lsample.index("CAS")]
    elif "CASID" in lsample:
        del lsample[lsample.index("CASID")]


    lcolor = ["blue_n", "red_n", "green_n"]

    for sample in lsample:
        prsub = pathFolder.createFolder(prout + str(sample) + "/")

        for CASID in list(chepg2.dAC50.keys()):
            if chepg2.dAC50[CASID][sample] == "NA" or chek293.dAC50[CASID][sample] == "NA": # have to look if open good
                continue

            if path.exists(prPNG + CASID + ".png"):
                copyfile(prPNG + CASID + ".png", prsub + CASID + ".png")

    #for color
    for color in lcolor:
        prsub = pathFolder.createFolder(prout + str(color) + "/")

        if verbose == 1: print(len(list(chepg2.dAC50.keys())))

        for CASID in list(chepg2.dAC50.keys()):
            if chepg2.dAC50[CASID]["cell_" + color] == "NA" or chek293.dAC50[CASID]["cell_" + color] == "NA"\
                    or chepg2.dAC50[CASID]["med_" + color] == "NA" or chek293.dAC50[CASID]["med_" + color] == "NA":
                continue

            if path.exists(prPNG + CASID + ".png"):
                copyfile(prPNG + CASID + ".png", prsub + CASID + ".png")



    # all active
    prsub = pathFolder.createFolder(prout + "all/")

    for CASID in list(chepg2.dAC50.keys()):
        flag = 1
        for color in lcolor:
            if chepg2.dAC50[CASID]["med_" + color] != "NA":
                continue

            if chek293.dAC50[CASID]["med_" + color] != "NA":
                continue

            if chepg2.dAC50[CASID]["cell_" + color] != "NA":
                continue

            if chek293.dAC50[CASID]["cell_" + color] != "NA":
                continue

            flag = 0
            break

        if flag == 1:
            if path.exists(prPNG + CASID + ".png"):
                copyfile(prPNG + CASID + ".png", prsub + CASID + ".png")


    runExternalSoft.crossVenn(cluc.pAC50, chepg2.pAC50, chek293.pAC50, prout)



def PCACross(pdesc, pAC50_hepg2, pAC50_hek293, nbNA, corval, maxQuantile, prCrossPCA):


    # output
    pdesc1D2Dclean = prCrossPCA + "descClean.csv"

    if not path.exists(pdesc1D2Dclean):

        if path.exists(pdesc) and path.getsize(pdesc) > 10:
            # preproc
            runExternalSoft.dataManager(pdesc, 0, corval, maxQuantile, nbNA, prCrossPCA)
        else:
            print(("Error ->", pdesc))

    runExternalSoft.drawPCACross(pdesc1D2Dclean, pAC50_hepg2, pAC50_hek293, prCrossPCA)



def VennTox21PubMed(pallIC50, prMain, prout):

    # create matrix smiles for Tox21
    ptox21Smiles = prout + "tox21SMI.csv"
    if not path.exists(ptox21Smiles):
        prSMI = prMain + "Desc/SMIclean/"
        ftox21Smiles = open(ptox21Smiles, "w")
        ftox21Smiles.write("CAS\tSMILES\n")
        lCASsmi = listdir(prSMI)
        for CASsmi in lCASsmi:
            fSMILES = open(prSMI + CASsmi, "r")
            linesSmiles = fSMILES.readlines()
            SMILES = linesSmiles[0].strip()
            fSMILES.close()
            CAS = CASsmi.split(".")[0]
            ftox21Smiles.write(CAS + "\t" + SMILES + "\n")
        ftox21Smiles.close()

    pPubChemSmiles = prout + "PubChemSMI.csv"
    if not path.exists(pPubChemSmiles):
        prSMIPubchem = prMain + "testing/SMI/"
        fPubChemSmiles = open(pPubChemSmiles, "w")
        fPubChemSmiles.write("ID\tSMILES\n")
        lIDsmi = listdir(prSMIPubchem)
        for IDsmi in lIDsmi:
            fSMILES = open(prSMIPubchem + IDsmi, "r")
            linesSmiles = fSMILES.readlines()
            SMILES = linesSmiles[0].strip()
            fSMILES.close()
            ID = IDsmi.split(".")[0]
            fPubChemSmiles.write(ID + "\t" + SMILES + "\n")
        fPubChemSmiles.close()


    for PUBMEDfold in listdir(prMain + "testing/"):
        if not search("old", PUBMEDfold) and PUBMEDfold != "data" and PUBMEDfold != "DESC" and PUBMEDfold != "QSARmodel" and PUBMEDfold != "SMI":
            prvenn = pathFolder.createFolder(prout + PUBMEDfold + "/")
            pIC50PubChem = prMain + "testing/" + PUBMEDfold + "/tableSmi.csv"
            runExternalSoft.vennPlotPubChem(ptox21Smiles, pallIC50, pPubChemSmiles, pIC50PubChem, prvenn)





