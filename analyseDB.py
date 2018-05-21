from os import path, listdir
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

    def computeDesc(self):


        pdesc1D2D = self.prDesc + "tableDesc1D2D"
        self.pdesc1D2D = pdesc1D2D

        prSMIclean = self.prDesc + "SMIclean/"
        pathFolder.createFolder(prSMIclean)
        self.pSMIclean = prSMIclean

        prDescbyCAS = self.prDesc + "DESCbyCAS/"
        pathFolder.createFolder(prDescbyCAS)
        self.prDescByCAS = prDescbyCAS


        if path.exists(pdesc1D2D) and path.getsize(pdesc1D2D) > 100:
            return pdesc1D2D
        else:
            fdesc1D2D = open(pdesc1D2D, "w")
            ldesc = chemical.getLdesc("1D2D")
            fdesc1D2D.write("CAS\t" + "\t".join(ldesc) + "\n")



        for pSMI in listdir(self.prSMI):
        #for pSMI in ["/home/borrela2/interference/spDataAnalysis/Desc/SMIclean/1212-72-2.smi"]: # to verify for one chem
            cas = pSMI.split("/")[-1].split(".")[0]
            print cas

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



    def computeFP(self, prFP):

        # set SMI after cleanning
        prSMIclean = self.prDesc + "SMIclean/"
        pathFolder.createFolder(prSMIclean)
        self.pSMIclean = prSMIclean

        #set FP
        self.prFP = prFP


        dFP = {}
        for pSMI in listdir(self.prSMI):
        #for pSMI in ["/home/borrela2/interference/spDataAnalysis/Desc/SMIclean/1212-72-2.smi"]: # to verify for one chem
            cas = pSMI.split("/")[-1].split(".")[0]
            print cas

            psmiles = self.prSMI + cas + ".smi"
            if path.exists(self.prSMI + cas + ".smi"):
                fsmiles = open(psmiles, "r")
                smiles = fsmiles.readlines()[0].strip()
                fsmiles.close()

                # chemical
                chem = chemical.chemical(cas, smiles)
                chem.prepareChem(prSMIclean)
                chem.computeFP()
                dFP[cas] = chem.FP

        # create matrix of of FP score similarity
        lmetric = to fisnish
        lcas = dFP.keys()
        imax = len(lcas)




                chem.computeFP(self.prFPbyCAS)

            # control if exit already
            if not path.exists(self.prDescByCAS + cas + ".txt"):

                psmiles = self.prSMI + cas + ".smi"
                if path.exists(self.prSMI + cas + ".smi"):
                    fsmiles = open(psmiles, "r")
                    smiles = fsmiles.readlines()[0].strip()
                    fsmiles.close()

                    # chemical
                    chem = chemical.chemical(cas, smiles)
                    chem.prepareChem(prSMIclean)
                    chem.computeFP(self.prFPbyCAS)
                    fff



                    chem.compute1D2DDesc(self.prDescByCAS)
                    err = chem.writeTablesDesc(self.prDescByCAS)#
                    if err == 1: chem.writelog(self.prlog)


        return


    def generatePNG(self):

        pathFolder.createFolder(self.prPNG)
        lnSMIs = listdir(self.pSMIclean)

        for nSMI in lnSMIs:
            runExternalSoft.molconvert(self.pSMIclean + nSMI, self.prPNG + nSMI[:-3] + "png")




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

        runExternalSoft.generateMainSOM(self.pdesc1D2Dclean, self.prAnalysis, sizeMap)
        pModel = self.prAnalysis + "SOMmodel.Rdata"
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



def VennCross(cluc, chepg2, chek293, prPNG, prout):


    if not "dresponse" in chepg2.__dict__:
        chepg2.responseCurves(drawn=0)

    if not "dresponse" in chek293.__dict__:
        chek293.responseCurves(drawn=0)

    lsample = chepg2.dresponse[chepg2.dresponse.keys()[0]].keys()
    lcolor = ["blue", "blue_n", "red", "green"]

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

        print len(chepg2.dresponse.keys()), "ddd"

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


