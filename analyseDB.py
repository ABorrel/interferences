from os import path, listdir
import pathFolder
from scipy import stats

import chemical
import runExternalSoft
import toolbox

class Descriptors:

    def __init__(self, prSMI, prDesc, prout, prlog):
        self.prSMI = prSMI
        self.prDesc = prDesc
        self.prlog = prlog
        self.prout = prout

    def computeDesc(self):


        pdesc1D2D = self.prDesc + "tableDesc1D2D"
        self.pdesc1D2D = pdesc1D2D

        prSMIclean = self.prDesc + "SMIclean/"
        pathFolder.createFolder(prSMIclean)
        self.pSMIclean = prSMIclean

        prDescbyCAS = self.prDesc + "DESCbyCAS/"
        pathFolder.createFolder(prDescbyCAS)


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


    def generatePNG(self):

        prPNG = self.prout + "PNG/"
        pathFolder.createFolder(prPNG)
        lnSMIs = listdir(self.pSMIclean)
        for nSMI in lnSMIs:
            runExternalSoft.molconvert(self.pSMIclean + nSMI, prPNG + nSMI[:-3] + "png")



    def setConstantPreproc(self, pAC50, corval, maxQuantile, prAnalysis):

        self.corval = corval
        self.maxQauntile = maxQuantile
        self.pAC50 = pAC50
        self.prAnalysis = prAnalysis


        # output
        paffclean = self.prAnalysis + "IC50Clean.csv"
        pdesc1D2Dclean = self.prAnalysis + "descClean.csv"

        if path.exists(paffclean) and path.exists(pdesc1D2Dclean):
            self.pAC50clean = paffclean
            self.pdesc1D2Dclean = pdesc1D2Dclean
            return 0

        elif path.exists(self.pdesc1D2D) and path.getsize(self.pdesc1D2D) > 10:
            # preproc
            runExternalSoft.dataManager(self.pdesc1D2D, self.pAC50, self.corval, self.maxQauntile, 1, self.prAnalysis)

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


