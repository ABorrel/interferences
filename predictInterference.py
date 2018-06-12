import pathFolder
import chemical
import runExternalSoft
import toolbox

from os import path, listdir




class predictor:


    def __init__(self, cDB, prclusters, lclust, prout, distMeth, aggMeth, verbose = 1):

        self.prout = pathFolder.createFolder(prout)
        self.prcluster = prclusters
        self.cDB = cDB
        self.verbose = verbose
        self.distMeth = distMeth
        self.aggMeth = aggMeth

        #load model
        self.loadCluster(lclust)




    def loadCluster(self, lclustConsidered =[]):

        lchanel = listdir(self.prcluster)
        dclustering = {}
        dCASCluster = {}
        for chanel in lchanel:
            dclustering[chanel] = {}
            lcell = listdir(self.prcluster + "/" + chanel + "/")

            for cell in lcell:
                dclustering[chanel][cell] = {}
                lclustering = listdir(self.prcluster + "/" + chanel + "/" + str(cell) + "/")

                if lclustConsidered != []:
                    lclustering = lclustConsidered

                for clustering in lclustering:

                    dclustering[chanel][cell][clustering]={}

                    # file in
                    pcluster = self.prcluster + chanel + "/" + str(cell) + "/" + str(
                        clustering) + "/enrichment_cluster.csv"
                    pCAScluster = self.prcluster + chanel + "/" + str(cell) + "/" + str(clustering) + "/" + str(
                        cell) + "_" + str(chanel) + "_cluster.csv"

                    # manage clustering computed
                    fcluster = open(pcluster, "r")
                    lclusters = fcluster.readlines()
                    fcluster.close()

                    lheader = lclusters[0].strip().split(",")
                    nbcluster = 0
                    for cluster in lclusters[1:]:
                        enrichCluster = cluster.strip().split(",")
                        nameclust = str(enrichCluster[1])
                        print nameclust
                        if int(nameclust) > nbcluster:
                            nbcluster = nameclust

                        dclust = {}
                        i = 1
                        imax = len(lheader)
                        while i < imax:
                            dclust[lheader[i].replace("\"", "")] = enrichCluster[i]
                            i += 1
                        dclustering[chanel][cell][clustering][nameclust] = dclust
                    dclustering[chanel][cell][clustering]["NB clusters"] = nbcluster
                    dclustering[chanel][cell][clustering]["files"] = [pcluster, pCAScluster]

                    #by compound
                    fCAScluster = open(pCAScluster, "r")
                    lCAScluster = fCAScluster.readlines()
                    fCAScluster.close()

                    for CASClust in lCAScluster[1:]:
                        elem = CASClust.strip().split(",")

                        CASID = str(elem[0].replace("\"", ""))
                        clust = str(elem[-1].replace("\"", ""))

                        if not CASID in dCASCluster:
                            dCASCluster[CASID] = {}
                        if not chanel in dCASCluster[CASID].keys():
                            dCASCluster[CASID][chanel] = {}
                        if not cell in dCASCluster[CASID][chanel].keys():
                            dCASCluster[CASID][chanel][cell] = {}

                        dCASCluster[CASID][chanel][cell][clustering] = clust

        self.dcluster = dclustering
        self.ChemClust = dCASCluster

        if self.verbose == 1:
            print self.ChemClust.keys()[1]
            print self.ChemClust[self.ChemClust.keys()[1]]
            print self.dcluster.keys()[1]
            print self.dcluster[self.dcluster.keys()[1]]

        return 0


    def predictlpSMI(self, lpsmi):

        for psmi in lpsmi:
            nameChemical = psmi.split("/")[-1][0:-4]
            filin = open(psmi, "r")
            smiles = filin.readlines()[0]
            smiles.strip()
            print smiles
            filin.close()

            self.predictSMI(nameChemical, smiles)




    def predictSMI(self, nameChemical, smiles):

        dpred = {}
        prresult = pathFolder.createFolder(self.prout + nameChemical + "/")
        chem = chemical.chemical(nameChemical, smiles)
        chem.prepareChem(prresult)
        chem.compute1D2DDesc(prresult)
        chem.writeTablesDesc(prresult)
        chem.computeFP(typeFP="All")

        for channel in self.dcluster:
            dpred[channel] = {}
            for cell in self.dcluster[channel].keys():
                dpred[channel][cell] = {}
                for typeDesc in self.dcluster[channel][cell].keys():
                    print channel, cell, typeDesc
                    if typeDesc == "Desc":
                        enrichment = runExternalSoft.findCluster(self.cDB.pdesc1D2Dclean, chem.pdesc,
                                                    self.dcluster[channel][cell][typeDesc]["files"][0],
                                                    self.dcluster[channel][cell][typeDesc]["files"][1],
                                                    self.distMeth, self.aggMeth)

                    else:
                        # generate FP
                        typeFP = typeDesc.split("-")[0][2:]
                        metric = typeDesc.split("-")[-1]
                        print typeFP, metric
                        dFP = {}
                        for CASID in self.cDB.dFP.keys():
                            dFP[CASID] = float(toolbox.computeSimilarityFP(self.cDB.dFP[CASID], chem.FP[typeFP], metric))
                        maxSim = max(dFP.values())
                        i = 0
                        imax = len(dFP.keys())
                        lCAS = dFP.keys()
                        while i < imax:
                            if float(dFP[lCAS[i]] == maxSim):
                                CASclose = lCAS[i]
                            i += 1

                        clusterfound = self.ChemClust[CASclose][channel][cell][typeFP]["Cluster"]
                        enrichment = self.dcluster[channel][cell][typeDesc][clusterfound]['Enrichment']
                    dpred[channel][cell][typeDesc] = enrichment


        print dpred


        pfilout = self.prout + "pred.csv"
        filout = open(pfilout, "w")
        filout.write("Clustering\tEnrichment\n")
        for k in dpred.keys():
            filout.write(str(k) + "\t" + str(dpred[k]) + "\n")
        filout.close()

        #runExternalSoft.generateResultFig(pfilout)