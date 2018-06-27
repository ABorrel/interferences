import pathFolder
import chemical
import runExternalSoft
import toolbox

from os import path, listdir
from re import search
from numpy import mean, std



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
            smiles = toolbox.loadSMILES(psmi)
            self.predictSMI(nameChemical, smiles, plot=1)





    def predictSMI(self, nameChemical, smiles, plot = 0, verbose = 0):

        dpred = {}
        prresult = pathFolder.createFolder(self.prout + nameChemical + "/")

        ppred = prresult + "pred"
        if path.exists(ppred):
            dpred = toolbox.loadMatrix(ppred)
            return dpred


        chem = chemical.chemical(nameChemical, smiles)
        chem.prepareChem(prresult)
        chem.compute1D2DDesc(prresult)
        chem.writeTablesDesc(prresult)
        chem.computeFP(typeFP="All")

        for channel in self.dcluster:
            for cell in self.dcluster[channel].keys():
                kpred = str(cell) + "_" + str(channel)
                dpred[kpred] = {}
                for typeDesc in self.dcluster[channel][cell].keys():
                    if verbose == 1:
                        print channel, cell, typeDesc
                        print self.dcluster[channel][cell].keys()
                    if search("Desc", typeDesc):
                        distMeth = typeDesc.split("-")[1]
                        aggMeth = typeDesc.split("-")[2]

                        enrichment = runExternalSoft.findCluster(self.cDB.pdesc1D2Dclean, chem.pdesc,
                                                    self.dcluster[channel][cell][typeDesc]["files"][0],
                                                    self.dcluster[channel][cell][typeDesc]["files"][1],
                                                    distMeth, aggMeth)

                    else:
                        # generate FP
                        typeFP = typeDesc.split("-")[0]
                        metric = typeDesc.split("-")[-1].split("_")[0]
                        metricAgg = typeDesc.split("-")[-1]
                        if verbose == 1: print typeFP, metric
                        dFP = {}
                        for CASID in self.cDB.dFP.keys():
                            if verbose == 1:
                                print self.cDB.dFP[CASID]
                                print chem.FP[typeFP]
                                print metric
                            dFP[CASID] = float(toolbox.computeSimilarityFP(self.cDB.dFP[CASID][typeFP], chem.FP[typeFP], metric))
                        maxSim = max(dFP.values())
                        i = 0
                        imax = len(dFP.keys())
                        lCAS = dFP.keys()
                        while i < imax:
                            if float(dFP[lCAS[i]] == maxSim):
                                CASclose = lCAS[i]
                            i += 1
                        if verbose == 1:
                            print CASclose
                            print channel, cell
                            print self.ChemClust[CASclose][channel][cell]

                        clusterfound = self.ChemClust[CASclose][channel][cell][str(typeFP) + "-" + str(metricAgg)]
                        enrichment = self.dcluster[channel][cell][typeDesc][clusterfound]['Enrichment']
                    dpred[kpred][typeDesc] = enrichment

        if plot == 1:
            self.writeResultBySMI(dpred, prresult)

        return dpred




    def writeResultBySMI(self, dpred, prresult):

        lheader = dpred[dpred.keys()[0]].keys()

        pfilout = prresult + "pred"
        filout = open(pfilout, "w")
        filout.write("Interferences" + "\t" + "\t".join(lheader) + "\n")
        for k in dpred.keys():
            filout.write(str(k))
            for h in lheader:
                try:
                    filout.write("\t" + str(dpred[k][h]))
                except:
                    filout.write("\tNA")
            filout.write("\n")
        filout.close()
        runExternalSoft.generateCardResult(pfilout)


    def summarizePredictor(self):

        psum = self.prout + "summary.txt"
        fsum = open(psum, "w")
        fsum.write("Channel\tCell\tDesc\tMSize\tSDsize\n")

        for channel in self.dcluster.keys():
            for cell in self.dcluster[channel].keys():
                for desc in self.dcluster[channel][cell].keys():
                    #print self.dcluster[channel][cell][desc]
                    lsize = []
                    for clust in self.dcluster[channel][cell][desc].keys():
                        try:lsize.append(float(self.dcluster[channel][cell][desc][clust]["size"]))
                        except:pass
                    fsum.write(str(channel) + "\t" + str(cell) + "\t" + str(desc) + "\t" + str(mean(lsize)) + "\t" + str(std(lsize)) + "\n")

        fsum.close()


    def validationPredictor (self, typeCellChannel, pAC50All):

        dAC50All = toolbox.loadMatrix(pAC50All)

        dCASact = {}
        dpredict = {}
        dCASact[typeCellChannel] = []
        for CASID in dAC50All.keys():# have to change
            if dAC50All[CASID][typeCellChannel] != "NA":
                dCASact[typeCellChannel].append(CASID)
            if not CASID in dpredict.keys():
                if not path.exists(self.cDB.prSMIclean + CASID + ".smi"):
                    continue
                else:
                    smiles = toolbox.loadSMILES(self.cDB.prSMIclean + CASID + ".smi")
                    dpredict[CASID] = self.predictSMI(CASID, smiles, plot=1)


        prval = pathFolder.createFolder(self.prout + "validation/" + typeCellChannel + "/")

        for typeAssay in dCASact.keys():
            channel = "_".join(typeAssay.split("_")[1:])
            cell = typeAssay.split("_")[0]
            kpred = str(cell) + "_" + str(channel)
            ldesc = dpredict[dpredict.keys()[0]][kpred]
            filout = open(prval + typeCellChannel, "w")
            filout.write("CASID" + "\t".join(typeCellChannel) + "\n")
            for CASID in dpredict.keys():
                filout.write(CASID)
                for desc in ldesc:
                    filout.write("\t" + str(dpredict[CASID][kpred][desc]))
                filout.write("\n")
            filout.close()
        return 0


