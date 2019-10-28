from multiprocessing.sharedctypes import class_cache

from nams import nams
import runExternalSoft
from os import path, listdir
import toolbox
from copy import deepcopy
import pathFolder


class MCSMatrix:

    def __init__(self, prSMI , prout):

        self.prSMI = prSMI
        self.prout = prout

    def computeMatrixMCS(self):

        pfiloutTanimoto = self.prout + "tanimoto.txt"
        pfiloutNBatomMax = self.prout + "maxAtom.txt"

        if path.exists(pfiloutTanimoto) and path.exists(pfiloutNBatomMax):
            print "in"
            #dMCSTanimoto = toolbox.loadMatrix(pfiloutTanimoto)
            #dMCSMax = toolbox.loadMatrix(pfiloutNBatomMax)
            #self.MCSTanimoto = dMCSTanimoto
            #self.MCSMax= dMCSMax

        else:
            lSMI = listdir(self.prSMI)
            lcas = [smi[:-4] for smi in lSMI]
            jmax = len(lcas)

            i = 0
            imax = jmax


            while i < imax:
                print i, imax
                print self.prout + str(lcas[i])
                if path.exists(self.prout + str(lcas[i])):
                    i += 1
                    continue

                dTemp = {}
                j = i
                while j < jmax:
                    print i, j
                    if not lcas[j] in dTemp.keys():
                        fsmii = open(self.prSMI + lcas[i] + ".smi", "r")
                        smii = fsmii.readlines()[0].strip()
                        fsmii.close()
                        fsmij = open(self.prSMI + lcas[j] + ".smi", "r")
                        smij = fsmij.readlines()[0].strip()
                        fsmij.close()

                        ltanimoto_max = get_Tanimoto(smii, smij)
                        dTemp[lcas[j]] = ltanimoto_max
                    j += 1
                i += 1

                # write temp matrix
                filout = open(self.prout + str(lcas[i]), "w")
                filout.write("CAS\tMCS\tMax atoms\n")
                for k in dTemp.keys():
                    filout.write(str(k) + "\t" + str(dTemp[k][0]) + "\t" + str(dTemp[k][1]) + "\n")
                filout.close()


            # merge matrix
            dMaxMCS = {}
            dTanimoto = {}
            for nfile in listdir(self.prout):
                if nfile[-3:] != "txt":
                    cmpdID1 = nfile
                    dMaxMCS[cmpdID1] = {}
                    dTanimoto[cmpdID1] = {}

                    filin = open(self.prout + cmpdID1, "r")
                    lscores = filin.readlines()
                    for score in lscores[1:]:
                        lelemscore = score.strip().split("\t")
                        cmpdID2 = lelemscore[0]
                        dTanimoto[cmpdID1][cmpdID2] = lelemscore[1]
                        dMaxMCS[cmpdID1][cmpdID2] = lelemscore[2]


            filoutTanimoto = open(pfiloutTanimoto, "w")
            filoutNBatomMax = open(pfiloutNBatomMax, "w")


            filoutTanimoto.write("\t".join(lcas) + "\n")
            filoutNBatomMax.write("\t".join(lcas) + "\n")

            for cmpdID1 in lcas:
                lwTanimoto = []
                lwMax = []
                for cmpdID2 in lcas:
                    try: lwTanimoto.append(str(dTanimoto[cmpdID1][cmpdID2]))
                    except: lwTanimoto.append(str(dTanimoto[cmpdID2][cmpdID1]))
                    try: lwMax.append(str(dMaxMCS[cmpdID1][cmpdID2][1]))
                    except: lwMax.append(str(dMaxMCS[cmpdID2][cmpdID1][1]))
                filoutTanimoto.write(cmpdID1 + "\t" + "\t".join(lwTanimoto) + "\n")
                filoutNBatomMax.write(cmpdID1 + "\t" + "\t".join(lwMax) + "\n")
            filoutTanimoto.close()
            filoutNBatomMax.close()


    def selectAnalogsMatrix(self, compoundID, tanimotocutoff = 0.8):

        if not "dMCSTanimoto" in dir(self):
            self.computeMatrixMCS()

        dTanimoto = self.MCSTanimoto[compoundID]
        daff = self.Aff
        dMax = {}

        print len(dTanimoto.keys())
        for cpID in dTanimoto.keys():
            if float(dTanimoto[cpID]) < tanimotocutoff:
                del dTanimoto[cpID]

        lcID = dTanimoto.keys()
        dTanimoto = {}


        for cID in lcID:
            dTanimoto[cID] = {}
            dMax[cID] = {}
            for cID2 in lcID:
                dTanimoto[cID][cID2] = self.MCSTanimoto[cID][cID2]
                dMax[cID][cID2] = self.MCSMax[cID][cID2]

        # write file
        pfiloutTanimoto = self.prout + compoundID + str(tanimotocutoff) + "_Tanimoto"
        pfiloutDMAX = self.prout + compoundID + str(tanimotocutoff) + "_MAX"
        filoutTanimoto = open(pfiloutTanimoto, "w")
        filoutDMAX = open(pfiloutDMAX, "w")

        filoutDMAX.write("\t".join(lcID) + "\n")
        filoutTanimoto.write("\t".join(lcID) + "\n")

        for cID in lcID:
            filoutTanimoto.write(str(cID) + "\t" + "\t".join([str(dTanimoto[cID][i]) for i in lcID]) + "\n")
            filoutDMAX.write(str(cID) + "\t" + "\t".join([str(dMax[cID][i]) for i in lcID]) + "\n")
        filoutDMAX.close()
        filoutTanimoto.close()

        pfiloutAff = self.prout + compoundID + "_Aff"
        filoutAff = open(pfiloutAff, "w")
        filoutAff.write("ID\tIC50\n")
        for cID in lcID:
            filoutAff.write(str(cID) + "\t" + daff[cID] + "\n")
        filoutAff.close()

        runExternalSoft.MatrixMCS(pfiloutTanimoto, pfiloutAff, pfiloutDMAX, 1)

        return


    def selectCluster(self, pfilecluster, prout=""):

        if prout == "":
            prout = self.prout + "-".join(pfilecluster[0:-4].split("_")[1:]).replace(".", "") + "/"
        else:
            prout = prout + "-".join(pfilecluster[0:-4].split("_")[1:]).replace(".","") + "/MCS/"
        pathFolder.createFolder(prout)

        if not "dMCSTanimoto" in dir(self):
            self.computeMatrixMCS()

        filecluster = open(pfilecluster, "r")
        llinesCluster = filecluster.readlines()
        filecluster.close()
        dcluster = {}

        for lineCluster in llinesCluster[1:]:
            cluster = lineCluster.strip().split(",")[-1].replace("\"", "")
            print cluster
            if not cluster in dcluster.keys():
                dcluster[cluster] = []
            compoundID = lineCluster.strip().split("\"")[1]
            dcluster[cluster].append(compoundID)

        dTanimoto = self.MCSTanimoto
        daff = self.Aff

        # reduce matrix
        for cluster in dcluster.keys():
            lcpd = dcluster[cluster]
            # tanimoto temp
            dtanimotocluster = {}
            daffcluster = deepcopy(daff)
            dMaxcluster = {}

            pfiloutclustertanimoto = prout + str(cluster) + "_Tanimoto"
            pfiloutclusterDmax = prout + str(cluster) + "_MAX"
            pfiloutclusterAff= prout + str(cluster) + "_aff"


            for cpID in lcpd:
                dtanimotocluster[cpID] = {}
                dMaxcluster[cpID] = {}
                for cpID2 in lcpd:
                    try:
                        dtanimotocluster[cpID][cpID2] = dTanimoto[cpID][cpID2]
                        dMaxcluster[cpID][cpID2] = self.MCSMax[cpID][cpID2]
                    except:
                        del dMaxcluster[cpID]
                        del dtanimotocluster[cpID]
                        break

            # write file
            filoutTanimoto = open(pfiloutclustertanimoto, "w")
            filoutDMAX = open(pfiloutclusterDmax, "w")
            filoutAff = open(pfiloutclusterAff, "w")

            filoutDMAX.write("\t".join(lcpd) + "\n")
            filoutTanimoto.write("\t".join(lcpd) + "\n")


            for cID in lcpd:
                if not cID in dtanimotocluster.keys():
                    continue
                filoutTanimoto.write(str(cID) + "\t" + "\t".join([str(dtanimotocluster[cID][i]) for i in lcpd]) + "\n")
                filoutDMAX.write(str(cID) + "\t" + "\t".join([str(dMaxcluster[cID][i]) for i in lcpd]) + "\n")
            filoutDMAX.close()
            filoutTanimoto.close()


            filoutAff.write("ID\tIC50\n")
            for cID in lcpd:
                try: filoutAff.write(str(cID) + "\t" + daff[cID] + "\n")
                except:continue
            filoutAff.close()

            runExternalSoft.MatrixMCS(pfiloutclustertanimoto, pfiloutclusterAff, pfiloutclusterDmax, 3)# option 3 to compute the difference of affinity

        self.clusters = dcluster



    def selectAnalogs(self, compoundID, cutoffMCS = 0.8):

        lout = []

        pfilout = self.prout + "analogs_" + str(compoundID) + "_" + str(cutoffMCS) + ".txt"
        pfiloutglobal = self.prout + "MCSglobal_" + str(compoundID) + "_" + str(cutoffMCS) + ".txt"
        if path.exists(pfiloutglobal):
            filin = open(pfiloutglobal, "r")
            lcpd = filin.readlines()
            filin.close()

            filout = open(pfilout, "w")
            filout.write("CMPD_CHEMBLID\tCANONICAL_SMILES\tTanimoto\tMaxAtom\n")
            for linecpd in lcpd:
                lelem = linecpd.strip().split("\t")
                dcpd = {}
                dcpd["CMPD_CHEMBLID"] = lelem[0]
                dcpd["CANONICAL_SMILES"] = lelem[1]
                dcpd["Tanimoto"] = lelem[2]
                dcpd["MaxAtom"] = lelem[3]


                if dcpd["Tanimoto"] >= cutoffMCS:
                    filout.write(str(dcpd["CMPD_CHEMBLID"]) + "\t" + str(dcpd["CANONICAL_SMILES"]) + "\t" +
                                str(dcpd["Tanimoto"]) + "\t" + str(dcpd["MaxAtom"]) + "\n")
                    lout.append(dcpd)

            filout.close()
            self.lanalogs = lout
            return lout


        filout = open(pfilout, "w")
        filoutglobal = open(pfiloutglobal, "w")
        filout.write("CMPD_CHEMBLID\tCANONICAL_SMILES\tTanimoto\tMaxAtom\n")
        filoutglobal.write("CMPD_CHEMBLID\tCANONICAL_SMILES\tTanimoto\tMaxAtom\n")
        for cp in self.sdata:
            if cp["CMPD_CHEMBLID"] == compoundID:
                smileref = cp["CANONICAL_SMILES"]
                break

        i = 0
        nbcp = len(self.sdata)
        while i < nbcp:
            try: tanimoto = get_Tanimoto(smileref, self.sdata[i]["CANONICAL_SMILES"])
            except:
                i += 1
                continue

            filoutglobal.write(str(self.sdata[i]["CMPD_CHEMBLID"]) + "\t" + str(self.sdata[i]["CANONICAL_SMILES"]) + "\t"
                         + str(tanimoto[0]) + "\t" + str(tanimoto[1]) + "\n")

            if tanimoto[0] >= cutoffMCS:
                filout.write(str(self.sdata[i]["CMPD_CHEMBLID"]) + "\t" + str(self.sdata[i]["CANONICAL_SMILES"]) + "\t"
                             + str(tanimoto[0]) + "\t" + str(tanimoto[1]) + "\n")
                lout = deepcopy(self.sdata[i])
            print i, nbcp
            i += 1

        filout.close()
        filoutglobal.close()

        self.lanalogs = lout
        return lout

def get_Tanimoto(smile1, smile2):
    ms = nams.Nams()

    #print(smile1, smile2)
    mol_t1 = ("smi", smile1)
    mol_t2 = ("smi", smile2)
    try:mol1, mol_info1 = ms.get_mol_info(mol_t1[0], mol_t1[1])
    except:return ["NA", "NA"]
    try:mol2, mol_info2 = ms.get_mol_info(mol_t2[0], mol_t2[1])
    except:return ["NA", "NA"]

    #print mol_info1.keys()
    #print mol_info2.keys()
    # similarity combination
    sim12, d_atoms12 = ms.get_similarity(mol_info1, mol_info2)
    sim21, d_atoms21 = ms.get_similarity(mol_info2, mol_info1)
    sim11, d_atoms11 = ms.get_similarity(mol_info1, mol_info1)
    sim22, d_atoms22 = ms.get_similarity(mol_info2, mol_info2)

    #print sim11, d_atoms11

    #test similarity
    sizeMCS = len(d_atoms12.keys())
    nbatomdiff = max([len(d_atoms11.keys()), len(d_atoms22.keys())]) - len(d_atoms12.keys())
    text = str(sizeMCS) + "-" + str(nbatomdiff)
    #print nbatomdiff

    ##### TEST FOR NAMS SCORE #######
    #################################

    #print len(d_atoms12.keys()), len(d_atoms11.keys()), len(d_atoms22.keys())

    #ks=d_atoms12.keys()
    #ks.sort()

    #Print the atomic alignment between the molecules and its similarity score
    #for k in ks:
    #    print "\t%5d (%3s)  -%3d  (%3s) --> %6.2f" % (k[0], mol1.atoms[k[0]-1].OBAtom.GetType(),
    #                                                  k[1], mol2.atoms[k[1]-1].OBAtom.GetType(),
    #                                                  d_atoms12[k])


    # sim11
    #ks=d_atoms11.keys()
    #ks.sort()

    #print("###############################")
    #Print the atomic alignment between the molecules and its similarity score
    #for k in ks:
    #    print "\t%5d (%3s)  -%3d  (%3s) --> %6.2f" % (k[0], mol1.atoms[k[0]-1].OBAtom.GetType(),
    #                                                  k[1], mol2.atoms[k[1]-1].OBAtom.GetType(),
    #                                                  d_atoms11[k])



    # based on a Jaccard score
    #print sim12, sim21, sim11, sim22
    #print d_atoms12
    score12 = sim12 / (sim11 + sim22 - sim12)
    score21 = sim21 / (sim11 + sim22 - sim21)
    #print score12

    return [score12, text]


#a = "/home/borrel/Yue_project/result/Pi_LSR_LGDsimilarity/KS-5_1B38/cycle-CON_1N3_3ULI.smi"
#b = "/home/borrel/Yue_project/result/Pi_LSR_LGDsimilarity/KS-5_1B38/onlyC_PM1_1PYE.smi"

#get_Tanimoto(a,b)