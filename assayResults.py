import pathFolder
import runExternalSoft
import loadDB

from os import path
from numpy import mean

class assays:
    def __init__(self, pfilin, prout, prlog):
        self.pfilin = pfilin
        self.name = pfilin.split("/")[-1].split(".")[0]
        self.loadAssay()
        proutSP = prout + pfilin.split("/")[-1].split(".")[0] + "/"
        pathFolder.createFolder(proutSP)

        self.proutSP = proutSP
        self.prout = prout
        self.prlog = prlog


    def loadAssay(self):

        lout = []
        filin = open(self.pfilin, "r")
        llines = filin.readlines()
        filin.close()

        lheader =llines[0].strip().split("\t")

        for chemical in llines[1:]:
            lval = chemical.strip().split("\t")
            dtemp={}
            i = 0
            while i < len(lval):
                dtemp[lheader[i]] = lval[i]
                i += 1
            lout.append(dtemp)

        self.lchem = lout
        return lout

    def writeTable(self, lheaders, pfilout):

        filout = open(pfilout, "w")
        filout.write("\t".join(lheaders) + "\n")

        for chem in self.lchem:
            if chem["CAS"] == "":
                continue
            lw = []
            for header in lheaders:
                if chem[header] != "":
                    lw.append(str(chem[header]))
                else:
                    lw.append("NA")
            filout.write("\t".join(lw) + "\n")

        filout.close()


    def writeAC50(self):

        dAC50 = {}
        lsample = []

        pfilout = self.proutSP + "AC50_sample"
        if path.exists(pfilout):
            self.pAC50 = pfilout
            return pfilout

        for chem in self.lchem:
            if chem["AC50"] == "":
                chem["AC50"] = "NA"
            CAS = chem["CAS"]
            if not CAS in dAC50.keys():
                dAC50[CAS] = {}
            dAC50[CAS][chem["SAMPLE_DATA_TYPE"]] = chem["AC50"]
            if not chem["SAMPLE_DATA_TYPE"] in lsample:
                lsample.append(chem["SAMPLE_DATA_TYPE"])

        filout = open(pfilout, "w")
        filout.write("CAS\t" + "\t".join(lsample) + "\n")
        for casID in dAC50.keys():
            if casID == "":
                continue
            lw = []
            for sample in lsample:
                if sample in dAC50[casID].keys():
                    lw.append(str(dAC50[casID][sample]))
                else:
                    lw.append("NA")
            filout.write(str(casID) + "\t" + "\t".join(lw) + "\n")
        filout.close()
        self.pAC50 = pfilout

        return pfilout


    def combineAC50(self):

        dAC50 = {}
        lsample = []

        pfilout = self.proutSP + "AC50_combine"
        if path.exists(pfilout):
            self.pAC50 = pfilout
            return pfilout

        for chem in self.lchem:
            if chem["AC50"] == "":
                chem["AC50"] = "NA"
            CAS = chem["CAS"]
            if not CAS in dAC50.keys():
                dAC50[CAS] = {}
            dAC50[CAS][chem["SAMPLE_DATA_TYPE"]] = chem["AC50"]
            if not chem["SAMPLE_DATA_TYPE"] in lsample:
                lsample.append(chem["SAMPLE_DATA_TYPE"])

        filout = open(pfilout, "w")
        filout.write("CAS\tIC50\n")
        for casID in dAC50.keys():
            if casID == "":
                continue

            lM = []
            for sample in lsample:
                if sample in dAC50[casID].keys():
                    if dAC50[casID][sample] != "NA":
                        lM.append(float(dAC50[casID][sample]))
            if lM == []:
                M = "NA"
            else:
                M = mean(lM)

            filout.write(str(casID) + "\t" + str(M) + "\n")
        filout.close()
        self.pAC50 = pfilout




    def corAssay(self, assay2):

        prcor = pathFolder.createFolder(self.prout + "cor/")

        #merge assays
        dAC50 = {}

        for chem in self.lchem:
            if chem["AC50"] != "":
                if not chem["CAS"] in dAC50.keys():
                    dAC50[chem["CAS"]] = {}
                    dAC50[chem["CAS"]][chem["SAMPLE_DATA_TYPE"]] = {}

                    dAC50[chem["CAS"]][chem["SAMPLE_DATA_TYPE"]][self.name] = chem["AC50"]
                    dAC50[chem["CAS"]][chem["SAMPLE_DATA_TYPE"]][assay2.name] = "NA"
                else:
                    if not chem["SAMPLE_DATA_TYPE"] in dAC50[chem["CAS"]].keys():
                        dAC50[chem["CAS"]][chem["SAMPLE_DATA_TYPE"]] = {}
                        dAC50[chem["CAS"]][chem["SAMPLE_DATA_TYPE"]][self.name] = chem["AC50"]
                        dAC50[chem["CAS"]][chem["SAMPLE_DATA_TYPE"]][assay2.name] = "NA"


        for chem2 in assay2.lchem:
            if chem2["AC50"] != "":
                if chem2["CAS"] in dAC50.keys():

                    if chem2["SAMPLE_DATA_TYPE"] in dAC50[chem2["CAS"]].keys():
                        dAC50[chem2["CAS"]][chem2["SAMPLE_DATA_TYPE"]][assay2.name] = chem2["AC50"]
                    else:
                        dAC50[chem2["CAS"]][chem2["SAMPLE_DATA_TYPE"]] = {}
                        dAC50[chem2["CAS"]][chem2["SAMPLE_DATA_TYPE"]][assay2.name] = chem2["AC50"]
                        dAC50[chem2["CAS"]][chem2["SAMPLE_DATA_TYPE"]][self.name] = "NA"

                else:
                    dAC50[chem2["CAS"]] = {}
                    dAC50[chem2["CAS"]][chem2["SAMPLE_DATA_TYPE"]] = {}
                    dAC50[chem2["CAS"]][chem2["SAMPLE_DATA_TYPE"]][assay2.name] = chem2["AC50"]
                    dAC50[chem2["CAS"]][chem2["SAMPLE_DATA_TYPE"]][self.name] = "NA"


        # open different file
        dfiles = {}
        pfiloutGlobal = prcor + str(self.name + "_" + assay2.name)
        dfiles["global"] = open(pfiloutGlobal, "w")
        dfiles["global"].write("CAS\tAC50_1\tAC50_2\n")

        for casID in dAC50.keys():
            for sample in dAC50[casID].keys():
                if not sample in dfiles.keys():
                    psample = prcor + str(self.name + "_" + assay2.name + "_" + str(sample))
                    dfiles[sample] = open(psample, "w")
                    dfiles[sample].write("CAS\tAC50_1\tAC50_2\n")
                dfiles[sample].write(str(casID) + "\t" + str(dAC50[casID][sample][self.name]) + "\t" + str(dAC50[casID][sample][assay2.name]) + "\n")
            dfiles["global"].write(str(casID) + "\t" + str(dAC50[casID][sample][self.name]) + "\t" + str(dAC50[casID][sample][assay2.name]) + "\n")

        for filin in dfiles.keys():
            pfilin = dfiles[filin].name
            dfiles[filin].close()
            runExternalSoft.corplotR(pfilin)



    def AC50Distribution(self):

        prIC50 = self.proutSP + "histIC50/"
        pathFolder.createFolder(prIC50)

        # load table
        pAC50 = self.writeAC50()

        # run hist plot
        runExternalSoft.plotAC50(pAC50, prIC50, self.name.split("-")[1])



    def responseCurves(self):

        prresponse = self.proutSP + "responseCurve/"
        pathFolder.createFolder(prresponse)

        dresponse = {}
        for chem in self.lchem:
            casID = chem["CAS"]
            if casID == '':
                continue
            typein = chem["SAMPLE_DATA_TYPE"]

            if not casID in dresponse.keys():
                dresponse[casID] = {}
            dresponse[casID][typein] = {}
            dresponse[casID][typein]["DATA"] = []
            dresponse[casID][typein]["CONC"] = []

            if chem["AC50"] != "":
                dresponse[casID][typein]["AC50"] = chem["AC50"]
            else:
                dresponse[casID][typein]["AC50"] = "NA"

            if chem["CURVE_CLASS2"] != "":
                dresponse[casID][typein]["CURVE_CLASS2"] = chem["CURVE_CLASS2"]
            else:
                dresponse[casID][typein]["CURVE_CLASS2"] = "NA"

            i = 0
            while i < 16:
                kdata = "DATA" + str(i)
                kconc = "CONC" + str(i)
                if not kdata in chem.keys() or chem[kdata] == "":
                    dresponse[casID][typein]["DATA"].append("NA")
                else:
                    dresponse[casID][typein]["DATA"].append(chem["DATA" + str(i)])

                if not kconc in chem.keys() or chem[kconc] == "":
                    dresponse[casID][typein]["CONC"].append("NA")
                else:
                    dresponse[casID][typein]["CONC"].append(chem["CONC" + str(i)])

                i += 1

        # compute response curves
        for CASID in dresponse.keys():
            pCASout = prresponse + str(CASID)
            print pCASout
            filout = open(pCASout, "w")
            filout.write("CONC\tDATA\tFluorophores\tCurveType\n")

            i = 0
            while i < 16:
                for sample in dresponse[CASID].keys():
                    filout.write(str(dresponse[CASID][sample]["CONC"][i]) + "\t" + str(
                        dresponse[CASID][sample]["DATA"][i]) + "\t" + str(sample) + "\t" + str(
                        dresponse[CASID][sample]["CURVE_CLASS2"]) + "\n")
                i += 1
            filout.close()

        pAC50 = self.writeAC50()
        runExternalSoft.plotResponsiveCurve(prresponse, pAC50, self.proutSP)



    #not used
    def cor3assays(self, cassay2, cassay3):

        self.corAssay(cassay2)
        self.corAssay(cassay3)
        cassay2.corAssay(cassay3)



    def extractChemical(self, pSDFTox21):

        prSMI = self.prout + "SMI/"
        pathFolder.createFolder(prSMI)

        prSDF = self.prout + "SDF/"
        pathFolder.createFolder(prSDF)

        # load DB
        db = loadDB.sdfDB(pSDFTox21, "CASRN", self.prout)
        db.parseAll()

        # extract chemical
        lnotfind = []
        for chem in self.lchem:
            # print chem.keys()
            cas = chem["CAS"]
            if cas == "":
                continue


            flag = 0
            for cpdDB in db.lc:
                if cpdDB["CASRN"] == cas:
                    # sdf
                    pfilSDF = prSDF + str(cas) + ".sdf"
                    if not path.exists(pfilSDF):
                        filsdf = open(pfilSDF, "w")
                        filsdf.write(cpdDB["sdf"])
                        filsdf.close()

                    #Smile
                    pfilSMI = prSMI + str(cas) + ".smi"
                    if cpdDB["SMILES"] != "":
                        if not path.exists(pfilSMI) and cpdDB["SMILES"] != "":
                            filSMI = open(pfilSMI, "w")
                            filSMI.write(cpdDB["SMILES"])
                            filSMI.close()
                        flag = 1
                        break

            if flag == 0 and not cas in lnotfind:
                lnotfind.append(cas)


        logfile = open(self.prlog + self.name + "-extract.log", "w")
        logfile.write("\n".join(lnotfind))
        logfile.close()

        self.prSMI = prSMI
        self.prSDF = prSDF
