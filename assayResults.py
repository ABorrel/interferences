import pathFolder
import runExternalSoft
from os import path

class assays:
    def __init__(self, pfilin, prout):
        self.pfilin = pfilin
        self.name = pfilin.split("/")[-1].split(".")[0]
        self.loadAssay()
        self.prout = prout


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

        pfilout = self.prout + "AC50_sample"
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








