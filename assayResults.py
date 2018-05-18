import pathFolder
import runExternalSoft
import loadDB

from os import path
from numpy import mean, std
from PIL import Image, ImageFont, ImageDraw
from math import log10
from shutil import copyfile

font = ImageFont.truetype("OpenSans-Regular.ttf", size=30)


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
            #chemical without CAS -> PTC124 in luc
            if dtemp["CAS"] != "":
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


    def writeAC50(self, filtercurvefit = 1):

        dAC50 = {}
        lsample = []

        if filtercurvefit == 1:
            pfilout = self.proutSP + "AC50_sample_curve"
        else:
            pfilout = self.proutSP + "AC50_sample"

        if path.exists(pfilout):
            self.pAC50 = pfilout
            self.loadAC50()
            return pfilout

        if filtercurvefit == 1:
            if not "dresponse" in self.__dict__:
                self.responseCurves(drawn=0)


        for chem in self.lchem:
            if chem["AC50"] == "":
                chem["AC50"] = "NA"
            CAS = chem["CAS"]

            if filtercurvefit == 1:
                #print CAS, "CAS"
                #print self.dresponse[CAS]
                #print chem["SAMPLE_DATA_TYPE"]
                #print self.dresponse[CAS][chem["SAMPLE_DATA_TYPE"]]
                if float(self.dresponse[CAS][chem["SAMPLE_DATA_TYPE"]]["CURVE_CLASS2"]) >= 4:
                    chem["AC50"] = "NA"
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
        self.dAC50 = dAC50
        self.pAC50 = pfilout

        return pfilout


    def loadAC50(self):
        self.dAC50 = {}
        filin = open(self.pAC50, "r")
        llinecas = filin.readlines()
        filin.close()

        ltypeAC50 = llinecas[0].strip().split("\t")[1:]

        for linecas in llinecas[1:]:
            lac50 = linecas.strip().split("\t")
            cas = lac50[0]
            self.dAC50[cas] = {}
            i = 1
            while i < len(lac50):
                self.dAC50[cas][ltypeAC50[i-1]] = lac50[i]
                i += 1


    def combineAC50(self):

        dAC50out = {}

        pfilout = self.proutSP + "AC50_combine"
        if path.exists(pfilout) and path.getsize(pfilout) > 50:
            self.pAC50 = pfilout
            self.loadAC50()
            return pfilout


        if not "dAC50" in self.__dict__:
                self.writeAC50()

        lsample = self.dAC50[self.dAC50.keys()[0]].keys()

        filout = open(pfilout, "w")
        filout.write("CAS\tIC50\n")
        for casID in self.dAC50.keys():
            if casID == "":
                continue
            lM = []
            for sample in lsample:
                if sample in self.dAC50[casID].keys():
                    if self.dAC50[casID][sample] != "NA":
                        lM.append(float(self.dAC50[casID][sample]))
            if lM == []:
                M = "NA"
            else:
                M = mean(lM)
            dAC50out[casID] = {}
            dAC50out[casID]["set1"] = M
            filout.write(str(casID) + "\t" + str(M) + "\n")
        filout.close()
        self.pAC50 = pfilout
        self.dAC50 = dAC50out



    def corAC50(self):

        if not "pAC50" in self.__dict__:
            self.writeAC50()

        pcor = self.proutSP + "corAC50/"
        pathFolder.createFolder(pcor)

        dtypecurve = {}
        #define different type of fluo
        ltypefluo = []
        for chem in self.lchem:
            typefluo = chem["SAMPLE_DATA_TYPE"]
            if not typefluo in ltypefluo:
                ltypefluo.append(typefluo)

        for chem in self.lchem:
            casID = chem["CAS"]
            if not casID in dtypecurve.keys():
                dtypecurve[casID] = {}
                for typefluo in ltypefluo:
                    dtypecurve[casID][typefluo] = "NA"
            dtypecurve[casID][typefluo] = chem["CURVE_CLASS2"]

        pcurve = pcor + "curve"
        fcurve = open(pcurve, "w")
        fcurve.write("CAS" + "\t" + "\t".join(ltypefluo) + "\n")

        for casID in dtypecurve.keys():
            fcurve.write(casID + "\t" + "\t".join([str(dtypecurve[casID][k]) for k in ltypefluo]) + "\n")
        fcurve.close()


        runExternalSoft.corAC50(self.pAC50, pcurve, pcor)


        return


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



    def responseCurves(self, drawn=0):

        prresponse = self.proutSP + "responseCurve/"
        pathFolder.createFolder(prresponse)

        self.prresponse = prresponse

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
            if path.exists(pCASout) and path.getsize(pCASout) > 10:
                continue
            #print pCASout
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

        # draw plot
        if drawn == 1:
            pAC50 = self.writeAC50(filtercurvefit=0)
            runExternalSoft.plotResponsiveCurve(prresponse, pAC50, self.proutSP)

        self.dresponse = dresponse



    def crossResponseCurves(self, cAssays, prout):

        self.responseCurves(drawn=0)
        cAssays.responseCurves(drawn=0)

        self.writeAC50(filtercurvefit=0)
        cAssays.writeAC50(filtercurvefit=0)


        runExternalSoft.crossResponseCurve(self.prresponse, cAssays.prresponse, self.pAC50, cAssays.pAC50, prout)

        return


    def barplotCurveClass(self, prout):

        if not "dresponse" in self.__dict__:
            self.responseCurves(drawn=0)

        dfile = {}
        for CASID in self.dresponse.keys():
            for sample in self.dresponse[CASID].keys():
                if not sample in dfile.keys():
                    dfile[sample] = open(prout + str(sample), "w")
                    dfile[sample].write("CASID\tCurve\n")

                if self.dresponse[CASID][sample]["AC50"] != "NA" and float(self.dresponse[CASID][sample]["CURVE_CLASS2"]) < 4 :
                    dfile[sample].write(str(CASID) + "\t" + str(self.dresponse[CASID][sample]["CURVE_CLASS2"]) + "\n")

        for sample in dfile.keys():
            pfile = dfile[sample].name
            dfile[sample].close()
            runExternalSoft.barplotCurve(pfile)


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



    def rankingTop(self, nrank, prpng, prresult, delcurveinact=0):

        if not "dresponse" in self.__dict__:
            self.responseCurves(drawn=0)
        if not "dAC50" in self.__dict__:
            self.writeAC50()

        dval = {}
        for CASID in self.dAC50.keys():
            for sample in self.dAC50[CASID].keys():
                if not sample in dval.keys():
                    dval[sample] = []
                print sample, CASID
                dval[sample].append(self.dAC50[CASID][sample])
        for sample in dval.keys():
            dval[sample] = sorted(dval[sample])


        for sample in dval.keys():
            prsample = prresult + str(sample) + "/"
            pathFolder.createFolder(prsample)

            lcasflag = []
            rank = 1
            for val in dval[sample]:
                if val =="NA":
                    break
                for cas in self.dAC50.keys():
                    if cas in lcasflag:
                        continue
                    elif self.dAC50[cas][sample] == val:

                        if not sample in self.dresponse[cas].keys():
                            curve = self.dresponse[cas]["set1"]["CURVE_CLASS2"]
                        else:
                            curve = self.dresponse[cas][sample]["CURVE_CLASS2"]

                        if delcurveinact == 1 and float(curve) >= 4:
                            continue
                        pimageout = prsample + str(rank) + "_" + str(cas) + ".png"
                        pcaspng = prpng + cas + ".png"
                        if not path.exists(pcaspng) or path.getsize(pcaspng) < 10:
                            print "not found", pcaspng
                            continue
                        else:
                            writeLine = "CAS: " + str(cas) + "\nRANK: " + str(rank) + "\nAC50: " + str(val) + "\nCurve: " + str(curve)

                            img = Image.open(pcaspng)
                            imgnew = Image.new("RGBA", (580, 775), (250, 250, 250))
                            imgnew.paste(img, (0,0))
                            draw = ImageDraw.Draw(imgnew)
                            draw.text((10, 600), str(writeLine.split("\n")[0]), (0, 0, 0), font=font)
                            draw.text((10, 625), str(writeLine.split("\n")[1]), (0, 0, 0), font=font)
                            draw.text((10, 650), str(writeLine.split("\n")[2]), (0, 0, 0), font=font)
                            draw.text((10, 675), str(writeLine.split("\n")[3]), (0, 0, 0), font=font)
                            imgnew.save(pimageout)

                        rank = rank + 1
                        lcasflag.append(cas)

                if rank > nrank:
                    break



    def summarize(self, prout):

        # table nbChem/active/inactive/M and SD AC50 + M and SD log10
        if not "dresponse" in self.__dict__:
            self.responseCurves(drawn=0)
        if not "dAC50" in self.__dict__:
            self.writeAC50()



        dsum = {}
        for sample in self.dresponse[self.dresponse.keys()[0]].keys():
            dsum[sample] = {}
            dsum[sample]["act"] = []
            dsum[sample]["inact"] = []

        for CASID in self.dresponse.keys():
            for sample in dsum.keys():
                if self.dresponse[CASID][sample]["AC50"] == "NA":
                    dsum[sample]["inact"].append(CASID)
                elif float(self.dresponse[CASID][sample]["CURVE_CLASS2"]) >=4 :
                    dsum[sample]["inact"].append(CASID)
                else:
                    dsum[sample]["act"].append(float(self.dresponse[CASID][sample]["AC50"]))


        pfilout = prout + "summarize.txt"
        filout = open(pfilout, "w")
        filout.write("Raw\tNb Chemical\tNb active\tNb inactive\tMean(AC50)\tSD(AC50)\tMean(-(logAC50))\tSD(-log(AC50))\n")
        for sample in dsum.keys():
            filout.write(str(sample) + "\t" + str(len(dsum[sample]["act"])+len(dsum[sample]["inact"])) + "\t" +
                         str(len(dsum[sample]["act"])) + "\t" +
                         str(len(dsum[sample]["inact"])) + "\t" + str(mean(dsum[sample]["act"])) + "\t" +
                         str(std(dsum[sample]["act"])) + "\t" +
                         str(mean([-log10(x) for x in dsum[sample]["act"]])) + "\t" +
                         str(std([-log10(x) for x in dsum[sample]["act"]])) + "\n")
        filout.close()



    def drawVennPlot(self, pranalysis, prPNG):

        if not "dresponse" in self.__dict__:
            self.responseCurves(drawn=0)

        lsample = self.dresponse[self.dresponse.keys()[0]].keys()

        i = 0
        imax = len(lsample)
        while i < imax-1:
            j = i + 1
            while j < imax:
                sample1 = lsample[i]
                sample2 = lsample[j]

                prsubpng = pathFolder.createFolder(pranalysis + str(sample1) + "-" + str(sample2) + "/")
                for CASID in self.dresponse.keys():
                    if self.dresponse[CASID][sample1]["AC50"] == "NA" or self.dresponse[CASID][sample2]["AC50"] == "NA":
                        continue
                    if float(self.dresponse[CASID][sample1]["CURVE_CLASS2"]) >= 4 or float(self.dresponse[CASID][sample2]["CURVE_CLASS2"]) >= 4:
                        continue

                    if path.exists(prPNG + CASID + ".png"):
                        copyfile(prPNG + CASID + ".png", prsubpng + CASID + ".png")
                j += 1
            i += 1
        runExternalSoft.vennPlot(self.pAC50, pranalysis)


    def createPCA(self, pdesc1D2D, pAC50, corval, maxQuantile, prPCA):

        # output
        pdesc1D2Dclean = prPCA + "descClean.csv"

        if not path.exists(pdesc1D2Dclean):

            if path.exists(pdesc1D2D) and path.getsize(pdesc1D2D) > 10:
                # preproc
                runExternalSoft.dataManager(pdesc1D2D, 0, corval, maxQuantile, prPCA)
            else:
                print "Error ->", pdesc1D2D

        runExternalSoft.drawPCA(pdesc1D2Dclean, pAC50, prPCA)
