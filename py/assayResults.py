import pathFolder
import runExternalSoft
import loadDB
import toolbox
import cytox

from os import path
from numpy import mean, std
from PIL import Image, ImageFont, ImageDraw
from math import log10
from shutil import copyfile
from re import search

font = ImageFont.truetype("OpenSans-Regular.ttf", size=30)


class assays:
    def __init__(self, pfilin, curvecutoff, effcutoff, curvePositive, curveNegative, prcytox, prout, prlog):
        self.pfilin = pfilin
        self.prcytox = prcytox
        self.name = pfilin.split("/")[-1].split(".")[0]
        self.loadAssay()
        proutSP = prout + pfilin.split("/")[-1].split(".")[0] + "/"
        pathFolder.createFolder(proutSP)

        self.proutSP = proutSP
        self.prout = prout
        self.prlog = prlog

        self.curveCutoff = curvecutoff
        self.curvePositive = curvePositive
        self.curveNegative = curveNegative

        self.effcutoff = effcutoff

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


    def writeAC50(self, filtercurvefit = 1, filterefficacy = 1, filterburst = 1, combine = 0):

        dAC50 = {}
        lsample = []

        pfilout = self.proutSP + "AC50_sample"
        peff = pathFolder.createFolder(self.proutSP + "efficient/") + "eff"
        if filtercurvefit == 1:
            pfilout = pfilout + "_curve"
            peff = peff + "_curve"
        if filterefficacy == 1:
            pfilout = pfilout + "_eff"
            peff = peff + "_eff"
        if filterburst == 1:
            pfilout = pfilout + "_burst"
            peff = peff + "_burst"
        if combine == 1:
            pfilout = pfilout + "_combine"
            peff = peff + "_combine"


        if path.exists(pfilout) and path.exists(peff):
            self.pAC50 = pfilout
            self.loadAC50()
            return pfilout

        #if filtercurvefit == 1:
        if not "dresponse" in self.__dict__:
            self.responseCurves(drawn=0)

        if filterburst == 1:
            dcytox = cytox.parsepdf(self.prcytox, self.prout)

        for chem in self.lchem:
            if chem["AC50"] == "":
                chem["AC50"] = "NA"
            CAS = chem["CAS"]
            if filtercurvefit == 1:
                if abs(float(self.dresponse[CAS][chem["SAMPLE_DATA_TYPE"]]["CURVE_CLASS2"])) >= self.curveCutoff:
                    chem["AC50"] = "NA"
                if self.curveNegative == 0 and float(self.dresponse[CAS][chem["SAMPLE_DATA_TYPE"]]["CURVE_CLASS2"]) < 0.0:
                    chem["AC50"] = "NA"
                if self.curvePositive == 0 and float(self.dresponse[CAS][chem["SAMPLE_DATA_TYPE"]]["CURVE_CLASS2"]) > 0.0:
                    chem["AC50"] = "NA"

            if filterefficacy == 1:
                # for luciferase

                if abs(float(self.dresponse[CAS][chem["SAMPLE_DATA_TYPE"]]["EFFICACY"])) < self.effcutoff:
                    chem["AC50"] = "NA"

            if not CAS in dAC50.keys():
                dAC50[CAS] = {}
            dAC50[CAS][chem["SAMPLE_DATA_TYPE"]] = chem["AC50"]
            if not chem["SAMPLE_DATA_TYPE"] in lsample:
                lsample.append(chem["SAMPLE_DATA_TYPE"])

        if combine == 1:
            dAC50 = self.combineIC50(dAC50)


        if filterburst == 1:
            self.filterCytox(dcytox, dAC50)


        lsamples = dAC50[dAC50.keys()[0]].keys()
        filout = open(pfilout, "w")
        feff = open(peff, "w")
        filout.write("CAS\t" + "\t".join(lsamples) + "\n")
        feff.write("CAS\t" + "\t".join([str(s + "\t" + s + "_eff") for s in lsamples]) + "\n")
        for casID in dAC50.keys():
            if casID == "":
                continue
            lw = []
            lweff = []
            for sample in lsamples:
                if sample in dAC50[casID].keys():
                    lw.append(str(dAC50[casID][sample]))
                    lweff.append(str(dAC50[casID][sample]))
                else:
                    lw.append("NA")
                    lweff.append("NA")
                if sample == "IC50":
                    leff = []
                    for setIC50 in self.dresponse[casID].keys():
                        leff.append(float(self.dresponse[casID][setIC50]["EFFICACY"]))
                    lweff.append(str(abs(mean(leff))))
                else:
                    lweff.append(str(abs(float(self.dresponse[casID][sample]["EFFICACY"]))))
            filout.write(str(casID) + "\t" + "\t".join(lw) + "\n")
            feff.write(str(casID) + "\t" + "\t".join(lweff) + "\n")
        filout.close()
        feff.close()

        runExternalSoft.plotAC50VSEff(peff)

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


    def combineIC50(self, dAC50):

        dAC50out = {}

        lsample = dAC50[dAC50.keys()[0]].keys()

        for casID in dAC50.keys():
            if casID == "":
                continue
            lM = []
            for sample in lsample:
                if sample in dAC50[casID].keys():
                    if dAC50[casID][sample] != "NA":
                        lM.append(float(dAC50[casID][sample]))
            if len(lM) < 3:
                M = "NA"
            else:
                M = mean(lM)
            dAC50out[casID] = {}
            dAC50out[casID]["IC50"] = M

        return dAC50out



    def corAC50(self):

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

        prAC50 = self.proutSP + "histAC50/"
        pathFolder.createFolder(prAC50)

        # run hist plot
        runExternalSoft.plotAC50(self.pAC50, prAC50, self.name.split("-")[1])



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

            if chem["EFFICACY"] != "":
                dresponse[casID][typein]["EFFICACY"] = chem["EFFICACY"]
            else:
                dresponse[casID][typein]["EFFICACY"] = 0.0

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

        # compute response curve
        self.dresponse = dresponse
        # draw plot
        if drawn == 1:
            runExternalSoft.plotResponsiveCurve(prresponse, self.pAC50, self.proutSP)




    def crossResponseCurves(self, cAssays, prout):

        self.responseCurves(drawn=0)
        cAssays.responseCurves(drawn=0)

        self.writeAC50(filtercurvefit=0, filterburst=0)
        cAssays.writeAC50(filtercurvefit=0, filterburst=0)

        runExternalSoft.crossResponseCurve(self.prresponse, cAssays.prresponse, self.pAC50, cAssays.pAC50, prout)

        # reimplement with filters
        self.writeAC50()
        cAssays.writeAC50()


    def barplotCurveClass(self, prout):

        if not "dresponse" in self.__dict__:
            self.responseCurves(drawn=0)

        dfile = {}
        dfile["all"] = open(prout + "all", "w")
        dfile["all"].write("CASID\tCurves\tAff\n")
        for CASID in self.dresponse.keys():
            for sample in self.dresponse[CASID].keys():
                if not sample in dfile.keys():
                    dfile[sample] = open(prout + str(sample), "w")
                    dfile[sample].write("CASID\tCurves\tAff\n")

                if self.dresponse[CASID][sample]["AC50"] == "NA":
                    continue
                if abs(float(self.dresponse[CASID][sample]["CURVE_CLASS2"])) >= self.curveCutoff:
                    continue
                if self.curveNegative == 0 and float(self.dresponse[CASID][sample]["CURVE_CLASS2"]) < 0.0:
                    continue
                if self.curvePositive == 0 and float(self.dresponse[CASID][sample]["CURVE_CLASS2"]) > 0.0:
                    continue

                dfile[sample].write(str(CASID) + "\t" + str(self.dresponse[CASID][sample]["CURVE_CLASS2"]) + "\t" + str(self.dresponse[CASID][sample]["AC50"]) + "\n")

                # case all color
                if search("_n", sample):
                    dfile["all"].write(str(CASID) + "\t" + str(self.dresponse[CASID][sample]["CURVE_CLASS2"]) + "\t" + str(self.dresponse[CASID][sample]["AC50"]) + "\n")

                #case all set
                dfile["all"].write(str(CASID) + "\t" + str(self.dresponse[CASID][sample]["CURVE_CLASS2"]) + "\t" + str(self.dresponse[CASID][sample]["AC50"]) + "\n")


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

                        if delcurveinact == 1:
                            if abs(float(curve)) >= self.curveCutoff:
                                continue
                            if self.curvePositive == 0 and float(curve) > 0:
                                continue
                            if self.curveNegative == 0 and float(curve) < 0:
                                continue

                        pcaspng = prpng + str(cas) + ".png"
                        if not path.exists(pcaspng):
                            continue
                        else:
                            pimageout = prresult + str(cas) + ".png"
                        writeLine = "CAS: " + str(cas) + "\nRANK: " + str(rank) + "\nAC50: " + str(
                            val) + "\nCurve: " + str(curve)

                        img = Image.open(pcaspng)
                        imgnew = Image.new("RGBA", (580, 775), (250, 250, 250))
                        imgnew.paste(img, (0, 0))
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

        dsum = {}
        lsample = self.dAC50[self.dAC50.keys()[0]].keys()
        if "CASID" in lsample:
            del lsample[lsample.index("CASID")]
        for sample in lsample:
            if sample == "CASID":
                continue
            dsum[sample] = {}
            dsum[sample]["act"] = []
            dsum[sample]["inact"] = []

        for CASID in self.dAC50.keys():
            for sample in lsample:
                if self.dAC50[CASID][sample] == "NA":
                    dsum[sample]["inact"].append(CASID)
                else:
                    dsum[sample]["act"].append(float(self.dAC50[CASID][sample]))


        pfilout = prout + "summarize_" + str(self.pAC50.split("/")[-1])

        filout = open(pfilout, "w")
        filout.write("Raw\tNb Chemical\tNb active\tNb inactive\tMean(AC50)\tSD(AC50)\tMean(-(logAC50))\tSD(-log(AC50))\n")
        for sample in dsum.keys():
            filout.write("%s\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\n"% (sample,
                         len(dsum[sample]["act"])+len(dsum[sample]["inact"]),
                             len(dsum[sample]["act"]), len(dsum[sample]["inact"]), mean(dsum[sample]["act"]),
                             std(dsum[sample]["act"]), mean([-log10(x) for x in dsum[sample]["act"]]),
                             std([-log10(x) for x in dsum[sample]["act"]])))
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
                    if self.dresponse[CASID][sample1]["AC50"] == "NA":
                        continue
                    if abs(float(self.dresponse[CASID][sample1]["CURVE_CLASS2"])) >= self.curveCutoff:
                        continue
                    if self.curveNegative == 0 and float(self.dresponse[CASID][sample2]["CURVE_CLASS2"]) < 0:
                        continue
                    if self.curvePositive == 0 and float(self.dresponse[CASID][sample2]["CURVE_CLASS2"]) > 0:
                        continue

                    if path.exists(prPNG + CASID + ".png"):
                        copyfile(prPNG + CASID + ".png", prsubpng + CASID + ".png")
                j += 1
            i += 1
        runExternalSoft.vennPlot(self.pAC50, pranalysis)


    def createPCA(self, pdesc1D2D, pAC50, corval, maxQuantile, nbNA, prPCA):

        # output
        pdesc1D2Dclean = prPCA + "descClean.csv"

        if not path.exists(pdesc1D2Dclean):

            if path.exists(pdesc1D2D) and path.getsize(pdesc1D2D) > 10:
                # preproc
                runExternalSoft.dataManager(pdesc1D2D, 0, corval, maxQuantile, nbNA, prPCA)
            else:
                print "Error ->", pdesc1D2D

        runExternalSoft.drawPCA(pdesc1D2Dclean, pAC50, prPCA)



    def createMDS(self, pdesc1D2D, pAC50, corval, maxQuantile, prMDS):

        # output
        pdesc1D2Dclean = prMDS + "descClean.csv"

        if not path.exists(pdesc1D2Dclean):

            if path.exists(pdesc1D2D) and path.getsize(pdesc1D2D) > 10:
                # preproc
                runExternalSoft.dataManager(pdesc1D2D, 0, corval, maxQuantile, prMDS)
            else:
                print "Error ->", pdesc1D2D

        runExternalSoft.drawMDS(pdesc1D2Dclean, pAC50, prMDS)


    def filterCytox(self, dcytox, dAC50):


        linter = list(set(dcytox.keys()) & set(dAC50.keys()))
        print len(linter), "-> l662"

        for CASinter in linter:
            for channel in dAC50[CASinter].keys():
                if dAC50[CASinter][channel] != "NA":
                    if float(dAC50[CASinter][channel]) > float(dcytox[CASinter]["CytoxMin"]): # or cytoxMin
                        dAC50[CASinter][channel] = "NA"





def histogramAC50(pAC50All, prhist):

    runExternalSoft.histGlobal(pAC50All, prhist)




def mergeAssays(cluc, chepg2, chek293):

    pfilout = cluc.prout + "AC50_all"
    if path.exists(pfilout):
        return pfilout
    filout = open(pfilout, "w")
    lheader = ["CASID", "Luc_IC50", "hepg2_med_blue", "hepg2_med_green", "hepg2_med_red", "hepg2_med_blue_n",
               "hepg2_med_green_n", "hepg2_med_red_n", "hek293_med_blue", "hek293_med_green", "hek293_med_red",
               "hek293_med_blue_n", "hek293_med_green_n", "hek293_med_red_n", "hepg2_cell_blue", "hepg2_cell_green",
               "hepg2_cell_red", "hepg2_cell_blue_n", "hepg2_cell_green_n", "hepg2_cell_red_n", "hek293_cell_blue",
               "hek293_cell_green", "hek293_cell_red", "hek293_cell_blue_n", "hek293_cell_green_n", "hek293_cell_red_n"]

    lchannel = ["med_blue", "med_green", "med_red", "med_blue_n", "med_green_n", "med_red_n", "cell_blue", "cell_green",
                "cell_red", "cell_blue_n", "cell_green_n", "cell_red_n"]
    filout.write("\t".join(lheader) + "\n")

    lCAS = list(set(cluc.dAC50.keys() + chepg2.dAC50.keys()))
    lCAS = list(set(lCAS + chek293.dAC50.keys()))
    for casID in lCAS:
        #print chepg2.dAC50[casID]
        #print chek293.dAC50[casID]
        if not casID in chepg2.dAC50.keys():
            chepg2.dAC50[casID] = {}
            for i in lchannel:
                chepg2.dAC50[casID][i] = "NA"
        if not casID in cluc.dAC50.keys():
            cluc.dAC50[casID] = {}
            cluc.dAC50[casID]["IC50"] = "NA"
        if not casID in chek293.dAC50.keys():
            chek293.dAC50[casID] = {}
            for i in lchannel:
                chek293.dAC50[casID][i] = "NA"


        linew = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" \
                "\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (casID, cluc.dAC50[casID]["IC50"],
                                                                            chepg2.dAC50[casID]["med_blue"],
                                                                            chepg2.dAC50[casID]["med_green"],
                                                                            chepg2.dAC50[casID]["med_red"],
                                                                            chepg2.dAC50[casID]["med_blue_n"],
                                                                            chepg2.dAC50[casID]["med_green_n"],
                                                                            chepg2.dAC50[casID]["med_red_n"],
                                                                            chek293.dAC50[casID]["med_blue"],
                                                                            chek293.dAC50[casID]["med_green"],
                                                                            chek293.dAC50[casID]["med_red"],
                                                                            chek293.dAC50[casID]["med_blue_n"],
                                                                            chek293.dAC50[casID]["med_green_n"],
                                                                            chek293.dAC50[casID]["med_red_n"],
                                                                            chepg2.dAC50[casID]["cell_blue"],
                                                                            chepg2.dAC50[casID]["cell_green"],
                                                                            chepg2.dAC50[casID]["cell_red"],
                                                                            chepg2.dAC50[casID]["cell_blue_n"],
                                                                            chepg2.dAC50[casID]["cell_green_n"],
                                                                            chepg2.dAC50[casID]["cell_red_n"],
                                                                            chek293.dAC50[casID]["cell_blue"],
                                                                            chek293.dAC50[casID]["cell_green"],
                                                                            chek293.dAC50[casID]["cell_red"],
                                                                            chek293.dAC50[casID]["cell_blue_n"],
                                                                            chek293.dAC50[casID]["cell_green_n"],
                                                                            chek293.dAC50[casID]["cell_red_n"])
        filout.write(linew)
    filout.close()

    return pfilout




def ChemByCurve(cassay, ppng, prout):


    if not "dresponse" in dir(cassay):
        cassay.responseCurves(drawn=0)

    for CASID in cassay.dresponse.keys():
        for condition in cassay.dresponse[CASID].keys():
            curveclass = cassay.dresponse[CASID][condition]['CURVE_CLASS2']
            AC50 = cassay.dresponse[CASID][condition]['AC50']
            prcurve = prout + condition + "/" + curveclass + "/"
            if not path.exists(prcurve):
                pathFolder.createFolder(prcurve)

            writeLine = ["CAS: " + str(CASID), "AC50: " + str(AC50), "Curve: " + str(curveclass)]
            pcaspng = ppng + CASID + ".png"
            if not path.exists(pcaspng):
                continue

            pimageout = prcurve + CASID + ".png"
            try:img = Image.open(pcaspng)
            except: continue
            imgnew = Image.new("RGBA", (580, 775), (250, 250, 250))
            imgnew.paste(img, (0, 0))
            draw = ImageDraw.Draw(imgnew)
            draw.text((10, 600), str(writeLine[0]), (0, 0, 0), font=font)
            draw.text((10, 625), str(writeLine[1]), (0, 0, 0), font=font)
            draw.text((10, 650), str(writeLine[2]), (0, 0, 0), font=font)
            imgnew.save(pimageout)

