from re import search
from os import listdir

import pathFolder
import runExternalSoft


class sdfDB:

    def __init__(self, psdf, kname, prout):
        self.psdf = psdf
        self.prout = prout
        self.name = kname


    def parseAll (self):

        lout = []

        filin = open(self.psdf, "r")
        handle_read = filin.read()
        filin.close()
        #print handle_read

        l_compound = handle_read.split("$$$$\n")

        for compound in l_compound:
            dcompound = {}
            dcompound["sdf"] = compound + "$$$$"
            llines = compound.split("\n")
            i = 0
            nblines = len(llines)
            while i < nblines:
                if search(">  <", llines[i]):
                    kin = llines[i].split(">  <")[1]
                    kin = kin.split(">")[0]
                    valuek = llines[i+1].strip()
                    dcompound[kin] = valuek
                i += 1
            if self.name in dcompound.keys():
                lout.append(dcompound)


        self.lc=lout
        print len(lout)

    def writeTable(self, filname):

        if not "lc" in dir(self):
            self.parseAll()

        ptable = self.prout + str(filname)
        filout = open(ptable, "w")

        lheader = [self.name]
        for compound in self.lc:
            for k in compound.keys():
                if not k in lheader:
                    lheader.append(k)

        lheader.remove("sdf")
        filout.write("\t".join(lheader) + "\n")

        for compound in self.lc:
            try: filout.write(compound[self.name])
            except:
                continue
            for h in lheader[1:]:
                if h in compound.keys():
                    filout.write("\t" + str(compound[h]))
                else:
                    filout.write("\tNA")
            filout.write("\n")
        filout.close()


    def renameHeader(self):
        if not "lc" in dir(self):
            self.parseAll()
        for compound in self.lc:
            sdfin = compound["sdf"]
            sdfin = sdfin.split("\n")
            sdfin[0] = compound[self.name]
            compound["sdf"] = "\n".join(sdfin)


    def splitSDF(self):

        if not "lc" in dir(self):
            self.parseAll()

        prSDF = self.prout + "cpdsdf/"
        pathFolder.createFolder(prSDF)
        self.prsdf = prSDF

        if len(self.lc) == len(listdir(prSDF)):
            return
        else:
            for compound in self.lc:
                pfilout = self.prsdf + compound[self.name] + ".sdf"
                filout = open(pfilout, "w")
                filout.write(compound["sdf"])
                filout.close()


    def drawMolecules(self, prpng):

        if not "prsdf" in dir(self):
            self.splitSDF()

        pathFolder.createFolder(prpng)
        self.prpng = prpng

        if len(listdir(self.prsdf)) == len(listdir(self.prpng)):
            return

        lfsdf = listdir(self.prsdf)
        for sdf in lfsdf:
            runExternalSoft.molconvert(self.prsdf + sdf, self.prpng + sdf[:-3] + "png")



    def writeTableSpecific(self, lkin, nametable):

        if not "lc" in dir(self):
            self.parseAll()

        pfilout = self.prout + nametable + ".csv"
        filout = open(pfilout, "w")

        # header table
        filout.write("ID" + "\t" + "\t".join(lkin) + "\n")

        lwritejs = []

        for compound in self.lc:
            try:
                filout.write(compound[self.name] + "\t")
                linejs = "\"" + compound[self.name] + "\":"
            except: continue
            lwrite = []
            for kin in lkin:
                if kin in compound.keys():
                    if kin == "SYNONYMS":
                        lwrite.append(str(compound[kin]).lower())
                    else:
                        lwrite.append(str(compound[kin]))
                else:
                    lwrite.append("NA")

            filout.write("\t".join(lwrite) + "\n")
            linejs = linejs + "[" + ",".join(["\"" + i + "\"" for i in lwrite]) + "]"
            lwritejs.append(linejs)

        filout.close()

