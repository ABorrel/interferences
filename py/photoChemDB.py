from os import path

import toolbox



class photoChem():


    def __init__(self, pDB, lassays, prout):
        self.pDB = pDB
        self.DB = toolbox.loadMatrixToDict(self.pDB, sep="\t")
        self.prout = prout
        self.lassays = lassays


    def crossSpectrumAssays(self):

        for cassays in self.lassays:
            name = cassays.name.split("-")[2]
            print name
            lsample = ["cell_blue_n", "cell_green_n", "cell_red_n", "med_blue_n", "med_green_n", "med_red_n"]

            pfilout = self.prout + "Spectrum_" + name + ".csv"
            filout = open(pfilout, 'w')
            filout.write("CAS\tWavelength\t" + "\t".join(lsample) + "\n")

            lwave = []
            llw = []
            for chemID in self.DB.keys():
                chem = self.DB[chemID]
                CASID = chem["Structure"].split("_")[1]
                Abs = chem["Wavelength"]
                try:
                    lw = "%s\t%s\t%s\n" % (CASID, Abs, "\t".join([str(cassays.dAC50[CASID][sample]) for sample in lsample]))
                    llw.append(lw)
                    lwave.append(Abs)
                except:
                    pass

            lwave.sort()
            for Abs in lwave:
                i = 0
                imax = len(llw)
                while i < imax:
                    if float(llw[i].split("\t")[1]) == float(Abs):
                        filout.write(llw[i])
                        del llw[i]
                        imax = imax - 1
                    else:
                        i = i + 1

            filout.close()


    def importDescriptors(self, prDesc = "/home/borrela2/interference/Desc/DESCbyCAS/"):

        ddesc = {}
        dwave = {}
        for chemID in self.DB.keys():
            chem = self.DB[chemID]
            CASID = chem["Structure"].split("_")[1]
            Abs = chem["Wavelength"]

            pdescin = prDesc + CASID + ".txt"
            if path.exists(prDesc + CASID + ".txt"):
                dtemp = toolbox.loadMatrixToDict(pdescin)
                ddesc.update(dtemp)
                dwave[CASID] = Abs

        self.dwave = dwave

        pdesc = self.prout + "descMat"
        fildesc = open(self.prout + "descMat", "w")
        ldesc = ddesc[ddesc.keys()[0]].keys()
        print ldesc
        print CASID
        del ldesc[ldesc.index("CAS")]
        fildesc.write("ID," + ",".join(ldesc) + ",Aff\n")
        for CASID in ddesc.keys():
            lw = []
            for desc in ldesc:
                if desc in ddesc[CASID].keys():
                    lw.append(str(ddesc[CASID][desc]))
                else:
                    lw.append("NA")
            fildesc.write("%s,%s,1\n"%(CASID, ",".join(lw)))
        fildesc.close()
        return pdesc


