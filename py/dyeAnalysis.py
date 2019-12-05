from os import  path
import toolbox




class DYE:
    def __init__(self, pDYE, lcAssays, prout):

        self.prout = prout
        self.pDYE = pDYE
        self.lassays = lcAssays
        self.dDye = toolbox.loadMatrix(self.pDYE, sep = ",")


    def crossDyeAssays(self):


        for cassays in self.lassays:
            name = cassays.name.split("-")[2]
            print(name)
            lsample =  ["cell_blue_n", "cell_green_n", "cell_red_n", "med_blue_n", "med_green_n", "med_red_n"]

            pfilout = self.prout + "DYE_" + name + ".csv"
            filout = open(pfilout,'w')
            filout.write("CAS\tColor\t" + "\t".join(lsample) + "\n")

            for dye in list(self.dDye.keys()):
                print(self.dDye[dye]["casrn"])
                try: filout.write("%s\t%s\t%s\n"%(self.dDye[dye]["casrn"], self.dDye[dye]["color"], "\t".join([str(cassays.dAC50[self.dDye[dye]["casrn"]][sample]) for sample in lsample])))
                except: pass


            filout.close()

    def importDescriptors(self, prDesc = "/home/borrela2/interference/Desc/DESCbyCAS/"):

        ddesc = {}
        dcolor = {}
        for chemID in list(self.dDye.keys()):
            chem = self.dDye[chemID]
            CASID = chem["casrn"]
            color = chem["color"]

            pdescin = prDesc + CASID + ".txt"
            if path.exists(prDesc + CASID + ".txt"):
                dtemp = toolbox.loadMatrixToDict(pdescin)
                ddesc.update(dtemp)
                dcolor[CASID] = color

        self.dcolor = dcolor

        pdesc = self.prout + "descMat"
        fildesc = open(self.prout + "descMat", "w")
        ldesc = list(ddesc[list(ddesc.keys())[0]].keys())
        print(ldesc)
        print(CASID)
        del ldesc[ldesc.index("CAS")]
        fildesc.write("ID," + ",".join(ldesc) + ",Aff\n")
        for CASID in list(ddesc.keys()):
            lw = []
            for desc in ldesc:
                if desc in list(ddesc[CASID].keys()):
                    lw.append(str(ddesc[CASID][desc]))
                else:
                    lw.append("NA")
            fildesc.write("%s,%s,1\n" % (CASID, ",".join(lw)))
        fildesc.close()
        return pdesc