
import toolbox




class DYE:
    def __init__(self, pDYE, lcAssays, prout):

        self.prout = prout
        self.pDYE = pDYE
        self.lassays = lcAssays


    def crossDyeAssays(self):


        dDYE = toolbox.loadMatrix(self.pDYE, sep = ",")
        print dDYE

        for cassays in self.lassays:
            name = cassays.name.split("-")[2]
            print name
            lsample =  ["cell_blue_n", "cell_green_n", "cell_red_n", "med_blue_n", "med_green_n", "med_red_n"]

            pfilout = self.prout + "DYE_" + name + ".csv"
            filout = open(pfilout,'w')
            filout.write("CAS\tColor\t" + "\t".join(lsample) + "\n")

            for dye in dDYE.keys():
                print dDYE[dye]["casrn"]
                try: filout.write("%s\t%s\t%s\n"%(dDYE[dye]["casrn"], dDYE[dye]["color"], "\t".join([str(cassays.dAC50[dDYE[dye]["casrn"]][sample]) for sample in lsample])))
                except: pass


            filout.close()