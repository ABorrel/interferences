from urllib.request import urlretrieve
from os import path
from shutil import move
from re import search
from random import shuffle

import runExternalSoft
import toolbox
#import chemical


def formatPubChemTable(pfilin, PRPUBCHEM, prout, update = 0):

    pfilout = prout + "tableSmi.csv"
    if path.exists(pfilout) and update == 0:
        return pfilout
    else:
        filout = open(pfilout, "w")
        filout.write("ID\tSMILES\tActive\n")
        dchem = toolbox.loadMatrix(pfilin, sep=",")
        #print dchem.keys()

        for chemID in list(dchem.keys()):
            cpubmed = Chem(chemID, PRPUBCHEM)
            SMILES = cpubmed.getSMILE()
            if search("Error", SMILES):
                continue


            if chemID == "RESULT_IS_ACTIVE_CONCENTRATION" or chemID == "RESULT_UNIT":
                continue
            #print SMILES

            
            # add filter 
            #print dchem[chemID]["Log of MaxDyeEquivalency"], dchem[chemID]["PUBCHEM_ACTIVITY_OUTCOME"]
            if "Log of MaxDyeEquivalency" in list(dchem[chemID].keys()):
                if float(dchem[chemID]["Log of MaxDyeEquivalency"]) < -7:
                    dchem[chemID]["PUBCHEM_ACTIVITY_OUTCOME"] = "Inactive"
            filout.write("%s\t%s\t%s\n"%(chemID, SMILES, dchem[chemID]["PUBCHEM_ACTIVITY_OUTCOME"]))
        filout.close()
    return pfilout



def computeDesc(passay, PRDESC, PRSMI, prout, nbfile = 1, update = 0):


    # by pass
    pdescout = prout + "descMat"
    paff = prout + "aff.txt"
    if path.exists(pdescout) and update == 0 and nbfile == 1:
        return pdescout
    elif path.exists(pdescout) and update == 0 and nbfile == 2 and path.exists(paff):
        return [pdescout, paff]


    dchem = toolbox.loadMatrix(passay)
    lchemID = list(dchem.keys())
    try:lchemID.remove("RESULT_UNIT")
    except:pass
    shuffle(lchemID)
    i = 0
    nbi = len(lchemID)
    while i < nbi:
        if search("error", dchem[lchemID[i]]["SMILES"].lower()): # case of the table is computed before
            del dchem[lchemID[i]]
            del lchemID[i]
            nbi = nbi - 1
            continue

        if dchem[lchemID[i]]["Active"] == "Inconclusive" or search("Error", dchem[lchemID[i]]["SMILES"]):
            del dchem[lchemID[i]]
            del lchemID[i]
            nbi = nbi -1
            continue

        # compute descriptors
        cchem = chemical.chemical(lchemID[i], dchem[lchemID[i]]["SMILES"])
        cchem.prepareChem(PRSMI)
        if search("error", cchem.log.lower()):
            del dchem[lchemID[i]]
            del lchemID[i]
            nbi = nbi - 1
            continue

        cchem.compute1D2DDesc(PRDESC)
        if search("error", cchem.log.lower()):
            del dchem[lchemID[i]]
            del lchemID[i]
            nbi = nbi - 1
            continue



        cchem.computeOpera(update=update)
        if search("error", cchem.log.lower()):
            del dchem[lchemID[i]]
            del lchemID[i]
            nbi = nbi -1
            i = i - 1
            continue



        cchem.writeTablesDesc(PRDESC, update=update)

        i = i + 1

    if nbfile == 1:
        fildesc = open(pdescout, "w")
        ldesc = chemical.getLdesc("1D2D", 1) + chemical.getLdesc("Opera", 0)
        fildesc.write("ID," + ",".join(ldesc) + ",Aff" +"\n")

        for chemID in lchemID:
            print(chemID)
            if dchem[chemID]["Active"] == "Active":
                aff = 1
            else:
                aff = 0
            pdesc = PRDESC + chemID + ".txt"
            if path.exists(pdesc):
                ddesc = toolbox.loadMatrix(pdesc)
                lval = []
                for desc in ldesc:
                    if not desc in list(ddesc[chemID].keys()):
                        lval.append("NA")
                    else:
                        lval.append(str(ddesc[chemID][desc]))

                fildesc.write(chemID + "," + ",".join(lval) + "," + str(aff) + "\n")
        fildesc.close()
        return pdescout

    else:

        fildesc = open(pdescout, "w")
        paff = prout + "aff.txt"
        filaff = open(paff, "w")
        ldesc = chemical.getLdesc("1D2D", 1) + chemical.getLdesc("Opera", 0)
        fildesc.write("ID," + ",".join(ldesc) + "\n")
        filaff.write("ID\tAff\n")

        for chemID in lchemID:
            print(chemID)
            if dchem[chemID]["Active"] == "Active":
                aff = 1
            else:
                aff = 0
            pdesc = PRDESC + chemID + ".txt"
            if path.exists(pdesc):
                ddesc = toolbox.loadMatrix(pdesc)
                lval = []
                for desc in ldesc:
                    if not desc in list(ddesc[chemID].keys()):
                        lval.append("NA")
                    else:
                        lval.append(str(ddesc[chemID][desc]))

                fildesc.write(chemID + "," + ",".join(lval) + "," + str(aff) + "\n")
                filaff.write(chemID + "\t" + str(aff) + "\n")
        fildesc.close()
        filaff.close()

        return [pdescout, paff]







class Chem:
    def __init__(self, pubchemID, prbase):
        self.pubchemID = pubchemID
        self.prbase = prbase


    def getSDF(self):

        psdf = self.prbase + self.pubchemID + ".sdf"
        if path.exists(psdf):
            self.psdf = psdf

        else:
            request = (
                        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/%s/record/SDF/?version=5&record_type=2d&response_type=save&response_basename=Structure2D_SID_%s" % (
                    self.pubchemID, self.pubchemID))

            presult = urlretrieve(request)
            if path.getsize(presult[0]) < 1000:
                request = (
                        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/record/SDF/?record_type=2d&response_type=save&response_basename=Structure2D_CID_%s" % (
                    self.pubchemID, self.pubchemID))
                presult = urlretrieve(request)

            # move file
            if path.exists(presult[0]):
                move(presult[0], psdf)
                self.psdf = psdf
            else:
                print("Error: PubChemRequest")
                return "None"



    def getSMILE(self):

        if not "psdf" in self.__dict__:
            self.getSDF()


        SMILES = runExternalSoft.babelConvertSDFtoSMILE(self.psdf)
        print(SMILES)

        if SMILES != "":
            self.smi = SMILES
            return SMILES
        else:
            print("Error SMILES")
            return "Error: SMILES generation"




