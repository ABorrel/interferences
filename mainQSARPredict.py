import toolbox
import QSARpredictor
import pubmed
import pathFolder
import chemical

from os import path, listdir
from re import search


PRMAIN = "/home/borrela2/interference/"
PRTESTING = "/home/borrela2/interference/testing/"
PRMODELS = PRTESTING + "QSARmodel/"
PRPUBCHEM = PRMAIN + "PUBCHEM/"

PRDESC = PRTESTING + "DESC/"
PRSMI = PRTESTING + "SMI/"



def formatPubChemTable(pfilin, prout):

    pfilout = prout + "tableSmi.csv"
    if path.exists(pfilout):
        return pfilout
    else:
        filout = open(pfilout, "w")
        filout.write("ID\tSMILES\tActive\n")
        dchem = toolbox.loadMatrix(pfilin, sep=",")
        #print dchem.keys()

        for chemID in dchem.keys():
            cpubmed = pubmed.Chem(chemID, PRPUBCHEM)
            SMILES = cpubmed.getSMILE()

            if chemID == "RESULT_IS_ACTIVE_CONCENTRATION" or chemID == "RESULT_UNIT":
                continue
            print SMILES

            filout.write("%s\t%s\t%s\n"%(chemID, SMILES, dchem[chemID]["PUBCHEM_ACTIVITY_OUTCOME"]))

        filout.close()

    return pfilout




def computeDesc(passay, PRDESC, PRSMI):


    dchem = toolbox.loadMatrix(passay)
    lchemID = dchem.keys()
    for chemID in dchem.keys():
        print dchem[chemID]
        cchem = chemical.chemical(chemID, dchem[chemID]["SMILES"])
        cchem.prepareChem(PRSMI)
        cchem.compute1D2DDesc(PRDESC)
        cchem.computeOpera()
        cchem.writeTablesDesc(PRDESC)

        ddd

        if search("error", cchem.log.lower()):
            del lchemID[lchemID.index(chemID)]


    ldesc = chemical.getLdesc()








for assay in listdir("/home/borrela2/interference/testing/data/"):

    dPubChemBioassays = "/home/borrela2/interference/testing/data/" + assay
    prout = pathFolder.createFolder(PRTESTING + assay[4:-4] + "/")
    passay = formatPubChemTable(dPubChemBioassays, prout)

    computeDesc(passay, PRDESC, PRSMI)











