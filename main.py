from idlelib import run
from os import path

import assayResults
import pathFolder
import liganddescriptors
import loadDB
import runExternalSoft


def main(pfilin, pSDFTox, lheader, pTox, prresult, desc1D2D = 1):

    #1 load data
    cchem = assayResults.assays(pfilin)

    cchem.writeTable(lheader, pTox)

    #2 control we have chemical
    db = loadDB.sdfDB(pSDFTox, "CASRN", prresult + "dbToxCas/")
    db.parseAll()



    lno = []
    dcas = {}
    for chem in cchem.lchem:
        #print chem.keys()
        cas = chem["CAS"]

        flag = 0
        for cpdDB in db.lc:
            if cpdDB["CASRN"] == cas:
                flag = 1
                break
        if flag == 0:
            lno.append(cas)
        else:
            if not cas in dcas.keys():
                dcas[cas] = {}
                dcas[cas]["db"] = cpdDB
                dcas[cas]["hts"] = chem

    # number not found in toxcast
    print lno

    #3 compute desc RDKIT
    prdesc = prresult + "DESC/"
    pathFolder.createFolder(prdesc)


    ddesc = {}
    ddesc["1D2D"] = prdesc + "desc.csv"

    if desc1D2D == 1:
        plog = prdesc + "log.txt"
        flog = open(plog, "w")
        for cas in dcas.keys():
            #dcas[cas]["db"]["SMILES"] = dcas[cas]["db"]["Structure_SMILES"]# to format
            print cas
            desc = liganddescriptors.Descriptors(dcas[cas]["db"], flog, prdesc, ddesc, namek="CASRN")
            desc.prepareChemical()
            if desc.log == "ERROR":
                continue
            else:
                desc.get_descriptor1D2D()
                desc.writeTablesDesc("1D2D")
        flog.close()


    #4 R analysis
    prAnalysis = prresult + "statAnalysis/"
    pathFolder.createFolder(prAnalysis)

    runExternalSoft.Rstat(ddesc["1D2D"], pTox, prAnalysis, "0.9", "80")



def mainCor(passay1, passay2, passay3, prout):

    # 1 load data
    cchem1 = assayResults.assays(passay1, prout)
    cchem2 = assayResults.assays(passay2, prout)
    cchem3 = assayResults.assays(passay3, prout)

    cchem1.corAssay(cchem3)
    cchem1.corAssay(cchem2)
    cchem2.corAssay(cchem3)



def responseCurves(passay, prout):

    prresponse = prout + "responseCurve/"
    pathFolder.createFolder(prresponse)

    #load table
    cchem = assayResults.assays(passay, prout)

    dresponse = {}
    for chem in cchem.lchem:
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
        filout.write("CONC\tDATA\tFluorophores\n")

        i = 0
        while i < 16:
            for sample in dresponse[CASID].keys():
                filout.write(str(dresponse[CASID][sample]["CONC"][i]) + "\t" + str(dresponse[CASID][sample]["DATA"][i]) + "\t" + str(sample) + "\n")
            i += 1
        filout.close()

    pAC50 = cchem.writeAC50()
    runExternalSoft.plotResponsiveCurve(prresponse, pAC50, prout)



def AC50Distribution(passay, prout):


    prIC50 = prout + "histIC50/"
    pathFolder.createFolder(prIC50)

    # load table
    cchem = assayResults.assays(passay, prout)
    pAC50 = cchem.writeAC50()

    runExternalSoft.plotAC50(pAC50, prIC50, cchem.name.split("-")[1])





#########
# MAIN  #
#########


pSDFToxCast= "/home/borrela2/ToxCast_release20151019/DSSTox_ToxCastRelease_20151019.sdf"
pSDFTox21 = "/home/borrela2/Tox21/TOX21SL.sdf"

pluc = "/home/borrela2/interference/data/luc/tox21-luc-biochem-p1/tox21-luc-biochem-p1.txt"
phek293 = "/home/borrela2/interference/data/luc/tox21-spec-hek293-p1/tox21-spec-hek293-p1.txt"
phepg2 = "/home/borrela2/interference/data/luc/tox21-spec-hepg2-p1/tox21-spec-hepg2-p1.txt"

prresults = "/home/borrela2/interference/spDataAnalysis/"

# folder specific for assays
prSphek293 = pathFolder.createFolder(prresults + phek293.split("/")[-1].split(".")[0] + "/")
prSpluc = pathFolder.createFolder(prresults + pluc.split("/")[-1].split(".")[0] + "/")
prSphepg2 = pathFolder.createFolder(prresults + phepg2.split("/")[-1].split(".")[0] + "/")


# plot correlation # -> not used
####################
#mainCor(pluc, phek293, phepg2, prresults)


# plot response curves #
########################
#responseCurves(pluc, prSpluc)
#responseCurves(phek293, prSphek293)
#responseCurves(phepg2, prSphepg2)



# IC50 hist #
#############
#AC50Distribution(pluc, prSpluc)
#AC50Distribution(phepg2, prSphepg2)
#AC50Distribution(phek293, prSphek293)



####### MOLECULAR DESCRIPTOR #
##############################

# luciferase #
##############
prresultLuc = "/home/borrela2/interference/luc-biochem/"
prresultLuc = pathFolder.createFolder(prresultLuc)

#main(pluc, pSDFTox21, ["CAS", "AC50", "ASSAY_OUTCOME"], prresultLuc + "AC50.csv", prresultLuc)




# hek293 #
##########
prresultHek293 = "/home/borrela2/interference/hek293/"
prresultHek293 = pathFolder.createFolder(prresultHek293)

#main(phek293, pSDFTox21, ["CAS", "AC50", "ASSAY_OUTCOME"], prresultHek293 + "AC50.csv", prresultHek293)


# hepg2 #
#########
prresultHepg2 = "/home/borrela2/interference/hepg2/"
prresultHepg2 = pathFolder.createFolder(prresultHepg2)

#main(phepg2, pSDFTox21, ["CAS", "AC50", "ASSAY_OUTCOME"], prresultHepg2 + "AC50.csv", prresultHepg2)





