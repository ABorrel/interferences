from os import path

import assayResults
import pathFolder
import liganddescriptors
import loadDB



def main(pfilin, pSDFToxCast, prresult):

    #1 load data
    lchem = assayResults.loadTableAssays(pfilin)
    print len(lchem)


    #2 control we have chemical
    db = loadDB.sdfDB(pSDFToxCast, "Substance_CASRN", prresult + "dbToxCas/")
    db.parseAll()

    lno = []
    lcas = []
    for chem in lchem:
        #print chem.keys()
        cas = chem["CAS"]
        lcas.append(cas)

        flag = 0
        for cpdDB in db.lc:
            if cpdDB["Substance_CASRN"] == cas:
                flag = 1
                break
        if flag == 0:
            lno.append(cas)
        else:
            if not cas in lcas:
                lcas.append(cas)

    # number not found in toxcast
    print len(lno)
    print lno


    # compute desc 1D - 2D


    #3 compute desc RDKIT

    liganddescriptors.





    #4 R clustering





#########
# MAIN  #
#########


pSDFToxCast= "/home/borrela2/ToxCast_release20151019/DSSTox_ToxCastRelease_20151019.sdf"

pluc = "/home/borrela2/interference/data/luc/tox21-luc-biochem-p1/tox21-luc-biochem-p1.aggregrated.txt"
phek293 = "/home/borrela2/interference/data/luc/tox21-spec-hek293-p1/tox21-spec-hek293-p1.txt"
phepg2 = "/home/borrela2/interference/data/luc/tox21-spec-hepg2-p1/tox21-spec-hepg2-p1.txt"



# luciferase #
##############
prresultLuc = "/home/borrela2/interference/luc-biochem/"
prresultLuc = pathFolder.createFolder(prresultLuc)

main(pluc, pSDFToxCast, prresultLuc)


# hek293 #
##########

prresultHek293 = "/home/borrela2/interference/hek293/"
prresultHek293 = pathFolder.createFolder(prresultHek293)

#main(phek293, pSDFToxCast, prresultHek293)




