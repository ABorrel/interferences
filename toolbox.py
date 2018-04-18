from os import system, path
from shutil import copy
from copy import deepcopy

def selectMinimalEnergyLigPrep(psdfin, psdfout):

    # case of only one
    filin = open(psdfin, "r")
    readfile = filin.read()
    filin.close()

    lsdf = readfile.split("$$$$\n")[:-1]


    if len(lsdf) == 1:
        copy(psdfin, psdfout)

    else:
        #find with the lower energy
        lenergy = []
        for sdfin in lsdf:
            energy = sdfin.split("> <r_lp_Energy>\n")[-1].split("\n")[0]
            print energy
            lenergy.append(float(energy))

        # take minimal energy
        ibest = lenergy.index(min(lenergy))
        print ibest
        filout = open(psdfout, "w")
        filout.write(lsdf[ibest] + "$$$$\n")
        filout.close()

    return psdfout



def renameHeaderSDF(pfilin):
    """Rename header with name file"""
    namesdf = pfilin.split("/")[-1].split(".")[0]
    filin = open(pfilin, "r")
    llines = filin.readlines()
    filin.close()
    llines[0] = str(namesdf) + "\n"

    filout = open(pfilin, "w")
    filout.write("".join(llines))
    filout.close()

import multiprocessing
import time

def timeFunction(funct, mol):

    manager = multiprocessing.Manager()
    lout = manager.list()

    p = multiprocessing.Process(target=funct, args=(mol, lout))
    p.start()
    time.sleep(2)

    if p.is_alive():
        p.terminate()
        p.join()
        return "ERROR"
    else:
        p.join()
        #print lout
        return lout[0]


from scipy import stats
from numpy import delete
def rankList(lin):

    liNA =  [i for i,x in enumerate(lin) if x == 'NA']
    linWithoutNA = deepcopy(lin)
    linWithoutNA = delete(linWithoutNA, liNA).tolist()


    linWithoutNA = [float(i) for i in linWithoutNA]
    n = len(linWithoutNA)

    lrank = stats.rankdata(linWithoutNA)
    lrank = [int(abs(i-n)) for i in lrank]

    for i in liNA:
        lrank.insert(i, "NA")

    return lrank
















