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
            print(energy)
            lenergy.append(float(energy))

        # take minimal energy
        ibest = lenergy.index(min(lenergy))
        print(ibest)
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


def computeSimilarityFP(FP1, FP2, typeMetric ):
    from rdkit import DataStructs

    if typeMetric == 'Tanimoto':
        return DataStructs.FingerprintSimilarity(FP1, FP2, metric=DataStructs.TanimotoSimilarity)
    elif typeMetric == "Dice":
        return DataStructs.DiceSimilarity(FP1, FP2)
    elif typeMetric == "Cosine":
        return DataStructs.FingerprintSimilarity(FP1, FP2, metric=DataStructs.CosineSimilarity)
    elif typeMetric == "Sokal":
        return DataStructs.FingerprintSimilarity(FP1, FP2, metric=DataStructs.SokalSimilarity)
    elif typeMetric == "Russel":
        return DataStructs.FingerprintSimilarity(FP1, FP2, metric=DataStructs.RusselSimilarity)
    elif typeMetric == "RogotGoldberg":
        return DataStructs.FingerprintSimilarity(FP1, FP2, metric=DataStructs.RogotGoldbergSimilarity)
    elif typeMetric == "AllBit":
        return DataStructs.FingerprintSimilarity(FP1, FP2, metric=DataStructs.AllBitSimilarity)
    elif typeMetric == "Kulczynski":
        return DataStructs.FingerprintSimilarity(FP1, FP2, metric=DataStructs.KulczynskiSimilarity)
    elif typeMetric == "McConnaughey":
        return DataStructs.FingerprintSimilarity(FP1, FP2, metric=DataStructs.McConnaugheySimilarity)
    elif typeMetric == "Asymmetric":
        return DataStructs.FingerprintSimilarity(FP1, FP2, metric=DataStructs.AsymmetricSimilarity)
    elif typeMetric == "BraunBlanquet":
        return DataStructs.FingerprintSimilarity(FP1, FP2, metric=DataStructs.BraunBlanquetSimilarity)





def loadMatrix(pmatrixIn, sep = "\t"):

    filin = open(pmatrixIn, "r")
    llinesMat = filin.readlines()
    filin.close()

    dout = {}
    line0 = formatLine(llinesMat[0])
    line1 = formatLine(llinesMat[1])
    lheaders = line0.split(sep)
    lval1 = line1.split(sep)

    # case where R written
    if len(lheaders) == (len(lval1) -1):
        lheaders = ["ID"] + lheaders


    i = 0
    while i < len(lheaders):
        if lheaders[i] == "":
            lheaders[i] = "ID"
        i += 1

    i = 1
    imax = len(llinesMat)
    while i < imax:
        lineMat = formatLine(llinesMat[i])
        lvalues = lineMat.split(sep)
        kin = lvalues[0]
        dout[kin] = {}
        j = 0
        if len(lvalues) != len(lheaders):
            print(lvalues)
            print(lheaders)
        jmax = len(lheaders)
        while j < jmax:
            try:dout[kin][lheaders[j]] = lvalues[j]
            except:pass
            j += 1
        i += 1

    return dout


def loadMatrixToDict(pmatrixIn, sep ="\t"):

    filin = open(pmatrixIn, "r")
    llinesMat = filin.readlines()
    filin.close()

    dout = {}
    line0 = formatLine(llinesMat[0])
    line1 = formatLine(llinesMat[1])
    lheaders = line0.split(sep)
    lval1 = line1.split(sep)

    # case where R written
    if len(lheaders) == (len(lval1)-1):
        lheaders.append("val")

    i = 1
    imax = len(llinesMat)
    while i < imax:
        lineMat = formatLine(llinesMat[i])
        lvalues = lineMat.split(sep)
        kin = lvalues[0]
        dout[kin] = {}
        j = 0
        if len(lvalues) != len(lheaders):
            print(lineMat)
            print(llinesMat[i])
            print(lvalues)
            print("Error => nb element", i)
            print(len(lvalues))
            print(len(lheaders))
            ddd

        jmax = len(lheaders)
        while j < jmax:
            dout[kin][lheaders[j]] = lvalues[j]
            j += 1
        i += 1

    return dout




def formatLine(linein):

    linein = linein.replace("\n", "")
    linenew = ""

    imax = len(linein)
    i = 0
    flagchar = 0
    while i < imax:
        if linein[i] == '"' and flagchar == 0:
            flagchar = 1
        elif linein[i] == '"' and flagchar == 1:
            flagchar = 0

        if flagchar == 1 and linein[i] == ",":
            linenew = linenew + " "
        else:
            linenew = linenew + linein[i]
        i += 1

    linenew = linenew.replace('\"', "")
    return linenew



def writeMatrix(ddesc, pdescAct, sep = "\t"):


    filout = open(pdescAct, "w")
    lheader = list(ddesc[list(ddesc.keys())[0]].keys())

    # put header in first

    if "CAS" in lheader:
        del lheader[lheader.index("CAS")]
        lheader = ["CAS"] + lheader
    elif "CASID" in lheader:
        del lheader[lheader.index("CASID")]
        lheader = ["CASID"] + lheader
    else:
        lheader = ["CASID"] + lheader
        for casID in list(ddesc.keys()):
            ddesc[casID]["CASID"] = casID


    filout.write(sep.join(lheader) + "\n")
    for casID in list(ddesc.keys()):
        lval = [str(ddesc[casID][i]) for i in lheader]
        filout.write(sep.join(lval) + "\n")
    filout.close()



def loadSMILES(psmi):


    filin = open(psmi, "r")
    smiles = filin.readlines()[0]
    smiles.strip()
    filin.close()

    return smiles




