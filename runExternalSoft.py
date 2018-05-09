from os import system, path, remove, chdir, getcwd
from re import search
from time import sleep


def runRCMD(cmd):

    chdir("./../Rscripts/")
    print cmd
    system(cmd)
    chdir("./../fluo/")


def runRQSARModeling(cmd):

    workdir = getcwd()
    chdir("/home/borrela2/QSARPR/source/")
    print(cmd)
    system(cmd)
    chdir(workdir)

def babelConvertSDFtoSMILE(sdfread, clean_smi=0, rm_smi=1):

    tempsdf = open("tempsdf.sdf", "w")
    tempsdf.write(sdfread)
    tempsdf.close()

    psmile = "tempsmile.smi"

    cmd_convert = "babel tempsdf.sdf " + psmile + " 2>/dev/null"
    system(cmd_convert)

    try : filin = open (psmile, "r")
    except : return "0"
    l_Fline = filin.readlines ()
    filin.close ()
    try : smile = l_Fline[0].split ("\t")[0]
    except : return "0"

    # rewrite path in filout
    if clean_smi == 1:
        filout = open (psmile, "w")
        filout.write (str (smile))
        filout.close ()

    if rm_smi == 1:
        system("rm " + psmile)


    return smile



def Rstat(pfilin1D2D, paff, prout, valcor = 0.9, maxquantile=80):

    cmdStat = "./descAnalysis.R " + str(pfilin1D2D) + " " + str(paff) + " " + str(prout) + " " + str(valcor) + " " + str(maxquantile)
    runRCMD(cmdStat)

    return




def corplotR(pfilin):

    cmdCorplot = "./corplot.R " + str(pfilin)
    runRCMD(cmdCorplot)



def plotAC50(pAC50, prout, typeAssays):

    cmdhist = "./distributions.R " + str(pAC50) + " " + str(prout) + " " + str(typeAssays)
    runRCMD(cmdhist)



def molconvert(pfilin, pfilout= ""):
    """Convert with black background"""
    if pfilout == "":
        pfilout = pfilin[:-3] + "png"

    if path.exists(pfilout):
        return pfilout
    #cmdconvert = "molconvert \"png:w500,Q100,#00000000\" " + pfilin + " -o " + pfilout  # for transparent background
    cmdconvert = "molconvert \"png:w500,Q100\" " + pfilin + " -o " + pfilout
    system(cmdconvert)
    return pfilout



def plotResponsiveCurve(prresponse, pAC50, prout):

    cmd = "./responseCurves.R " + str(prresponse) + " " + str(pAC50) + " " + str(prout)
    runRCMD(cmd)



def dataManager(pdesc, pAC50, corval, maxQauntile, prout):

    cmd = "./preprocData.R " + str(pdesc) + " " + str(pAC50) + " " + str(corval) + " " + str(maxQauntile) + " 1 " + str(prout)
    runRCMD(cmd)

    pfilout = prout + "DescClean.txt"

    if path.exists(pfilout):
        return pfilout
    else:
        return 0


def prepDataQSAR(pdesc, pAC50, prout, valcor, maxquantile, splitratio, typeAff="All"):

    cmd = "./QSARsPrep.R " + pdesc + " " + pAC50 + " " + prout + " " + str(valcor) + " " + str(maxquantile) + " " + str(splitratio) + " 1 " + typeAff
    runRQSARModeling(cmd)



def QSARReg(ptrain, ptest, pcluster, prresult, nbCV):

    cmd = "./QSARsReg.R " + ptrain + " " + ptest + " " + pcluster + " " + prresult + " " + str(nbCV) + " >" + prresult + "RegResults"
    runRQSARModeling(cmd)


def clustering(pdesc, pAC50, prresult, dist, aggreg, clusteringMeth, optClusterMeth):


    pcluster = prresult + "cluster.csv"
    if path.exists(pcluster):
        if pAC50 == "0":
            return pcluster
        else:
            cmd = "./clusterAnalysis.R " + pdesc + " " + pAC50 + " " + pcluster + " " + prresult + " " + dist + " " + clusteringMeth + " " + aggreg + " " + optClusterMeth
    else:
        cmd = "./clusterAnalysis.R " + pdesc + " " + pAC50 + " 0 " + prresult + " " + dist + " " + clusteringMeth + " " + aggreg + " " + optClusterMeth

    runRCMD(cmd)

    if path.exists(pcluster):
        return pcluster
    else:
        return 0


def drawEnrichSOM(pdesc1D2Dclean, pAC50, prSOM):


    cmdSOM = "./SOMaps.R " + str(pdesc1D2Dclean) + " " + pAC50 + " " + str(prSOM)
    runRCMD(cmdSOM)

    return



def CrossClusterIC50(pdesc, pAC50, pclust, prout):


    cmd = "./clusterIC50.R " + pdesc + " " + pAC50 + " " + pclust + " " + str(prout)
    runRCMD(cmd)



def crossA50s(pdesc, pAC50luc, pAC50hepg, pAC50hek, prout):

    cmd = "./crossAC50bycluster.R " + pdesc + " " + pAC50luc + " " + pAC50hepg + " " + pAC50hek + " " + prout
    runRCMD(cmd)



def corAC50(pAC50, pcurve, prresult):
    cmd = "./corAC50.R " + str(pAC50) + " " + str(pcurve) + " " + prresult
    runRCMD(cmd)


def barplotCurve(pfilin):
    cmd = "./barplotCurve.R " + str(pfilin)
    runRCMD(cmd)