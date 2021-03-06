from os import system, path, remove, chdir, getcwd, listdir
from re import search
from time import sleep

OPERA = "/home/borrela2/softwares/OPERA/OPERA_CLi_Linux/run_OPERA.sh"
MATLAB = "/usr/local/MATLAB/MATLAB_Runtime/v94"
PRSOURCE = "/home/borrela2/interference/sources/fluo/"
PADEL = "/home/borrela2/softwares/padel/PaDEL-Descriptor.jar"


def runRCMD(cmd, out = 0):

    chdir("./../Rscripts/")
    print cmd
    if out == 0:
        system(cmd)
        output = 0
    else:
        import subprocess
        output = subprocess.check_output(cmd, shell=True)
    chdir("./../fluo/")
    return output


def runRQSARModeling(cmd):

    workdir = getcwd()
    chdir("/home/borrela2/development/QSARPR/source/")
    print(cmd)
    system(cmd)
    chdir(workdir)

def babelConvertSDFtoSMILE(sdfread, clean_smi=1, rm_smi=1):

    if path.exists(sdfread):
        psmile = sdfread[:-3] + "smi"
        tempsdf = sdfread

    else:
        tempsdf = open("tempsdf.sdf", "w")
        tempsdf.write(sdfread)
        tempsdf.close()
        psmile = "tempsmile.smi"

    cmd_convert = "babel " + tempsdf + " " + psmile + " 2>/dev/null"
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


def crossResponseCurve(prresponse1, prresponse2, pAC501, pAC502, prout):

    cmd = "./crossResponseCurves.R " + str(prresponse1) + " " + str(prresponse2) + " " + str(pAC501) + " " + str(pAC502) + " " + str(prout)
    runRCMD(cmd)


def dataManager(pdesc, pAC50, corval, maxQauntile, nbNA, prout):


    pfiloutDesc = prout + "descClean.csv"
    pfiloutAff = prout + "AC50Clean.csv"
    if path.exists(pfiloutDesc): # control only desc because AC50 can be 0
        return [pfiloutDesc, pfiloutAff]
    else:
        cmd = "./preprocData.R " + str(pdesc) + " " + str(pAC50) + " " + str(corval) + " " + str(maxQauntile) + " 0 " + str(prout) + " " + str(nbNA)
        runRCMD(cmd)
        if path.exists(pfiloutDesc):
            return [pfiloutDesc, pfiloutAff]

    return 1


def prepDataQSAR(pdesc, pAC50, prout, valcor, maxquantile, splitratio, nbNA, logAff = "0", typeAff="All"):

    cmd = "./QSARsPrep.R " + pdesc + " " + pAC50 + " " + prout + " " + str(valcor) + " " + str(maxquantile) + " " + str(splitratio) + " " + str(logAff) + " " + typeAff + " " + str(nbNA)
    runRQSARModeling(cmd)



def QSARReg(ptrain, ptest, pcluster, prresult, nbCV):

    cmd = "./QSARsReg.R " + ptrain + " " + ptest + " " + pcluster + " " + prresult + " " + str(nbCV) + " >" + prresult + "perf.txt"
    runRQSARModeling(cmd)


def QSARClass(ptrain, ptest, prresult, nbCV):

    cmd = "./QSARsClass.R " + ptrain + " " + ptest + " 0 " + prresult + " " + str(nbCV) + " > " + prresult + "perf.txt"
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



def clusteringSimilarityMatrix(pmatrix, prresult, aggreg, clusteringMeth, optClusterMeth):

    pcluster = prresult + "cluster.csv"
    if path.exists(pcluster) and path.getsize(pcluster) > 10:
            return pcluster
    else:
        cmd = "./clusterFromSimilarity.R " + pmatrix + " " + str(clusteringMeth) + " " + str(aggreg) + " " + str(optClusterMeth) + " " + str(prresult)
        runRCMD(cmd)

    if path.exists(pcluster) and path.getsize(pcluster) > 10:
            return pcluster
    else:
        return 1




def drawEnrichSOM(pdesc1D2Dclean, pAC50, pmodel, prSOM):

    cmdSOM = "./SOMaps.R " + str(pdesc1D2Dclean) + " " + pAC50 + " " + str(pmodel) + " " + str(prSOM)
    runRCMD(cmdSOM)

    return


def TtestDesc(pdesc1D2Dclean, pAC50All, presult):

    cmd = "./comparisonActInact.R " + str(pdesc1D2Dclean) + " " + str(pAC50All) + " " + presult
    runRCMD(cmd)


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


def vennPlot(dAC50, pranalysis):
    cmd = "./VennDiagram.R " + str(dAC50) + " " + str(pranalysis)
    runRCMD(cmd)



def crossVenn(plucAC50, phepg2AC50, phek293AC50, prout):

    cmd = "./VennDiagramCross.R " + str(plucAC50) + " " + str(phepg2AC50) + " " + str(phek293AC50) + " " + str(prout)
    runRCMD(cmd)


def generateMainSOM(pdesc, prout, sizeMap, pSOMmodel):

    cmd = "./generateSOMModel.R " + pdesc + " " + prout + " " + str(sizeMap) + " " + pSOMmodel
    runRCMD(cmd)



def drawPCA(pdesc1D2Dclean, pAC50, prPCA):

    cmd = "./PCAAnalysis.R " + str(pdesc1D2Dclean) + " " + str(pAC50) + " " + str(prPCA)
    runRCMD(cmd)


def drawPCACross(pdesc1D2Dclean, pAC50_hepg2, pAC50_hek293, prCrossPCA):

    cmd = "./PCAAnalysisCross.R " + str(pdesc1D2Dclean) + " " + str(pAC50_hepg2) + " " + str(pAC50_hek293) + " " + str(prCrossPCA)
    runRCMD(cmd)


def drawPCACross2Set(pdesc1, paff, pdesc2, prout):

    cmd = "./PCACross2Set.R " + str(pdesc1) + " " + str(paff) + " " + str(pdesc2) + " " + str(prout)
    runRCMD(cmd)


def drawMDS(pdesc1D2Dclean, pAC50, prMDS):

    cmd = "./MDSAnalysis.R " + str(pdesc1D2Dclean) + " " + str(pAC50) + " " + str(prMDS)
    runRCMD(cmd)


def enrichmentCluster(pclusters, pallAC50, prCluster):

    cmd = "./enrichmentCluster.R " + str(pallAC50) + " " + str(pclusters) + " " + str(prCluster)
    runRCMD(cmd)

def enrichmentIndex(pdesc, prresult, pAC50All, methCluster, methDist, methAgg):

    if methDist == None:
        methDist = "None"

    cmd = "./enrichmentIndex.R " + pdesc + " " + pAC50All + " " + prresult + " " + methCluster + " " + methDist + " " + methAgg
    runRCMD(cmd)

def preciseEnrichmentIndex(pdesc, presult, pAC50All, ptable, methCluster, methDist, methAgg):

    if methDist == None:
        methDist = "None"

    cmd = "./enrichmentOptimal.R " + pdesc + " " + pAC50All + " " + ptable + " " + presult + " " + methCluster + " " \
          + methDist + " " + methAgg
    runRCMD(cmd)

def finalClustering(pmatrixIn, pAC50Full, pcluster, channel, distMeth, aggMeth, prtemp, verbose = 0):

    if verbose == 1:
        print pmatrixIn, ": matrix in"
        print pAC50Full, ": AC50 full"
        print pcluster, ": cluster"
        print channel, "chanell"
        print distMeth, "dist meth"
        print aggMeth, "aggregation"
        print prtemp, "path out"


    cmd = "./finalCluster.R " + pmatrixIn + " " + str(pAC50Full) + " " + pcluster + " " + channel + " " + str(distMeth) + " " + aggMeth + " " + prtemp
    runRCMD(cmd)


def findCluster(pdescAll, pdescChem, penrichment, pcluster, distMeth, aggMeth):

    cmd = "./findBestCluster.R " + str(pdescAll) + " " + str(pdescChem) + " " + str(penrichment) + " " + str(pcluster) + " " + str(distMeth) + " " + str(aggMeth)
    out = runRCMD(cmd, 1)

    out = out.strip().split(" ")[-1]
    return out

def generateCardResult(pfilout):

    cmd = "./cardResult.R " + str(pfilout)
    runRCMD(cmd)



def histGlobal(pAC50All, prhist):

    cmd = "./histAff.R " + str(pAC50All) + " " + str(prhist)
    runRCMD(cmd)


def visualizeActive(pdesc, pAC50, distmeth, aggType, prout):

    cmd = "./visuActive.R " + str(pdesc) + " " + str(pAC50) + " " + str(distmeth) + " " + str(aggType) + " " + str(prout)
    runRCMD(cmd)


def createImportanceTable(pmodel, ML, ptrain, prout):

    cmd = "./createImportanceTable.R " + str(pmodel) + " " + str(ML) + " " + ptrain + " " + prout
    runRQSARModeling(cmd)

    pimportance = prout + "ImportanceDesc"
    if path.exists(pimportance):
        return pimportance
    else:
        return 1


def runImportanceDesc(pimportance, nb):

    cmd = "./importancePlot.R " + str(pimportance) + " " + str(nb)
    runRQSARModeling(cmd)



def plotAC50VSEff(pfilin):

    cmd = "./plotEffVSAC50.R " + pfilin
    runRCMD(cmd)



def runPadel(prin=""):
    """Input include a folder of sdf file"""
    pfilout = prin + "tem.csv"
    if path.exists(pfilout) and path.getsize(pfilout) > 50:
        return pfilout

    if prin == "":
        return "ERROR - Padel Input"
    else:
        cmd = "java -jar " + PADEL + " -maxruntime 10000 -2d -dir " + str(prin) + " -file " + pfilout
        print cmd
        system(cmd)

    return pfilout




def runOPERA(psdf, p2Ddesc, prtemp):

    lfstart = [path.basename(psdf), "tem.csv"]

    if len (listdir(prtemp)) == 2:
        ppred = prtemp + path.basename(psdf)[:-3] + "csv"
        cmd = "%s %s -d %s -o %s -a -x" % (OPERA, MATLAB, p2Ddesc, ppred)
        print cmd
        chdir(prtemp)
        system(cmd)
        chdir(PRSOURCE)

    lfilout = []

    lfile = listdir(prtemp)
    for filpred in lfile:
        if not filpred in lfstart and not search("^\.", filpred):
            lfilout.append(prtemp + filpred)

    return lfilout


def plotAC50VSProb(pfilin, prout):

    cmd = "./qualityPred.R " + pfilin + " " + prout
    runRCMD(cmd)



def predictModel(pdesc, pmodelR, ML, proutbyRmodel):


    pfilout = proutbyRmodel + "perf_" + ML + "_" + pmodelR.split("/")[-1][:-6] + ".csv"
    cmd = "./predictfromModel.R " + pdesc + " " + pmodelR + " " + ML + " " + pfilout
    if not path.exists(pfilout):
        runRQSARModeling(cmd)

    if path.exists(pfilout):
        return pfilout
    else:
        return "ERROR"


def qualityPred(pfsum):

    pfilout = pfsum + "_perf.csv"
    if path.exists(pfilout):
        return pfilout

    cmd = "./perfProb.R " + pfsum + " " + pfilout
    runRQSARModeling(cmd)

    if path.exists(pfilout):
        return pfilout
    else:
        return "ERROR"


def vennPlotPubChem(pTox21Smiles, pTox21IC50, pPubChemSmile, pPubChemIC50, prout):

    cmd = "./VennDiagramTox21PubChem.R " + pTox21Smiles + " " + pTox21IC50 + " " + pPubChemSmile + " " + pPubChemIC50 + " " + prout
    runRCMD(cmd)

