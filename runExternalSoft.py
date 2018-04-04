from os import system, path, remove
from re import search
from time import sleep

LIGPREP = "/opt/schrodinger2017-3/ligprep"
PADEL = "/home/aborrel/softwares/padel/PaDEL-Descriptor.jar"


def runLigprep(psmilin, forcefield="OPLS3", stereoisoster=1):

    """Maybe fix more option"""

    if forcefield == "OPLS3":
        bff = "16"
    else:
        bff = "14"

    cmd = LIGPREP + " -ismi " + psmilin + " -osd " + psmilin[0:-4] + ".sdf" + " -bff " + str(bff) + " -epik -s " + str(stereoisoster) + " -WAIT -NJOBS 3"

    print cmd
    system(cmd)

    # control if file exist
    if not path.exists(psmilin[0:-4] + ".sdf"):
        return "Ligprep ERROR"
    else:
        try:remove("tem.log")
        except:pass
        try:remove("tem-dropped.smi")
        except: pass
        try:remove("tem-dropped-indices.txt")
        except: pass

    return psmilin[0:-4] + ".sdf"

def runPadel(prin=""):
    """Input include a folder of sdf file"""
    if prin == "":
        return "ERROR - Padel Input"
    else:
        cmd = "java -jar " + PADEL + " -maxruntime 10000 -3d -dir " + str(prin) + " -file " + prin + "tem.desc"
        print cmd
        system(cmd)

    return prin + "tem.desc"


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

    print cmdStat
    system(cmdStat)

    return




def corplotR(pfilin):

    cmdCorplot = "./corplot.R " + str(pfilin)
    print cmdCorplot
    system(cmdCorplot)



def plotAC50(pAC50, prout, typeAssays):

    cmdhist = "./distributions.R " + str(pAC50) + " " + str(prout) + " " + str(typeAssays)
    print cmdhist
    system(cmdhist)



def molconvert(pfilin, pfilout= ""):
    """Convert with black background"""
    if pfilout == "":
        pfilout = pfilin[:-3] + "png"

    if path.exists(pfilout):
        return pfilout
    cmdconvert = "molconvert \"png:w500,Q100,#00000000\" " + pfilin + " -o " + pfilout
    system(cmdconvert)
    return pfilout



def plotResponsiveCurve(prresponse, pAC50, prout):

    cmd = "./responseCurves.R " + str(prresponse) + " " + str(pAC50) + " " + str(prout)
    print cmd
    system(cmd)


