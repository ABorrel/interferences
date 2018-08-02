from PyPDF2 import PdfFileReader
from re import search, findall
from os import path

import toolbox

def parsepdf(prcytox, prresult):

    ppdf = prcytox + "toxsci-15-0719-File012.pdf"
    ptable = prcytox + "toxsci-15-0719-File009.csv"

    pfilout = prresult + "cytox.csv"
    if path.exists(pfilout):
        dout = toolbox.loadMatrix(pfilout, sep="\t")
        return dout



    dtable = toolbox.loadMatrix(ptable, sep=",")
    lCASID = []
    for chem in dtable.keys():
        lCASID.append(dtable[chem]["CASRN"])

    dout = {}
    fpdf = open(ppdf, "rb")
    pdfReader = PdfFileReader(fpdf)
    nbpage = pdfReader.getNumPages()
    #nbpage = 32

    i = 0
    while i < nbpage:
        print i
        pageObj = pdfReader.getPage(i)
        pageText = pageObj.extractText()
        llines = pageText.split("\n")

        for line in llines:
            for CASID in lCASID:
                if search(CASID, line):
                    dout[CASID] = {}
                    CAStemp = CASID
                    break
            if search("cytotox min=", line):
                cytoxmin = findall(r"[-+]?\d*\.\d+|\d+", line)
                #print cytoxmin
                dout[CAStemp]["CytoxMin"] = cytoxmin[0]
        i += 1

    filout = open(pfilout, "w")
    filout.write("CAS\tCytoxMin\n")
    for CASID in dout.keys():
        filout.write(str(CASID) + "\t" + str(dout[CASID]["CytoxMin"]) + "\n")
    filout.close()

    return dout

