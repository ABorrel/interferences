





def loadTableAssays(pfilin):

    lout = []
    filin = open(pfilin, "r")
    llines = filin.readlines()
    filin.close()

    lheader =llines[0].strip().split("\t")

    for chemical in llines[1:]:
        lval = chemical.strip().split("\t")
        dtemp={}
        i = 0
        while i < len(lval):
            dtemp[lheader[i]] = lval[i]
            i += 1
        lout.append(dtemp)

    return lout



