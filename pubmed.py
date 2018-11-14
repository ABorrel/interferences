from urllib import urlretrieve
from os import path
from shutil import move

import runExternalSoft



class Chem:
    def __init__(self, pubchemID, prbase):
        self.pubchemID = pubchemID
        self.prbase = prbase


    def getSDF(self):

        psdf = self.prbase + self.pubchemID + ".sdf"
        if path.exists(psdf):
            self.psdf = psdf

        else:
            request = (
                        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/%s/record/SDF/?version=5&record_type=2d&response_type=save&response_basename=Structure2D_SID_%s" % (
                    self.pubchemID, self.pubchemID))

            presult = urlretrieve(request)
            if path.getsize(presult[0]) < 1000:
                request = (
                        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/record/SDF/?record_type=2d&response_type=save&response_basename=Structure2D_CID_%s" % (
                    self.pubchemID, self.pubchemID))
                presult = urlretrieve(request)

            # move file
            if path.exists(presult[0]):
                move(presult[0], psdf)
                self.psdf = psdf
            else:
                print "Error: PubChemRequest"
                return "None"



    def getSMILE(self):

        if not "psdf" in self.__dict__:
            self.getSDF()


        SMILES = runExternalSoft.babelConvertSDFtoSMILE(self.psdf)
        print SMILES

        if SMILES != "":
            self.smi = SMILES
            return SMILES
        else:
            print "Error SMILES"
            return "Error: SMILES generation"




