import pathFolder
import chemical

from os import path




class predictionInterference:


    def __init__(self, cDB, prclusters, prout):

        self.prout = pathFolder.createFolder(prout)
        self.loadCluster(prclusters)
        self.cDB = cDB



    def loadCluster(self):




        return





    def predictlpSMI(self, lpsmi):

        for psmi in lpsmi:
            nameChemical = psmi.split("/")[0][0:-4]
            filin = open(psmi, "r")
            smiles = filin.readlines()[0]
            smiles.strip()
            print smiles
            filin.close()

            self.predictSMI(nameChemical, smiles)




    def predictSMI(self, nameChemical, smiles):


        prresult = pathFolder.createFolder(self.prout + nameChemical + "/")
        chem = chemical.chemical(nameChemical, smiles)
        chem.prepareChem(prresult)
        chem.compute1D2DDesc(prresult)
        chem.writeTablesDesc(prresult)
        chem.computeFP(typeFP="All")





        return