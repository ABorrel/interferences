from os import path, listdir
from pdb import run
from re import search
from shutil import copyfile

import runExternalSoft
import pathFolder




class clustering:

    def __init__(self, pdesc, prout, prPNG, corval, maxQuantile):


        self.prout = prout
        self.corval = corval
        self.maxquantile = maxQuantile
        self.pdesc = pdesc
        self.prPNG = prPNG




    def createMainClustering(self, disttype = "euc", aggregtype = "ward.D2", clusterType = "hclust", optimalCluster = "gap_stat"):

        self.distmeth = disttype
        self.aggType = aggregtype
        self.clusterMeth = clusterType
        self.optimalNBclustMeth = optimalCluster

        self.prCluster = self.prout + str(self.clusterMeth) + "_" + str(self.distmeth) + "_" + str(aggregtype.replace(".", "")) + "_" + str(self.optimalNBclustMeth) + "/"
        pathFolder.createFolder(self.prCluster)


        # data preparation
        self.pdesclean = self.prCluster + "descClean.csv"

        if not path.exists(self.pdesclean):
            if path.exists(self.pdesc) and path.getsize(self.pdesc) > 10:
                # preproc
                runExternalSoft.dataManager(self.pdesc, 0, self.corval, self.maxquantile, self.prCluster)
            else:
                print "Error ->", self.pdesc


        pcluster = self.prCluster + "cluster.csv"
        if not path.exists(pcluster):
            #clustering -> first level
            pcluster = runExternalSoft.clustering(self.pdesclean, "0", self.prCluster, self.distmeth, self.aggType, self.clusterMeth, self.optimalNBclustMeth)

        #Clustering second level
        if pcluster != 0:
            self.createSecondaryClustering(pcluster)


        # create main cluster file
        pclustersFinal = self.prCluster + "clusterMain.csv"
        if not path.exists(pclustersFinal):
            fclustersFinal = open(pclustersFinal, "w")
            fclustersFinal.write("ID,Cluster1,Cluster2\n")

            fcluster1 = open(pcluster, "r")
            lchemCluster1 = fcluster1.readlines()
            fcluster1.close()

            dclust = {}
            for chemCluster1 in lchemCluster1[1:]:
                chemCluster1 = chemCluster1.strip().replace("\"", "").split(",")
                chemID = chemCluster1[0]
                clust = chemCluster1[1]
                dclust[chemID] = [clust]

            for fileCluster in listdir(self.prCluster):
                if search("Clust", fileCluster):
                    pclust2 = self.prCluster + fileCluster + "/cluster.csv"
                    if path.exists(pclust2):
                        fclust2 = open(pclust2, "r")
                        lchemCluster2 = fclust2.readlines()
                        fclust2.close()

                        for chemCluster2 in lchemCluster2[1:]:
                            chemCluster2 = chemCluster2.strip().replace("\"", "").split(",")
                            chemID = chemCluster2[0]
                            clust2 = chemCluster2[1]

                            dclust[chemID].append(clust2)

            #write main cluster
            for chemID in dclust.keys():
                if len(dclust[chemID]) == 1:
                    dclust[chemID].append("1")
                fclustersFinal.write(str(chemID) + "," + ",".join(dclust[chemID]) + "\n")
            fclustersFinal.close()

        self.pclusters = pclustersFinal


    def createSecondaryClustering(self, pClusters):

        fcluster = open(pClusters, "r")
        lchemicals = fcluster.readlines()
        fcluster.close()

        dclust = {}
        for chemical in lchemicals[1:]:
            chemical = chemical.strip().replace("\"", "")
            chemical = chemical.split(",")
            ID = chemical[0]
            cluster = chemical[1]

            if not cluster in dclust.keys():
                dclust[cluster] = []
            dclust[cluster].append(ID)


        fdesc = open(self.pdesclean, "r")
        lchemdesc = fdesc.readlines()
        fdesc.close()

        ddesc = {}
        for chemdesc in lchemdesc[1:]:
            ID = chemdesc.split(",")[0].replace("\"", "")
            ddesc[ID] = chemdesc


        #write cluster and chemical
        for cluster in dclust.keys():
            prcluster = self.prCluster + "Clust" + str(cluster) + "/"
            if not path.exists(prcluster + "cluster.csv"):
                pathFolder.createFolder(prcluster)
                pdesc = prcluster + "descClean.csv"
                fdesc = open(pdesc, "w")
                fdesc.write(lchemdesc[0])

                for chemID in dclust[cluster]:
                    fdesc.write(ddesc[chemID])
                fdesc.close()

                runExternalSoft.clustering(pdesc, "0", prcluster, self.distmeth, self.aggType, self.clusterMeth, self.optimalNBclustMeth)

    def applyMainClusters(self, pAC50, prout):

        prclusterApplied = prout + self.prCluster.split("/")[-2] + "/"
        pathFolder.createFolder(prclusterApplied)

        # first level of cluster
        runExternalSoft.CrossClusterIC50(self.pdesclean, pAC50, self.pclusters, prclusterApplied)

        dclust = {}
        fcluster = open(self.pclusters, "r")
        lchem = fcluster.readlines()
        fcluster.close()

        for chem in lchem[1:]:
            chem = chem.strip().replace("\"", "").split(",")
            chemID = chem[0]
            cluster = str(chem[1])

            if not cluster in dclust.keys():
                dclust[cluster] = []
            dclust[cluster].append(chemID)
        print dclust

        prclusterSub = prclusterApplied + "clusterSub/"
        pathFolder.createFolder(prclusterSub)

        for clust in dclust.keys():
            print clust
            prbyclust = prclusterSub + "clust" + str(clust) + "/"
            pathFolder.createFolder(prbyclust)

            #file to copy
            pdescsub = prbyclust + "descClean.csv"
            pclustsub = prbyclust + "cluster.csv"

            copyfile(self.prCluster + "Clust" + str(clust) + "/descClean.csv", pdescsub)
            if not path.exists(self.prCluster + "Clust" + str(clust) + "/cluster.csv"):
                fclustsub = open(pclustsub, "w")
                fclustsub.write('\"ID\",\"cluster\"\n')
                for chem in dclust[clust]: fclustsub.write('\"' + str(chem) + '\"' + "," + '\"1\"\n')
                fclustsub.close()
            else:
                copyfile(self.prCluster + "Clust" + str(clust) + "/cluster.csv", pclustsub)

            runExternalSoft.CrossClusterIC50(pdescsub, pAC50, pclustsub, prbyclust)



        return

    def corelAllAssays(self, cluc, chepg2, chek293):


        prInterfer = self.prout + "interfer/"
        pathFolder.createFolder(prInterfer)

        # cluster main
        fcluster = open(self.pclusters, "r")
        lclusters = fcluster.readlines()

        dclust = {}
        for clusters in lclusters:
            clusters = clusters.strip().split(",")

            chemID =clusters[0]
            clustname = str(clusters[1]) + "_" + str(clusters[2])
            if not clustname in dclust.keys():
                dclust[clustname] = []

            dclust[clustname].append(chemID)

        #open descriptors
        fdesc = open(self.pdesc, "r")
        lchemdesc = fdesc.readlines()
        fdesc.close()

        ddesc = {}
        for chemdesc in lchemdesc[1:]:
            chemID = chemdesc.split("\t")[0]
            ddesc[chemID] = chemdesc


        for cluster in dclust.keys():
            prclustsub = pathFolder.createFolder(prInterfer + str(cluster) + "/")

            #file descriptors
            pdesc = prclustsub + "desc.csv"
            filedesc = open(pdesc, "w")
            filedesc.write(lchemdesc[0])

            for chemical in dclust[cluster]:
                #png
                ppng = self.prPNG + chemical + ".png"
                if path.exists(ppng):
                    copyfile(ppng, prclustsub + chemical + ".png")
                #desc
                if chemical in ddesc.keys():
                    filedesc.write(ddesc[chemical])
            fdesc.close()

            runExternalSoft.crossA50s(pdesc, cluc.pAC50, chepg2.pAC50, chek293.pAC50, prclustsub)
        return






def createSOM(pdesc1D2D, pAC50, corval, maxQuantile, pModel, prSOM):



    # output
    pdesc1D2Dclean = prSOM + "descClean.csv"

    if not path.exists(pdesc1D2Dclean):

        if path.exists(pdesc1D2D) and path.getsize(pdesc1D2D) > 10:
            # preproc
            runExternalSoft.dataManager(pdesc1D2D, 0, corval, maxQuantile, prSOM)
        else:
            print "Error ->", pdesc1D2D


    runExternalSoft.drawEnrichSOM(pdesc1D2Dclean, pAC50, pModel, prSOM)

