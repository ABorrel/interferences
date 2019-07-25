#!/usr/bin/env Rscript
source("cardMatrix.R")

################
#     MAIN     #
################

args <- commandArgs(TRUE)
pfilin = args[1]

#pfilin = "/home/borrela2/interference/spDataAnalysis/predictions/test/pred"

din = read.csv(pfilin, sep = "\t", header = TRUE)
rownames(din) = din[,1]
din = din[,-1]
din = - din
calibrate = seq(0,15, by = (15/(dim(din)[1]-1))) 
din = cbind(din, calibrate)

cardMatrix(t(din), pfilin, "red")

