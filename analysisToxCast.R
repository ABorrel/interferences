library(ggplot2)





dToxCastAC50 = read.csv("./../data/ac50_Matrix_170614.csv", sep = ",")

dToxCastChemicals = read.csv("./../data/toxcast_chemicals_2017-06-14.csv")


# analyse master table assays #
###############################
dToxCastAssay = read.csv("./../data/ToxCast_assay_master_2017-06-14.csv", sep = ",")


print("assay_source_name")
lassaySource = dToxCastAssay[,which(colnames(dToxCastAssay)=="assay_source_name")]
print(unique(lassaySource))
print(length(unique(lassaySource)))

print("assay_design_type_sub")
lassaydesignsub = dToxCastAssay[,which(colnames(dToxCastAssay)=="assay_design_type_sub")]
print(unique(lassaydesignsub))
print(length(unique(lassaydesignsub)))

print("detection_technology_type")
lassaytech = dToxCastAssay[,which(colnames(dToxCastAssay)=="detection_technology_type")]
print(unique(lassaytech))
print(length(unique(lassaytech)))

print("detection_technology_type_sub")
lassaytechsub = dToxCastAssay[,which(colnames(dToxCastAssay)=="detection_technology_type_sub")]
print(unique(lassaytechsub))
print(length(unique(lassaytechsub)))

print("detection_technology")
lassaydetectechno = dToxCastAssay[,which(colnames(dToxCastAssay)=="detection_technology")]
print(unique(lassaydetectechno))
print(length(unique(lassaydetectechno)))



