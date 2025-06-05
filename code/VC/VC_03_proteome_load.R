# Only ITC and EC were used

options(warn = 0)
library(diann)
library(vsn)
library(sva)

# Loading proteome data ----
validN1U <- read.delim("~/Data/slope-models-adpkd/data/itc_oc_ec_result_N1nr.txt")

sel = grep("Protein.Group|^Norm_",colnames(validN1U), value = T)
ecDat = validN1U[,sel] # 731 1640
colnames(ecDat) = gsub("Norm_", "", colnames(ecDat))

rownames(ecDat) = ecDat$Protein.Group
ecDat$Protein.Group = NULL

ecDat = ecDat[,grep("ec_|itc_", colnames(ecDat))] # 731 1095

# to remove the ones that don't have annotation - dates are not matching
# and also pools are removed
stayclin = which(!is.na(metaclin2$eGFR) & !grepl("re", metaclin2$Batch)) 

metaclin2 = metaclin2[stayclin,] # 847  37

# selecting common samples between proteome and annotation file
commonsampss = intersect(colnames(ecDat), rownames(metaclin2)) # 844
ecDat = ecDat[,commonsampss] # 731 844
metaclin2 = metaclin2[commonsampss,] # 844  37

data.na.ecDat <- apply(ecDat,1,function(x) length(which(is.na(x)))/length(x))
data.s.ecDat <- ecDat[which(data.na.ecDat <= 0.8),] # if it has NA more than 80%, remove it 
# proteins=338 & samples=844
save(data.s.ecDat, file = paste0("~/Data/slope-models-adpkd/data/data_s_ecDat.RData"))
