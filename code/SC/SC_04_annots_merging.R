# preparing the proteome annotation (cannot) and merging it 
# with the main annotation file (annot.half) 

library(openxlsx)
setwd("~/Data/slope-models-adpkd")
library(openxlsx)
cannot <- read.xlsx("~/data/ADPKDComplete_v2.xlsx")

# removing some variables
cannot[,c("Monat", "CKD.Stage.G", "Control.vs..Affected", "eGFR", "Alter..in.Jahren.", "Jahr")] = NULL

# patients who don't have identifier ID
identifieridchange = which(cannot$identifier == cannot$Labornummer.Serum)
library(readxl)
Serumnumbers_oc <- read_excel("data/Serumnumbers_oc.xlsx")
cannot[identifieridchange,"identifier"] = Serumnumbers_oc$identifier[match(cannot[identifieridchange,"Labornummer.Serum"],
                                                     Serumnumbers_oc$`ID as in Plate`)]

# preparing the annotation data
cannot$Zuordnung <- gsub("Kontrolle","AB",cannot$Zuordnung) #added to change the names kontrolle to AB
cannot$Labornummer.Serum = as.character(cannot$Labornummer.Serum)
cannot$ID.sc <- gsub("xxx","xxx",cannot$ID.sc) # IDs are masked
cannot$Geschlecht = gsub('\\s+', '',cannot$Geschlecht)
cannot$eGFR_v2 = as.numeric(cannot$eGFR_v2)

# converting dates
cannot$Einschluss[grep("10", nchar(cannot$Einschluss))] = as.character(as.Date(cannot$Einschluss[grep("10", nchar(cannot$Einschluss))],format = '%d/%m/%Y'))
cannot$Einschluss[grep("5", nchar(cannot$Einschluss))] = as.character(as.Date(as.numeric(cannot$Einschluss[grep("5", nchar(cannot$Einschluss))]), origin = "1899-12-30"))

cannot$Geburtsdatum = as.Date(cannot$Geburtsdatum, format = '%d/%m/%Y')
cannot$Datum.Untersuchung = as.Date(cannot$Datum.Untersuchung, format = '%d/%m/%Y')
cannot$Einschluss = as.Date(cannot$Einschluss)

cannot$CKD = eGFRtoStageAB(cannot$eGFR_v2)
cannot$CKD_v2 = NULL

cannot$Mayo.class = gsub("^-$", NA, cannot$Mayo.class)
cannot = cannot[-grep("Typ 2",cannot$Mayo.class),] # removing mayo class 2

wbfactor = names(which(unlist(lapply(cannot, function(x) length(unique(na.omit(x))))) < 15))
cannot[,wbfactor] = lapply(cannot[,wbfactor], factor)

# pulling the dates from main annotation data since every sample 
# in the proteome were baseline
pullProtSampDate = annot.half %>% group_by(patient_id) %>% filter(examination == "baseline")
pullProtSampDate = pullProtSampDate[,c("patient_id","examination", "exam_date")]
cannot$protSamp = pullProtSampDate$exam_date[match(cannot$identifier, pullProtSampDate$patient_id)] # 495 20

cannot$matching_sc = paste0(cannot$identifier, "&", cannot$protSamp)
utershuch = cannot[!is.na(cannot$protSamp),]

merged.cannot = left_join(utershuch, annot.half, by = "matching_sc") # 329

library(plyr) 
merged.cannot = rbind.fill(merged.cannot, cannot[is.na(cannot$protSamp),])# 495 348
merged.cannot$batch = as.factor(gsub("_.*","",merged.cannot$ID.sc))

# removing the ones that don't have batch number
merged.cannot = merged.cannot[-which(is.na(merged.cannot$ID.sc)),] # 357 349

# assigning rownames
rownames(merged.cannot) = merged.cannot$ID.sc

merged.cannot = merged.cannot %>% mutate(merged_eGFR = coalesce(patient__eGFR,eGFR_v2))
merged.cannot = merged.cannot %>% mutate(merged_mayo = coalesce(adpkd_progression__mayo,Mayo.class))
merged.cannot = merged.cannot %>% mutate(merged_age = coalesce(patient__age,age))

merged.cannot$patient__CKD = as.factor(merged.cannot$patient__CKD)
merged.cannot = merged.cannot %>% mutate(merged_CKD = coalesce(patient__CKD,CKD))

merged.cannot$Geschlecht = gsub("m", "male", merged.cannot$Geschlecht)
merged.cannot$Geschlecht = gsub("w", "female", merged.cannot$Geschlecht)
merged.cannot = merged.cannot %>% mutate(merged_gender = coalesce(patient__gender,Geschlecht))
merged.cannot$merged_gender = as.factor(merged.cannot$merged_gender) #357 215

merged.cannot[,c("Mayo.class", "adpkd_progression__mayo", "patient__age", "age", 
                 "patient__gender", "Geschlecht", "eGFR_v2", "adpkd_progression__egfr", 
                 "patient__eGFR", "CKD", "patient__CKD", "adpkd_progression__ckd")] = NULL # 357 342

rem.cols.v = unique(c(which(lapply(merged.cannot, function(x){length(na.omit(x))}) < 2), # the ones that have only one value in the column # 29
                    which(lapply(merged.cannot, function(x){length(unique(na.omit(x)))}) < 2))) # the ones that have only one *unique* value in the column
rem.cols.all = c(grep("medical_examination__age|urine_spontaneous__tolvaptan_therapy",
                    colnames(merged.cannot)),rem.cols.v) # 33
merged.cannot = merged.cannot[,-rem.cols.all] # 357 309
