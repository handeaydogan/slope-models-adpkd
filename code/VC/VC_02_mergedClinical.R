# MERGING CLINICAL AND PROTEOME ANNOTATION DATA OF VALIDATION COHORTS ----
# | # Preparing the Main Matching Data ----
# getting proteome sample annotation 
protann <- read.delim("~/Data/slope-models-adpkd/data/itc_oc_ec_AN_IDs.txt") # 1644   18
protann$sample_name[grep(" ", protann$sample_name)] = gsub(" ", "", protann$sample_name[grep(" ", protann$sample_name)])
protann$RunDate = as.Date(gsub("Mai","May", protann$RunDate), 
                          format = "%y_%b_%d")
# changing the column name
protann$matchingclin = protann$X
protann$X = NULL

# matching ID and sample name is changed to "pool"
protann$matchingclin[protann$Type == "pool"] = "pool"
protann$sample_name[protann$Type == "pool"] = "pool"
rownames(protann) = protann$Name

# meta matching file is created to match the clinical 
# annotation with proteome annotation
matchinDF = protann # 1644   18

# extracting indexes of external cohort, the ones that don't have 
# patient ID (only have serum ID) and healthy cohort  
duplictn = c(grep("\\D", matchinDF$sample_name, invert = T), 
             grep("AB-", matchinDF$sample_name))

# transferring the sample names of duplictn into matchingclin column
matchinDF[duplictn,"matchingclin"] = matchinDF[duplictn,"sample_name"] # 1644   18

# getting the proteome annotation of internal-temporal cohort (ITC, project itc)
library(readxl)
cannotitc <- read_excel("~/Data/slope-models-adpkd/data/annotationitc.xlsx")
cannotitc = cannotitc[,c(2:5)] # 849   4
cannotitc$sampleDate = as.Date(cannotitc$sampleDate)
cannotitc$matchingclin = cannotitc$Iknum
colnames(cannotitc) = c("Iknum", "identifieritc", "batch", 
                         "sampleDateitc", "matchingclin") # 849   5

# some patients only had serum ID, preparing ID conversion df
library(readxl)
Serumnumbers_oc <- read_excel("~/Data/slope-models-adpkd/data/Serumnumbers_oc.xlsx") # 169   4
Serumnumbers_oc$`ID from TA` = NULL
colnames(Serumnumbers_oc) = c("matchingclin","sampleDateoc","identifieroc")
Serumnumbers_oc$matchingclin = as.character(Serumnumbers_oc$matchingclin) # 169   3

library(dplyr)
# merging ITC clinical annotation and serum IDs with matchinDF
matchinDF2 = left_join(matchinDF,cannotitc, by = "matchingclin") # 1644 22
matchinDF3 = left_join(matchinDF2,Serumnumbers_oc, by = "matchingclin") # 1644 24

# merging sampling dates
matchinDF4 = matchinDF3 %>% 
  mutate(merged_sampleDates = coalesce(sampleDateoc,sampleDateitc)) # 1644 25

# replacing serum IDs with identifier IDs
matchinDF4$identifiers = matchinDF4$sample_name
chnging = setdiff(grep("\\D",matchinDF4$identifiers,invert = T),
                  grep("xxxx",matchinDF4$identifiers)) # IDs are masked
matchinDF4$identifiers[chnging] = matchinDF4$identifieroc[chnging] # 1644   26

# fixing naming of one sample which has the seperator ":" 
# were used instead of "-"
rownames(matchinDF4) = matchinDF4$Name # 1644   26

# | # External Cohort Annotation ----
source("~/Data/slope-models-adpkd/code/VC/VC_01_ECclinical.R", echo = TRUE)
# only baseline measurements were in the proteome
# subsetting the clinical data for baseline 
cannot.ext.last = cannot.long[grep("_m01",cannot.long$tpID),] # 228  37

# column names are changed to make them compatible with naming in 
# screening cohort (SC) 
cannot.ext.last$before = cannot.ext.last$slope
cannot.ext.last$merged_eGFR = cannot.ext.last$patient__eGFR
cannot.ext.last$selectedmut = cannot.ext.last$first_mutation_category

cannot.ext3 = cannot.ext.last # 228  39
cannot.ext3$sample_name = cannot.ext3$protmatch

cannot.ext3$mayo5 = cannot.ext3$mayo

cannot.ext3$arterial_hypertension = NA
cannot.ext3$arterial_hypertension[which(cannot.ext3$hypertension_binary == "yes")] = 3
cannot.ext3$arterial_hypertension[which(cannot.ext3$hypertension_binary == "no")] = 0

# extracting relevant information
clinec = cannot.ext3[,c("sex","age","patient__eGFR","mayoClass","mayo5",
                        "selectedmut", "before","sample_name", 
                        "date_lab_serum", "arterial_hypertension", 
                        "count", "bTF", "tkv")]
clinec$fampos = NA
clinec$usymp = NA # 228  #15

# | #  Internal/Temporal Cohort Annotation -----
source("~/Data/slope-models-adpkd/code/00functions.R", echo = TRUE)
source("~/Data/slope-models-adpkd/code/SC/SC_01_annot.half.R", echo = TRUE)
source("~/Data/slope-models-adpkd/code/SC/SC_02_interventions.R", echo = TRUE)
source("~/Data/slope-models-adpkd/code/SC/SC_03_mergedCreasUpd.R", echo = TRUE)

# generating matching ID for proteome annotation of ITC to match it with 
# annot.half (the main ADPKD registry clinical annotation file)
cannotitc$matching_itc = paste0(cannotitc$identifieritc, "&", 
                                  cannotitc$sampleDateitc) # 849 6
examination_itc = cannotitc[!is.na(cannotitc$sampleDateitc),] # 831 6

annot.half$matching_itc = paste0(annot.half$patient_id, "&",
                                  annot.half$exam_date)
sematch = intersect(annot.half$crea.match.se, cannotitc$matching_itc) # 630
exammatch = intersect(annot.half$matching_itc, cannotitc$matching_itc) # 630
setdiff(exammatch,sematch)

# changing some dates to make them compatible 
annot.half[grep(paste(setdiff(sematch, exammatch), collapse = "|"), 
                annot.half$crea.match.se),"matching_itc"] = setdiff(sematch, 
                                                                     exammatch)
# merging it with annot.half
mergeditc = left_join(examination_itc, annot.half, 
                       by = "matching_itc") # 831 333

library(plyr)
# merging it with proteome annotation and preparing the data
mergeditc = rbind.fill(mergeditc, 
                        cannotitc[is.na(cannotitc$sampleDateitc),]) # 849 333
mergeditc$batchID = mergeditc$batch
mergeditc$batch = as.factor(gsub("_.*","",mergeditc$batchID))
rownames(mergeditc) = mergeditc$batchID
mergeditc$patient__CKD = as.factor(eGFRtoStageAB(mergeditc$patient__eGFR))
dim(mergeditc) # 849 334

library(stringr)
rownames(mergeditc) = str_sub(mergeditc$Iknum, 3, nchar(mergeditc$Iknum[1]))
mergeditc$Zuord = as.factor(str_sub(mergeditc$Iknum, 1, 2))
mergeditc$batch = as.factor(gsub("pool", "batch", mergeditc$batch))

# adjusting th levels of mayo class
mergeditc$mayoClass = mergeditc$adpkd_progression__mayo
mergeditc$mayoClass = as.factor(mergeditc$mayoClass)
levels(mergeditc$mayoClass) = c("1A" = "Less", "1B" = "Less", "1C" = "Mid",
                                 "1D" = "More", "1E" = "More")
mergeditc$mayo5 = mergeditc$adpkd_progression__mayo # 849 337

# adjusting the levels of PKD1 and PKD2 genotype
mergeditc$Class_ADPKD[which(is.na(mergeditc$Class_ADPKD))] = "none/others"
mergeditc$selectedmut = as.factor(mergeditc$Class_ADPKD)
levels(mergeditc$selectedmut) = list(unknown = "none/others", 
                                      PKD1vus = "PKD1 VUS",
                                      PKD2 = "PKD2", PKD2vus = "PKD2 VUS",
                                      PKD1NT = "PKD1 non-truncating",
                                      PKD1T = "PKD1 truncating")
# extracting relevant information
clinitc = mergeditc[,c("patient__gender", "patient__age", "patient__eGFR",
                         "mayoClass", "mayo5", "selectedmut", "before", "Iknum", 
                         "bTF", "arterial_hypertension", 
                         "tomography__kidney_volume",
                         "family__family_members_with_adpkd",
                         "urological_symptoms")] # 849  13

#Merging Clinical Data ----
metaclin = left_join(matchinDF4, clinec, by = "sample_name") # 1644 40 
metaclin = left_join(metaclin, clinitc, by = "Iknum") # 1644 52

# coalescing the common variables and renaming them
metaclin2 = metaclin %>% mutate(gender = coalesce(patient__gender,sex))
metaclin2 = metaclin2 %>% mutate(age = coalesce(patient__age,age))
metaclin2 = metaclin2 %>% mutate(mayo = coalesce(mayoClass.x,mayoClass.y))
metaclin2 = metaclin2 %>% mutate(mayo5 = coalesce(mayo5.x,mayo5.y))
metaclin2 = metaclin2 %>% 
  mutate(eGFR = coalesce(patient__eGFR.x,patient__eGFR.y))
metaclin2 = metaclin2 %>% 
  mutate(selectedmut = coalesce(selectedmut.x,selectedmut.y))
metaclin2 = metaclin2 %>% mutate(before = coalesce(before.x,before.y))
metaclin2 = metaclin2 %>% mutate(bTF = coalesce(bTF.x,bTF.y))
metaclin2 = metaclin2 %>% 
  mutate(hypertension = coalesce(arterial_hypertension.x,arterial_hypertension.y))
metaclin2 = metaclin2 %>% 
  mutate(postfam = coalesce(fampos,family__family_members_with_adpkd))
metaclin2 = metaclin2 %>% mutate(usyeni = coalesce(usymp,urological_symptoms))
metaclin2 = metaclin2 %>% 
  mutate(tkvmerged = coalesce(tkv,tomography__kidney_volume))
metaclin2 = metaclin2 %>% 
  mutate(merged_dates = coalesce(merged_sampleDates,date_lab_serum))
# 1644   64

# removing previous colmuns
metaclin2[,c("patient__gender", "patient__age", "patient__eGFR.x", "mayoClass.x", 
             "fampos", "selectedmut.x", "before.x", "Iknum", "bTF.x", "bTF.y",
             "sex","patient__eGFR.y","mayoClass.y", "selectedmut.y", "before.y", 
             "usymp", "sample_name", "mayo5.x", "mayo5.y", 
             "arterial_hypertension.x", "urological_symptoms", 
             "arterial_hypertension.y", "family__family_members_with_adpkd",
             "tkv", "tomography__kidney_volume", "merged_sampleDates", 
             "date_lab_serum")] = NULL # 1644   37

# changing the rownames
rownames(metaclin2) = metaclin2$Name #1644   37

# the ones which have duplicates and mazo class 2 will be removed
remwdos = grep(paste(unique(c(mayoclass2,dupsFound,
                              substring(gsub("\\.|-","",mayoclass2G),2))),
                     collapse = "|"),metaclin2$identifiers)

metaclin2 = metaclin2[-remwdos,]  # 1605   37
