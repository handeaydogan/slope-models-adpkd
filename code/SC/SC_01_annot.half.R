library(plyr)
library(dplyr)
library(purrr)
# Opening all files ----
setwd("~/Data/ADPKD_Cohort_ClinicalD") 
temp = list.files()
myfiles = lapply(temp, read.csv) # reading every csv file in the folder
myfiles1 = lapply(myfiles, function(x){ # empty strings to NAs
  replace(x, x == "", NA)
})

# include file name in column names
counter = 1
forcolnam = gsub("^_|__|_$","",gsub("[0-9]|.csv","",temp))
myfiles2 = lapply(myfiles1, function(x){
  colnames(x) <- paste0(forcolnam[counter],"__",colnames(x))
  counter <<- counter + 1
  return(x)
})

# Removing Keys ----
# selecting files contain two columns --> annotation files 
# (for some columns in other files)
levelsH = which(unlist(lapply(myfiles2, ncol)) == 2) 
myfiles3 = myfiles2[-levelsH]

#Fixing column names ----
patientIDs = which(unlist(lapply(myfiles3, function(x) any(grep("patient_id", 
                                                                colnames(x))))))
myfiles3[patientIDs] = lapply(myfiles3[patientIDs], function(x){
  colnames(x)[grep("patient_id", colnames(x))] = grep("__patient_id", 
                                                      gsub("__x_patient_id","__patient_id",
                                                           colnames(x)),value = T)
  return(x)
})

# Selecting files with "examination" ----
# selecting data frames containing columns named as their ID and examination
items = which(unlist(lapply(myfiles3, function(x) any(grep("examination", colnames(x))))))

# change some more column names before starting
myfiles3[items] = lapply(myfiles3[items], function(x){
  colnames(x)[grep("x_examination|__examination", colnames(x))] = 
    grep("__examination", gsub("x_examination","examination",colnames(x)),
         value = T)
  return(x)
})

library(reshape2)
# generating unique ids with their ID and examination (baseline, follow_up_1)
# and merging them (dataframes)
withExam = myfiles3[items]
withExam.list = lapply(withExam, function(x){
  x$matching <- paste0(x[,grep("patient_id", colnames(x))],
                       "&", x[,grep("x_examination|__examination", colnames(x))])
  return(x)
})
withExam.merged.all = Reduce(function(x, y) merge(x, y, by = "matching", all = T), withExam.list)

# adding patient id and examination
withExam.merged.all$patient_id = gsub("&.*", "", withExam.merged.all$matching)
withExam.merged.all$examination = gsub(".*&", "", withExam.merged.all$matching) # 3439  209

# Selecting files with "patient ID" ----
# the ones that don't have examination column (baseline or followup)
withoutexam = myfiles3[-items]
withoutexam = lapply(withoutexam, function(x){
  colnames(x)[grep("patient_id",colnames(x))] <- "patient_id"
  colnames(x)[grep("parent",colnames(x))] <- "parent"
  return(x)
})

haveID = which(lapply(withoutexam, function(x) any(grep("patient_id", 
                                                        colnames(x)))) == T)
haveparent = which(lapply(withoutexam, function(x) any(grep("parent", 
                                                            colnames(x)))) == T)
onlyParent = setdiff(haveparent,haveID)
onlyID = setdiff(haveID, haveparent)
bothPID = intersect(haveID, haveparent) 

# leaving family info out from merging
parentpstIDs.list = withoutexam[c(onlyID, bothPID)] 

# Merging examination and patient ID files (prev two sections) ----

# patient, fimly, healt_status, meta, mutation csv files were merged.
merged2 = c(list(withExam.merged.all),parentpstIDs.list[1:5]) 
mergedExtolvs = merged2 %>% reduce(left_join, by = "patient_id") 
#combined #3439  276 #1314 unique patients

# Fixing data related things ----
patterns.checking = names(which(sort(table(unlist(lapply(colnames(mergedExtolvs), function(x) gsub(".*__","", x))))) > 3))
patterns.checking = paste0("__", unique(gsub("^x_", "", patterns.checking)))

# 2:date, 4:id, 5:exam (4 and 5 all true)
cbind.data.frame(patterns.checking, do.call(rbind,lapply(patterns.checking, function(y){
  table(apply(mergedExtolvs[,grep(y, colnames(mergedExtolvs))], 1, 
              function(x) length(unique(na.omit(x))) == 1))
})))# all examination and ids are the same, so you can remove it

# removing the same columns
wbremoved = grep("__patient_id|__examination", colnames(mergedExtolvs))
extra.cols.examID = mergedExtolvs[,wbremoved]
mergedExtolvs = mergedExtolvs[,-wbremoved]

extra.cols.parent.ident.imperror = mergedExtolvs[,grep("parent|__ident|__import_errors", colnames(mergedExtolvs))]
annot.half = mergedExtolvs[,-grep("parent|__ident|__import_errors", colnames(mergedExtolvs))] # 3439  200

# changing dates
annot.half[,grep("date", colnames(annot.half))] = lapply(annot.half[,grep("date", colnames(annot.half))],
                                                       function(x){
                                                         if (any(na.omit(nchar(x)) == 10)) {
                                                           x = as.Date(x, format = "%Y-%m-%d")
                                                         }
                                                         else if (any(na.omit(nchar(x)) > 10)) {
                                                           x = as.Date(x, format = "%Y-%m-%d %H:%M:%OS")
                                                         }
                                                       })

#removing some columns before starting
rem.cols.v = unique(c(which(lapply(annot.half, function(x){length(na.omit(x))}) < 2), 
                      # the ones that have only one value in the column # 0
                    which(lapply(annot.half, function(x){length(unique(na.omit(x)))}) < 2)))
                      # the ones that have only one *unique* value in the column # 3

annot.half = annot.half[,-rem.cols.v]
dim(annot.half) # 3439  197

# assign classes
# # factors
factor.colnams = colnames(annot.half)[intersect(which(unlist(lapply(annot.half, function(x) length(unique(na.omit(x))))) < 15),
                                              which(unlist(lapply(annot.half, function(x) length(na.omit(x)))) > 15))]
factor.colnams2 = grep("no_members|__number_of", factor.colnams, value = T, invert = T)
annot.half[,factor.colnams2] = lapply(annot.half[,factor.colnams2], factor)

# # integers to numeric
annot.half[,which(unlist(lapply(annot.half, class)) == "integer")] = lapply(annot.half[,which(unlist(lapply(annot.half, class)) == "integer")], as.numeric) 
# 3439  197

# Mayo class 2 ----
library(readxl)
#additional information regarding mayo class, dialysis and Tolvaptan usage
mct_table <- read_excel("~/Data/slope-models-adpkd/data/Mayo2_Stand231123.xlsx")
mct_table = as.data.frame(mct_table)
rownames(mct_table) = mct_table$identifier
mct_table$`Mayo 2` = as.logical(mct_table$`Mayo 2`)
mayoclass2 = mct_table[which(mct_table$`Mayo 2` == TRUE),"identifier"]

# Duplication - Dates ----
dupsFound = "xxxxxxxxxxxxxxxx" # IDs are masked # duplicated

wbdel = unique(c(mayoclass2,dupsFound))
annot.half = annot.half[-grep(paste(wbdel, collapse = "|"), annot.half$patient_id),] #3318  197 

# Removing duplicated serum dates and patient followups with no dates----
datesFindNA = as.data.frame(annot.half[,c("patient_id","examination","matching",
                                          "haematology__date", "medical_examination__date",
                                          "serum__date", "tomography__date",
                                          "urine_collection__date", "urine_sediment__date",
                                          "urine_spontaneous__date", "urine_stix__date")])
rownames(datesFindNA) = datesFindNA$matching
datesFindNAlist = split(datesFindNA , f = datesFindNA$patient_id)
datesFindNAlist2 = lapply(datesFindNAlist, function(x) apply(x, 1, function(y) length(which(is.na(y)))))
# no date is available

alldateNA = as.character(unlist(lapply(datesFindNAlist2, function(x) names(which(x > 7)))))

# only one date available
onedateav = as.character(unlist(lapply(datesFindNAlist2, function(x) names(which(x == 7)))))
removalDupslist = list()
invisible(lapply(onedateav, function(names){
  identifierid = gsub("&.*", "", names)
  indxdatenotna = which(!is.na(datesFindNAlist[[identifierid]][names,]))
  namedatenotna = colnames(datesFindNAlist[[identifierid]][names,])[indxdatenotna]
  namedatenotna = setdiff(namedatenotna, c("patient_id", "examination","matching"))

  allserum = datesFindNAlist[[identifierid]][,"serum__date"]
  isdups = datesFindNAlist[[identifierid]][names,"serum__date"]
  dups = which(allserum == isdups)
  if (length(dups) > 1) {
    removalDupslist[[paste(names)]] <<- names
  }
}))
serumdups = as.character(unlist(removalDupslist)) # 378 sample were removed
dupsnnas = c(serumdups,alldateNA)

annot.dupsnnas = annot.half[grep(paste(dupsnnas, collapse = "|"), annot.half$matching),]
annot.half = annot.half[-grep(paste(dupsnnas, collapse = "|"), annot.half$matching),] 
# 2940  197  # 1273 unique patients

# Fixing height -----
annot.half$vitality__height_eski = annot.half$vitality__height
wbnasheight = unique(c(which(annot.half$vitality__height < 105), 
                       which(annot.half$vitality__height > 300))) # wb NAs
annot.half$vitality__height[wbnasheight] = NA

wbmedian = grep("xxxxxxxxxxxxxxxx", # IDs are masked
                annot.half$patient_id)
annot.half[wbmedian,] = annot.half[wbmedian,] %>% 
  group_by(patient_id) %>% 
  mutate(vitality__height = median(vitality__height,na.rm = T))

annot.half = annot.half %>% group_by(patient_id) %>% 
  dplyr::mutate(avheight = mean(vitality__height, na.rm = T)) #tibble

# TKV to htTKV ----
library(readr)
missingmayos <- read_csv("~/Data/slope-models-adpkd/data/missingmayos.csv")
tkvtransf = annot.half[grep("identifier", annot.half$patient_id),
                       c("patient_id","tomography__date", 
                         "medical_examination__date", 
                         "tomography__kidney_volume", 
                         "adpkd_progression__httkv",
                         "tomography__comment", "avheight",
                         "tomography__method")]

# selecting the ones that don't have TKV from the missingmayos variable
indxxxes = grep(paste(grep("identifier",missingmayos$Var1,value = T), collapse = "|"), 
                tkvtransf$patient_id)
tkvtransf = tkvtransf[indxxxes,]

# which one has a comment 
cmmt = which(!is.na(tkvtransf$tomography__comment))
splttedtxt = unlist(strsplit(tkvtransf$tomography__comment[cmmt], split = " "))

# changing TKV for that patient
tkvtransf$tomography__kidney_volume[cmmt] = sum(as.numeric(grep("^?[0-9.]+$",
                                                                splttedtxt, 
                                                                as.numeric, 
                                                                value = T)))
# (327+400)#(pi/6)*(323+400) # because method is is 1 (MRT)

# removing the ones that don't have either tomography date or kidney volume
tkvtransf = tkvtransf[-which(is.na(tkvtransf$tomography__date) | 
                               is.na(tkvtransf$tomography__kidney_volume)),]

# creating a matching key
tkvtransf$matchtkv = paste0(tkvtransf$patient_id, "&", tkvtransf$tomography__date)

# heights of all patients whose samples were used
library(readr)
allpatientsIDs_usedinAnalysis <- read_csv("~/Data/slope-models-adpkd/data/allpatientsIDs_usedinAnalysis.csv")
allpat = grep("identifier", allpatientsIDs_usedinAnalysis$x, value = T)

annot.half$tomography__kidney_volume_old = annot.half$tomography__kidney_volume
annot.half$tomography__kidney_volume = NULL
annot.half$adpkd_progression__httkv_old = annot.half$adpkd_progression__httkv
annot.half$matchtkv = paste0(annot.half$patient_id, "&", annot.half$tomography__date)

annot.half.upd = annot.half
annot.half.upd = left_join(annot.half.upd,
                           tkvtransf[,c("tomography__kidney_volume","matchtkv")], 
                           by = "matchtkv")

annot.half.upd$tomography__kidney_volume = coalesce(annot.half$tomography__kidney_volume_old,
                                                    annot.half.upd$tomography__kidney_volume)


annot.half.upd$httkvcalc = annot.half.upd$tomography__kidney_volume/(annot.half.upd$avheight/100)
annot.half.upd$adpkd_progression__httkv_old = annot.half.upd$adpkd_progression__httkv
annot.half.upd$adpkd_progression__httkv = NULL
annot.half.upd$adpkd_progression__httkv = coalesce(annot.half$adpkd_progression__httkv_old,
                                                   annot.half.upd$httkvcalc)

annot.half = annot.half.upd # 2940  203

# #n patient=1106
# paste0("npatient=", length(unique(grep("identifier", annot.half$patient_id, value = T))))
# #n observation=2630
# paste0("nobservation=", length((grep("identifier", annot.half$patient_id, value = T))))
# #n followup range
# paste0("nfollowupRange=[", max(table(grep("identifier", annot.half$patient_id, value = T))),"-", # 9
#        min(table(grep("identifier", annot.half$patient_id, value = T))), "]") # 1

