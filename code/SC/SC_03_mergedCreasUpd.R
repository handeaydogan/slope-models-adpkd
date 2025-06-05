
# uploading the externally measured serum creatinine values  
hist.crea <- read.csv("~/Data/ADPKD_Cohort_ClinicalD/_Creatinine_NocoDB_csv__202309191451.csv") 
hist.crea$date = as.Date(hist.crea$date ,format = "%d.%m.%y")
hist.crea$crea.match = paste0(hist.crea$id1, "&", hist.crea$date)
hist.crea$created_at = NULL
hist.crea$updated_at = NULL # 6225 6 # 727 unique patient

# removing wbdel (generated in SC_01_annothalf.R)
# # patients with duplicated dates and mayo class 2
hist.crea = hist.crea[-grep(paste(wbdel, collapse = "|"), hist.crea$id1),] # 6053 6 # 704 unique patient

# removing duplicated IDs & dates in hist.crea
hist.crea = hist.crea[-which(duplicated(hist.crea$crea.match)),] # 5958 6 # 704 unique patients

# identification mat for patients
identf.mat = unique(annot.half[,c("patient_id","patient__birthdate", "patient__gender")]) # 1273  3

library(dplyr)
hist.crea$patient_id = hist.crea$id1
hist.crea2 = left_join(hist.crea,identf.mat, by = "patient_id") # 5958   9

# creating matching for annot.half
annot.half$crea.match.med = ifelse(is.na(annot.half$medical_examination__date), NA,
                                 paste0(annot.half$patient_id, "&", # length(na.omit(annot.half$crea.match.med)) # 2932
                                        annot.half$medical_examination__date))
annot.half$crea.match.se = ifelse(is.na(annot.half$serum__date), NA,
                                paste0(annot.half$patient_id, "&", # length(na.omit(annot.half$crea.match.med)) # 2809
                                       annot.half$serum__date))
annot.half$matchingHistcrea = annot.half$crea.match.med # because this one has less mistakes
annot.half_mut = split(annot.half , f = annot.half$patient_id)

# trying to find out which has more duplication --> it was serum_date
dupsMedSe = list(); dupsSeMed = list()
invisible(lapply(annot.half_mut, function(x){
  dupsdupsmedse = duplicated(x$medical_examination__date) | duplicated(x$medical_examination__date, fromLast = TRUE)
  if (any(dupsdupsmedse)) {
    dupsMedSe[paste0(x$patient_id)] <<- "DUPS"
  }
  dupsdupssemed = duplicated(x$serum__date) | duplicated(x$serum__date, fromLast = TRUE)
  if (any(dupsdupssemed)) {
    dupsSeMed[paste0(x$patient_id)] <<- "DUPS"
  }
}))
reqchange = sort(unique(c(names(unlist(dupsMedSe)), names(unlist(dupsSeMed)))))

# changing some dates to make them compatible 
findreqchange = do.call(rbind.data.frame, annot.half_mut[reqchange])[,c("matching","serum__date","medical_examination__date","matchingHistcrea","crea.match.se")]
findreqchange$datedif = abs(findreqchange$serum__date - findreqchange$medical_examination__date)
changedat = findreqchange[which(findreqchange$datedif < 180 & findreqchange$datedif > 0 ),]
wbchanging = changedat[which(!is.na(match(changedat$crea.match.se, hist.crea2$crea.match))),"matching"]
annot.half[annot.half$matching == as.character(wbchanging),"matchingHistcrea"] = as.character(changedat[which(!is.na(match(changedat$crea.match.se, hist.crea2$crea.match))),"crea.match.se"])
# 2940  206

# detailed identification mat from annot.half
identf.mat = annot.half[,c("patient_id","patient__birthdate", "patient__gender", 
                           "serum__creatinine", "matching", 
                           "medical_examination__date", "matchingHistcrea",
                           "crea.match.med", "crea.match.se", 
                           "adpkd_progression__mayo","adpkd_progression__httkv",
                           "adpkd_progression__httkv_old",
                           grep("^medication__",colnames(annot.half),value = T))] # 2940   23

# generating matching for hist.crea
hist.crea2$matchingHistcrea = paste0(hist.crea2$id1, "&", as.Date(hist.crea2$date ,format = "%d.%m.%y"))
merged.creas = merge(hist.crea2, identf.mat, by = "matchingHistcrea", all = T) # 8124   32

merged.creas2 = merged.creas %>% mutate(merged_creat = coalesce(serum__creatinine,crea))
merged.creas2 = merged.creas2 %>% mutate(exam_date = coalesce(date,medical_examination__date))
merged.creas2[,c("crea", "serum__creatinine", "date", "medical_examination__date")] = NULL # 8124   30

# merging several columns (because its duplicated)
merging.colnames = names(table(gsub("\\.x|\\.y","", grep("\\.x|\\.y", colnames(merged.creas2), value = T))))
merged.creas3 = merged.creas2

for (i in merging.colnames) {
  ind.i = grep(i, colnames(merged.creas3))
  merged.creas3[,ncol(merged.creas3) + 1] <- coalesce(merged.creas3[,ind.i[1]],merged.creas3[,ind.i[2]])
  colnames(merged.creas3)[ncol(merged.creas3)] <- i  # Rename columns
  merged.creas3 = merged.creas3[,-ind.i]
} # 8124   27

# calculate age variable from date of birth
merged.creas3$patient__age = as.numeric(round(difftime(merged.creas3$exam_date,
                                                       merged.creas3$patient__birthdate,
                                                       units = "days")/365.25,2))
# calculate eGFR
merged.creas3$patient__eGFR = computeGFR_upd(merged.creas3,
                                             dob = "patient__age",
                                             cdate = NULL,
                                             creatinine = "merged_creat",
                                             gender = "patient__gender")
# calculate CKD stage
merged.creas3$patient__CKD = eGFRtoStageAB(merged.creas3$patient__eGFR) #8124 30

# merging merged.creas3 with data containing interventions and their dates
mysamps.annot = left_join(merged.creas3, tolv_used_p2, by = "patient_id") #8124  48 
mysamps.annot$adpkd_progression__mayo_org = mysamps.annot$adpkd_progression__mayo #original mayo groups

# split the merged.crea3 according to patient IDs
mysamps.annot <- split(mysamps.annot, f = merged.creas3$patient_id)

# adding dates of interventions (Tolvaptan,transplant,nephrectomy and dialysis)
# to create "before" slopes
invisible(lapply(tolv_used_p2$identifier, function(dnamx){
  if (!is.na(match(dnamx,names(mysamps.annot)))) {
    diffdates = mysamps.annot[[dnamx]]$`Start 45/15` - mysamps.annot[[dnamx]]$exam_date
    befores = which(diffdates > 0)
    
    mysamps.annot[[dnamx]]$tolv_upd <<- 2
    if (length(befores) > 0) {
      mysamps.annot[[dnamx]]$tolv_upd[befores] <<- 1
    }
    mysamps.annot[[dnamx]]$tolv.ever_upd <<- ifelse(any(mysamps.annot[[dnamx]]$tolv_upd == 2),"yes", "no")
  }
}))

# set patient's variables to 1, no and NA
# patients do not have any info on intervention were assumed to not used intervention
invisible(lapply(setdiff(names(mysamps.annot), tolv_used_p2$identifier), function(x){
  mysamps.annot[[x]]$tolv_upd <<- 1
  mysamps.annot[[x]]$tolv.ever_upd <<- "no"
}))

# generating new variables and keeping the old versions
mysamps.annot = lapply(mysamps.annot, function(x){
  x$medication__tolvaptan_old = x$medication__tolvaptan
  x$medication__tolvaptan = x$tolv_upd
  return(x)
})

# calculating the mayo class and keeping the old variable in different name
mysamps.annot2 = lapply(mysamps.annot, function(y){
  mayo_oo = as.vector(apply(y, 1, function(x){
    findMayo(mayoChart = newmayo, agem = "Age", agep = x["patient__age"],httkv = x["adpkd_progression__httkv"])
  }))
  y$adpkd_progression__mayo_org2 <- y$adpkd_progression__mayo
  y$adpkd_progression__mayo = NA
  y$adpkd_progression__mayo = mayo_oo
  
  return(y)
  })

# adding extra mayo class info
library(openxlsx)
cannot <- read.xlsx("~/Data/slope-models-adpkd/data/ADPKDComplete_v2.xlsx")
cannot$Mayo.class = gsub("^-$", NA, cannot$Mayo.class)
cannot = cannot[-grep("Typ 2",cannot$Mayo.class),] #remove mayo class 2 (n=1)

# one patient only had serum ID, which is converted to identifier ID
identifieridchange = which(cannot$identifier == cannot$Labornummer.Serum)
library(readxl)
Serumnumbers_oc <- read_excel("~/Data/slope-models-adpkd/data/Serumnumbers_oc.xlsx")
cannot[identifieridchange,"identifier"] = Serumnumbers_oc$identifier[match(cannot[identifieridchange,"Labornummer.Serum"],Serumnumbers_oc$`ID as in Plate`)]

extramayosfound = intersect(cannot$identifier[!is.na(cannot$Mayo.class)],
                            names(which(unlist(lapply(mysamps.annot2, function(x) 
                              all(is.na(x$adpkd_progression__mayo)))))))
extramayodf = cannot[grep(paste(extramayosfound, collapse = "|"), 
                          cannot$identifier),c("identifier", "Mayo.class")]
extramayodf$matching = paste0(extramayodf$identifier, "&", "baseline")
rownames(extramayodf) = extramayodf$identifier

mysamps.annot2 = lapply(mysamps.annot2, function(x){
  # mayo class filling from cannot SC will be done on this mayo groups
  x$adpkd_progression__mayo_org2filled = x$adpkd_progression__mayo 
  return(x)
})

sapply(extramayodf$identifier, function(x){
  indx = grep(extramayodf[x,"matching"], mysamps.annot2[[x]]$matching)
  mysamps.annot2[[x]][indx,"adpkd_progression__mayo_org2filled"] <<- extramayodf[x,"Mayo.class"]
  mysamps.annot2[[x]][indx,"adpkd_progression__mayo"] <<- extramayodf[x,"Mayo.class"]
})

mysamps.annot2 = lapply(mysamps.annot2, function(x){
  library(tidyr)
  x = x %>%
    #to fill NAs from top to down
    mutate(prev_val = (adpkd_progression__mayo)) %>%
    fill(prev_val, .direction = "down") %>%
    mutate(adpkd_progression__mayo = prev_val) %>%
    dplyr::select(-prev_val) %>%
    #to fill remaining NAs on top of the vector
    mutate(next_val = (adpkd_progression__mayo)) %>%
    fill(next_val, .direction = "up") %>%
    mutate(adpkd_progression__mayo = next_val) %>%
    dplyr::select(-next_val)
  return(x)
})

# if patient__eGFR is missing, remove that sample
removedsampleslope = list()
mysamps.annot2 = lapply(mysamps.annot2, function(x){
  remx = which(is.na(x$patient__eGFR))
  removedsampleslope[[paste(unique(x$patient_id))]] <<- ifelse(length(remx) > 0,x[remx,"matching"],NA)
  if (length(remx) > 0) {
    x = x[-remx,]
  }
  if (any(is.na(x$matchingHistcrea))) {
    x = x[-which(is.na(x$matchingHistcrea)),]
  }
  return(x)
}) # 1276

# after removing the ones that don't have name in previous piece of code,
# the ones that don't have any data is now being removed
mysamps.annot2[which(unlist(lapply(mysamps.annot2, nrow)) == 0)] = NULL # 1265

# finding how many samples from how many patients were removed due to interventions
extrainf = mysamps.annot2
lapply(names(extrainf), function(x){
  extrainf[[x]]$difdates <<- NA
  extrainf[[x]]$difdates_ext <<- NA

  extrainf[[x]]$difdates <<- as.numeric(extrainf[[x]]$exam_date - extrainf[[x]]$start4515orginal)
  extrainf[[x]]$difdates_ext <<- as.numeric(extrainf[[x]]$exam_date - extrainf[[x]]$DATES)

  extrainf[[x]]$ddcons <<- NA
  extrainf[[x]]$rems <<- NA

  extrainf[[x]]$ddcons[which(extrainf[[x]]$difdates >= 0)] <<- "tolvaptan"
  extrainf[[x]]$ddcons[which(extrainf[[x]]$difdates_ext >= 0)] <<- "others"

  indx = which(!is.na(extrainf[[x]]$ddcons))
  if (all(is.na(extrainf[[x]]$ddcons))) {
    extrainf[[x]]$rems <<- "no removal"
  }
  else if (extrainf[[x]]$ddcons[indx[1]] == "tolvaptan") {
    extrainf[[x]]$rems[indx] <<- "tolvaptan"
  }
  else if (extrainf[[x]]$ddcons[indx[1]] == "others") {
    extrainf[[x]]$rems[indx] <<- "others"
  }
  else{
    extrainf[[x]]$rems <<- NA
  }
})

removalscoll = do.call(rbind.data.frame, lapply(extrainf, function(x){
  cbind.data.frame(length(which(x$rems == "tolvaptan")),
                   length(which(x$rems == "others")))
}))

removalscoll$cmb = apply(removalscoll, 1, sum)

names(commtfix) = problems_mct # both defined in SC_02_medication_tolv
checthemso = rownames(removalscoll)[which(removalscoll$`length(which(x$rems == "others"))` != 0)]

paste0("n removed samples (total)=", sum(removalscoll$cmb))
paste0("from n patients=", length(which(removalscoll$cmb != 0)))

paste0("n removed samples (tolvaptan)=", sum(removalscoll[which(removalscoll$`length(which(x$rems == "tolvaptan"))` > 0),1]))
paste0("from n patients=", length(which(removalscoll$`length(which(x$rems == "tolvaptan"))` > 0)))


# slope calculation
merged.creas.final2 = lapply(mysamps.annot2, function(x){
  # ordering the samples according to date
  x = x[order(x$exam_date),]
  
  # whether a patient has used intervention
  x$tolv.ever = ifelse(any(na.omit(x$medication__tolvaptan) == 2), "yes", "no")

  # calculating slopes for each group
  # # assigning groups (whether they started and stopped and started again)
  x$groups <- as.factor(tuvele(x))
  
  # # calculating slopes for each group
  x.splt = split(x, f = x$groups)
  x <- do.call(rbind.data.frame, lapply(split(x, f = x$groups), get_coef_new))
  
  # to select which samples are before, on or after intervention
  x$forwd = rolling.mean(x$medication__tolvaptan)#, rev = T)
  x$revs = rev(rolling.mean(rev(x$medication__tolvaptan)))
  
  x$bTF = x$medication__tolvaptan == 1 & x$forwd == 1 & x$revs == 1 | x$medication__tolvaptan == 1 & x$forwd == 1 & x$revs == 2
  x$oTF = x$medication__tolvaptan == 2 & x$forwd == 2
  x$aTF = x$medication__tolvaptan == 1 & x$revs == 1 & x$forwd == 2
  
  # calculate the mean of the slopes
  x$before = mean(x$slope[which(x$bTF)], na.rm = T)
  x$ontolv = mean(x$slope[which(x$oTF)], na.rm = T)
  x$after = mean(x$slope[which(x$aTF)], na.rm = T)
  return(x)
})

merged.creas.final3 = do.call(rbind.data.frame, merged.creas.final2) # 8016   60
merged.creas.final3$medication__tolvaptan = as.factor(merged.creas.final3$medication__tolvaptan)
rownames(merged.creas.final3) = NULL

# which patients have double mayo class information?
# because several MRI measurements were performed for some patients at 
# different time points leading to more than one mayo class
doublemayos = merged.creas.final3 %>%  group_by(patient_id) %>%
  reframe(unique = paste(unique(adpkd_progression__mayo), collapse = " /////"),
          uniqcount = length(na.omit(unique(adpkd_progression__mayo))))

# extracting the first calculated mayo class per patient
wbcdoublemayo = merged.creas.final3 %>% group_by(patient_id) %>% 
  filter(!is.na(adpkd_progression__httkv)) %>% slice_min(exam_date, with_ties = F)

# saving the original mayo classes
merged.creas.final3$adpkd_progression__mayo_org2filled_duprem = merged.creas.final3$adpkd_progression__mayo

# replacing all mayo classes with the first calculated mayo class per patient
invisible(sapply(unique(merged.creas.final3$patient_id), function(x){
  if (!all(is.na(wbcdoublemayo[grep(x, wbcdoublemayo$patient_id),"adpkd_progression__mayo"]))) {
    selctmayo = rep(pull(wbcdoublemayo[grep(x, wbcdoublemayo$patient_id),],"adpkd_progression__mayo"), times = length(grep(x, merged.creas.final3$patient_id)))
    merged.creas.final3[grep(x, merged.creas.final3$patient_id),"adpkd_progression__mayo"] <<- selctmayo 
  }
}))

extra_egfr = merged.creas.final3[grep("identifier", merged.creas.final3$patient_id),] %>%
  group_by(patient_id) %>% reframe("usedbs" = length(which(bTF == T)))

paste0("n removed eGFR measurement (<3 measurement)=", sum(extra_egfr[which(extra_egfr$usedbs > 0 & extra_egfr$usedbs < 3),"usedbs"]))
# 601 eGFR measurement were removed

paste0("from n patients=", length(which(extra_egfr$usedbs > 0 & extra_egfr < 3)))
# from 477 patients

paste("############")

paste0("n eGFR measurement (total, used for before slope)=", 
       sum(extra_egfr[which(extra_egfr$usedbs > 2),"usedbs"]))
# 5371 eGFR measurement were used

paste0("from n patients=", length(which(extra_egfr$usedbs > 2))) 
# from 578 patients were used for slope calculation


selannoth = c("matching", "merged_creat", "patient__age", "patient__eGFR","exam_date",
              "patient__CKD", "tolv.ever", "groups", "slope", "rmse", "ontolv",
              "before", "after","bTF","oTF","aTF","adpkd_progression__mayo_org",
              "adpkd_progression__mayo_org2", "adpkd_progression__mayo_org2filled",
              "adpkd_progression__mayo_org2filled_duprem", "adpkd_progression__mayo")
annot.half2 = merged.creas.final3[,selannoth] # 8016   21
colnames(annot.half)[grep("_mayo", colnames(annot.half))] = "adpkd_progression__mayo_old"
annot.half2 = merge(annot.half, annot.half2, by = "matching", all.x = T)

annot.half = annot.half2 # 2940  226

# removing patients that dont have eGFR
annot.half = annot.half[-which(is.na(annot.half$patient__eGFR)),] 
# 2863  226

# n patient=1102
paste0("npatient=", length(unique(annot.half[grep("identifier",annot.half$patient_id),"patient_id"])))
# n observation=2577
paste0("nobservation=", length((annot.half[grep("identifier",annot.half$patient_id),"patient_id"])))
# n followup range
paste0("nfollowupRange=[", max(table(grep("identifier", annot.half$patient_id, value = T))),"-", # 9
       min(table(grep("identifier", annot.half$patient_id, value = T))), "]") # 1

paste0("##############")

# n patient=1102
paste0("npatient=", length(unique(merged.creas.final3[grep("identifier",merged.creas.final3$matchingHistcrea),"patient_id"])))
# n observation=7730
paste0("neGFRmeasurement=", length((merged.creas.final3[grep("identifier",merged.creas.final3$matchingHistcrea),"patient_id"])))
# n followup range
paste0("nfollowupRange=[", max(table((merged.creas.final3[grep("identifier",merged.creas.final3$matchingHistcrea),"patient_id"]))),"-", # 9
       min(table((merged.creas.final3[grep("identifier",merged.creas.final3$matchingHistcrea),"patient_id"]))), "]") # 1



mutations_ndb <- read.csv("~/Data/ADPKD_Cohort_ClinicalID_upd/clean_202310201230.csv")
mutations_ndb$patient_id = mutations_ndb$title
annot.half = merge(x = mutations_ndb, y = annot.half, all.y = T, by = "patient_id") # 2863 245

propkdall <- read.csv("~/Data/ADPKD_Cohort_ClinicalID_upd/automation_plots_202310231040.csv")
propkdall$automation_exam_date = propkdall$exam_date
propkdall$exam_date = NULL
propkdall$matching = paste(propkdall$patient_id,propkdall$examination, sep = "&")
propkdall[,c("patient_id", "examination")] = NULL
annot.half = merge(propkdall, annot.half, by = "matching", all.y = T) # 2863  325

annot.half$matching_SC = paste0(annot.half$patient_id, "&", annot.half$exam_date)

# patients who are using octreotid
octreotid = "xxxxxxxxxxxxxxxx" # IDs are masked
annot.half$octreotid = F
annot.half[grep(paste(octreotid, collapse = "|"), annot.half$patient_id, ignore.case = T),"octreotid"] = T

