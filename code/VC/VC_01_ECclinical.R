# CLINICAL DATA OF EXTERNAL VALIDATION ----

source("~/Data/slope-models-adpkd/code/00functions.R", echo = TRUE)
library(haven)
library(labelled)
library(readxl)

# | # Preparing the Data ----

# loading the external clinical data
external <- haven::read_sav("~/Data/slope-models-adpkd/data/DIPAK_Obs_Cologne.sav") # 696 932
external[grep("xxxxxxxxxxxxxxxx", external$stdnr),c("date_lab_serum_2","serum_creatinine_2")] = NA # IDs are masked

# df containing column name and their labels
lbels.ext = cbind.data.frame("names" = names(unlist(lapply(external, function(x) attr(x,"label")))),
                              "labels" = as.character(unlist(lapply(external, function(x) attr(x,"label")))))

# replacing empty cells with NA
emptyStrings = names(which(apply(external, 2, 
                                 function(x) length(which(x == "" | x == " ") > 0)) != 0))
external[,emptyStrings] = apply(external[,emptyStrings], 2, 
                                 function(x) ifelse(x == "", NA, x)) # 696 932

# removal of mayo class 2 (28 were atypical)
remmy2 = which(external$mayo_typical_atypical == 2)
mayoclass2G = external$stdnr[remmy2]
external = external[-remmy2,] # 668 932

# extracting medical and medication history 
medicalhistNmedic = external[,c(1:5,224, # 224 is hypertension
                                 grep("medic_|medical_history", 
                                      colnames(external)))] # 668 714 

# extracting required annotations (the variables used in the models) 
externalreq = external[,c(grep("serum_creatinine|date_lab_serum|mutation_category", 
                                 colnames(external), value = T),
                            "stdnr","sex", "dateofbirth", "intnr", "ethnicity", 
                            "hypertension_binary")] # 668 30

# | # | # Serum creatinine and eGFR calculations  ----
library(reshape2)
# "date_lab_serum" and "serum_creatinine" variables are 
# converted from wide to long format
externalreqM = melt(externalreq, id.vars = colnames(externalreq)[23:30]) 
# 14696    13

# extracting the follow-up number (1 to 11)
externalreqM$nums = gsub(".*serum_|.*creatinine_","",externalreqM$variable)

# adding "0"s in front of the numbers which are from 1-9
externalreqM$nums = ifelse(nchar(externalreqM$nums) < 2, 
                            paste0("0", externalreqM$nums), externalreqM$nums)

# merging the sample ID with follow-up number
externalreqM$tpID = paste0(externalreqM$stdnr, "_m", externalreqM$nums)

# column IDs are generated for "date_lab_serum" and "serum_creatinine"
externalreqM$colnam = gsub("_([0-9])|_([0-9])([0-9])","",externalreqM$variable) 
#14696    13

# now "date_lab_serum" and "serum_creatinine" are each in one column 
externalLong = dcast(externalreqM, 
                      tpID + stdnr + dateofbirth + sex + hypertension_binary + 
                        intnr + ethnicity + first_mutation_category + 
                        second_mutation_category~colnam) # 7348   11

rownames(externalLong) = externalLong$tpID
externalLong = as.data.frame(to_factor(externalLong))

# extracting the dates and converting it into date format
alldates = grep("date", colnames(externalLong))
externalLong[,alldates] = lapply(externalLong[,alldates], as.Date, 
                                  origin = "1970-01-01")
# converting the unit of serum creatinine
externalLong$serum_creatinine = as.numeric(externalLong$serum_creatinine)
externalLong$serum_creatinine = 0.0113*externalLong$serum_creatinine

# calculating age
externalLong$age = round(as.numeric(difftime(externalLong$date_lab_serum,
                                              externalLong$dateofbirth, 
                                              units = "days")/365.25), 
                          digits = 2)
# calculating eGFR
externalLong$patient__eGFR = computeGFR_upd(externalLong, dob = "age",
                                             creatinine = "serum_creatinine",
                                             gender = "sex") # 7348   13

# | # | # Interventions ----
# | # | # | # Transplant and Dialysis ----
library(reshape2)
# extracting only medical history for checking whether the 
# patient underwent kidney transplant 
transplant = medicalhistNmedic[,c(1,grep("medical_history_name|medical_history_startdate|medical_history_stopdate", 
                                         colnames(medicalhistNmedic)))]
# converting from wide to long format 
transplant = melt(transplant, id = "stdnr")

# column IDs are generated for "medical_history_name", 
# "medical_history_startdate", "medical_history_stopdate" 
transplant$colname = gsub("_([0-9]+)$","",transplant$variable)

# extracting the follow-up number
transplant$nm = gsub(".*_","",transplant$variable)

# merging the sample ID with follow-up number
transplant$rowname = paste0(transplant$stdnr, "#", transplant$nm)

# extracting only the patients who underwent kidney transplant and dialysis
# none of the patients underwent nephrectomy
transplantdone = transplant[grep("kidneytranspl|kidney transpl|dialysis", 
                                 transplant$value, ignore.case = T),"rowname"]
transplantdone = unique(transplantdone)

# extracting the full medical information regarding those interventions
transplant2 = transplant
transplant2 = transplant2[grep(paste(paste0("^",transplantdone,"$"), 
                                     collapse = "|"),transplant2$rowname),]

transplantbenz = dcast(transplant2, rowname+stdnr~colname)

# generating a column in the data for intervention
externalLong$transplant_dialysis = NA

#extracting those dates from transplantbenz and writing them on externalLong
invisible(lapply(transplantbenz$stdnr, function(x){
  externalLong[externalLong$stdnr == x,"transplant_dialysis"] <<- as.character(transplantbenz[transplantbenz$stdnr == x,"medical_history_startdate"])
})) #7348   14


# | # | # | # Tolvaptan and Octeroid  ----
library(reshape2)
# extracting medication usage
tolvdec = medicalhistNmedic[,c(1,grep("medic_name|medic_start_date|medic_stop_date|medic_baseline_use", 
                                      colnames(medicalhistNmedic)))]
# converting from wide to long format 
tolvdec2 = melt(tolvdec, id = "stdnr")

# column IDs are generated for "medic_name", "medic_baseline_use", 
# "medic_start_date", "medic_stop_date"
tolvdec2$colname = gsub("_([0-9]+)$","",tolvdec2$variable)

# extracting the follow-up number
tolvdec2$nm = gsub(".*_","",tolvdec2$variable)

# merging the sample ID with follow-up number
tolvdec2$rowname = paste0(tolvdec2$stdnr, "#", tolvdec2$nm)

# | # | # | # | # Octreoid ----
# extracting patients with octreoid usage
octreotidext = tolvdec2[grep("octreotide|lanreotid|pasireotid|Sandostatin|mytolac|somatolina|somatoline|signifor", 
                              tolvdec2$value, ignore.case = T),"rowname"]

# extracting IDs of those patients
octreotidext = unique(gsub("#.*","",octreotidext))
octreotidext = substring(gsub("\\.|-","",octreotidext),2)

# | # | # | # | # Tolvaptan ----
# | # | # | # | # | # Usage started after the study  ----
# extracting patients with tolvaptan usage
tolvusingp = tolvdec2[grep("tolvaptan|tolv", tolvdec2$value, ignore.case = T),"rowname"]
tolvusingp = unique(tolvusingp)

# extracting the full medical information regarding those interventions
tolvgen = tolvdec2[grep(paste(paste0("^",tolvusingp,"$"), collapse = "|"),tolvdec2$rowname),]  

# extracting IDs of those patients
tolvbenz = dcast(tolvgen, rowname~colname)
tolvbenz$stdnr = gsub("#.*","",tolvbenz$rowname)

# | # | # | # | # | # Usage started before the study  ----
library(dplyr)
#extracting IDs of those patients
usingfromstart = unique(tolvbenz[which(tolvbenz$medic_baseline_use == "1"),"stdnr"])

# assigning start dates for those patients (actual dates were not available)
dateFixingTolv = externalLong[grep(paste(usingfromstart,collapse = "|"),externalLong$stdnr),c("stdnr","date_lab_serum")]
dateFixingTolv = dateFixingTolv %>% group_by(stdnr) %>% 
  reframe(startDate = min(date_lab_serum,na.rm = T))
dateFixingTolv = as.data.frame(dateFixingTolv)
invisible(lapply(usingfromstart, function(x){
  tolvbenz[tolvbenz$stdnr == x,"medic_start_date"] <<- as.character(dateFixingTolv[dateFixingTolv$stdnr == x,"startDate"])
}))

# converting it to a date
tolvbenz$medic_start_date = as.Date(tolvbenz$medic_start_date)

# | # | # | # Merging intervention dates with samples ----
library(dplyr)

#taking the minimum date of start ()
tolvbenz = tolvbenz %>% group_by(stdnr) %>% 
  dplyr::mutate(startDate = min(medic_start_date,na.rm = T))

# extracting the start date information
tolv_end = tolvbenz[,c("stdnr", "startDate")]
tolv_end = as.data.frame(unique(tolv_end))

# writing the start date information into main annotation data
externalLong$tolvStart = NA
invisible(lapply(tolv_end$stdnr, function(x){
  externalLong[externalLong$stdnr == x,"tolvStart"] <<- as.character(tolv_end[tolv_end$stdnr == x,"startDate"])
}))

# merging the start dates of all interventions and naming it "tolvStartTranspl"
externalLong = externalLong %>% 
  mutate(tolvStartTranspl = coalesce(tolvStart,transplant_dialysis))

# generating a new variable "medication__tolvaptan" and 
# the patients who dont have a start date are assigned to "no"
externalLong$medication__tolvaptan = NA
externalLong[which(is.na(externalLong$tolvStartTranspl)),"medication__tolvaptan"] = "no"

# calculating the time difference between the serum creatinine measurement
# and intervention start date
externalLong$diff = as.numeric(difftime(externalLong$tolvStartTranspl,
                                         externalLong$date_lab_serum, 
                                         units = "days"))

# if a patient have time difference 
invisible(lapply(unique(externalLong$stdnr), function(x){ #patient-wise
  indx = which(externalLong$stdnr == x) # indexes of all samples 
  
  # indexes of the negative time difference were selected because 
  # if it is negative, then it means the measurement were collected after the 
  # start date of the intervention
  indxn = indx[which(externalLong$diff[indx] <= 0)]
  if (length(indxn) > 0) { # and if there are negative time difference
    # usage is assigned to "yes"
    externalLong[indxn,"medication__tolvaptan"] <<- "yes" 
  }
  # and if the are positive time difference
  indxp = indx[which(externalLong$diff[indx] > 0)] 
  if (length(indxp) > 0) { # usage is assigned to "no"
    externalLong[indxp,"medication__tolvaptan"] <<- "no"
  }
}))
# length(unique(externalLong$stdnr)) # 668 patient

# | # | # Subsetting clinical data according to proteome -----
library(stringr)
# only some patients had proteome data, so extracting clinical 
# annotation of those patients only
cannot = externalLong[grep("xxxxxxxxxxxxxxxx", externalLong$stdnr),] #2530 18 # IDs are masked
#230 patients

# generating proteome ID matches by extracting only numbers from the patient ID
cannot$protmatch = gsub("\\.", "", str_sub(cannot$stdnr, start = 2, end = 6))

# generating examination date variable 
cannot$exam_date = cannot$date_lab_serum

# removing samples which don't have intervention data or eGFR measurement
rem.cannot = which(is.na(cannot$patient__eGFR) | is.na(cannot$medication__tolvaptan))
cannot = cannot[-rem.cannot,] # 1481 20 # 228 patients

# | # | # Slope Calculation ----
cannot.splt = split(cannot, f = cannot$stdnr)
cannot.splt = lapply(cannot.splt, function(x){
  # ordering the samples according to date
  x = x[order(x$exam_date),]
  
  # whether a patient has used tolvaptan
  x$medication__tolvaptan_s = ifelse(is.na(x$medication__tolvaptan), NA,
                                     ifelse(x$medication__tolvaptan == "yes", 
                                            2, 1))
  
  # calculating slopes for each group
  # # assigning groups (whether they started and stopped and started again)
  x$groups <- as.factor(tuvele(x))
  
  # # calculating slopes for each group
  x.splt = split(x, f = x$groups)
  x <- do.call(rbind.data.frame, lapply(split(x, f = x$groups), get_coef_new))
  
  # to select which samples are before, on or after Tolvaptan
  x$forwd = rolling.mean(x$medication__tolvaptan_s)
  x$revs = rev(rolling.mean(rev(x$medication__tolvaptan_s)))
  
  x$bTF = x$medication__tolvaptan_s == 1 & x$forwd == 1 & x$revs == 1 | x$medication__tolvaptan_s == 1 & x$forwd == 1 & x$revs == 2
  x$oTF = x$medication__tolvaptan_s == 2 & x$forwd == 2
  x$aTF = x$medication__tolvaptan_s == 1 & x$revs == 1 & x$forwd == 2
  
  # calculate the mean of the slopes
  x$before = mean(x$slope[which(x$bTF)], na.rm = T)
  x$ontolv = mean(x$slope[which(x$oTF)], na.rm = T)
  x$after = mean(x$slope[which(x$aTF)], na.rm = T)
  return(x)
})


cannot.splt2 <- do.call(rbind.data.frame, cannot.splt)
cannot.splt2 = cannot.splt2 %>% group_by(stdnr) %>% 
  dplyr::mutate(count = length(na.omit(which(before == slope & bTF == T)))) 
# 1481 33 && # 228 unique patients

# removal due to tolvaptan
length(which(cannot.splt2$bTF == F & !is.na(cannot.splt2$tolvStart))) 
length(unique(cannot.splt2[which(cannot.splt2$bTF == F & !is.na(cannot.splt2$tolvStart)),]$stdnr))
# 258 eGFR measurement were removed from 55 patient

# removal due to dialysis, nephrectomy, transplants
length(which(cannot.splt2$bTF == F & !is.na(cannot.splt2$transplant_dialysis)))
# 0 

# removal due to < 3 eGFR measurement
lt3grn = cannot.splt2[which(cannot.splt2$bTF == T),] # 1223 33 - counted as bTF
lt3grn2 = cannot.splt2[which(cannot.splt2$bTF == T & !is.na(cannot.splt2$before)),] # 1214   33
# # 1223-1214 = 9 patient had one measurement

length(which(lt3grn2$count == 2)) # 20
length(unique(lt3grn2[which(lt3grn2$count == 2),]$stdnr)) # 10
# # 20 measurement from 10 patient (two measurement per patient)

lt3grn3 = lt3grn2[-which(lt3grn2$count < 3),] #1194   33
length(unique(lt3grn3$stdnr)) # reamining 180 patient

# | # | # htTKV calculation ----
# | # | # | # Preparing height information ----
# extracting height information
heightext = external[,c(grep("dateofbirth|stdnr|height", colnames(external)))] # 668 30

# converting height from m to cm 
heightext[,-c(1:2)] = lapply(heightext[,-c(1:2)], function(x){
  ifelse(is.na(x), NA, ifelse(x < 3, x*100, x))
}) 

# if height is < 100 cm or >300 cm, assigned them to NA
heightext <- heightext %>%
  mutate_if(is.numeric, ~ ifelse(. < 100, NA, .)) %>% 
  mutate_if(is.numeric, ~ ifelse(. > 300, NA, .))

# calculating patient-wise median height 
heightext$heightav = apply(heightext[,-c(1:2)], 1, median, na.rm = T)

# converting height from cm to m for htTKV calculation
heightext$heightavm = heightext$heightav/100

# merge the data with annotation data
externalHeight = left_join(external, heightext[,c("stdnr","heightavm")], 
                            by = "stdnr")

# | # | # | # Preparing TKV information ----
# extracting htTKV related measurements
mriRel = externalHeight[,grep("stdnr|sex|dateofbirth|intnr|ethnicity|mri|right|left|heightavm", 
                               colnames(externalHeight))]
# following measurements do not exist
mriRel$date_mri_4 = NULL
mriRel$date_mri_5 = NULL

library(reshape2)
# converting from wide to long format 
mriRel.m = melt(mriRel, id = c("stdnr","sex","dateofbirth","intnr",
                               "ethnicity","heightavm"))

# column IDs are generated for "date_mri", "right_kidney_volume", 
# "left_kidney_volume" 
mriRel.m$var = gsub("_(\\d+)$", "", mriRel.m$variable)

# extracting the follow-up number
mriRel.m$tp = gsub(".*_", "tp", mriRel.m$variable)

# merging the sample ID with follow-up number
mriRel.m$idtp = paste0(mriRel.m$stdnr, "&", mriRel.m$tp)

mriRel.m.d = dcast(idtp+stdnr+sex+dateofbirth+intnr+ethnicity+heightavm~var, 
                   data = mriRel.m)
mriRel.m.d$date_mri = as.Date(mriRel.m.d$date_mri, origin = "1970-01-01")

# calculating TKV
mriRel.m.d$tkv = as.numeric(mriRel.m.d$left_kidney_volume) +
  as.numeric(mriRel.m.d$right_kidney_volume)

# merging the sample ID with date of the MRI
mriRel.m.d$matching = paste0(mriRel.m.d$stdnr, "&", mriRel.m.d$date_mri)
mriRel.m.d$matching[grep("&NA$", mriRel.m.d$matching)] = NA

# | # | # | # Preparing htTKV information ----
cannot.splt5 = cannot.splt2

# merging the sample ID with date of the serum creatinine measurement
cannot.splt5$matching = paste0(cannot.splt5$stdnr, "&",
                               cannot.splt5$date_lab_serum) # 1481   33
cannot.splt5$matching[grep("&NA$", cannot.splt5$matching)] = NA

# merging the TKV information with annotation data
cannot.splt5 = left_join(cannot.splt5, 
                         mriRel.m.d[,c("tkv", "matching")], 
                         by = "matching") #1481   35
# removing the samples which don't have MRI date, 
# left/right kidney volume, and height
mriRel.m.d2 = mriRel.m.d[(!is.na(mriRel.m.d$date_mri) & !is.na(mriRel.m.d$left_kidney_volume) & !is.na(mriRel.m.d$right_kidney_volume) & !is.na(mriRel.m.d$heightavm)),]

# generating patient IDs
mriRel.m.d2$stdnr = gsub("&.*","",mriRel.m.d2$idtp)

# calculating age at the time of MRI measurement
mriRel.m.d2$ageR = round(as.numeric(difftime(mriRel.m.d2$date_mri,
                                             mriRel.m.d2$dateofbirth,
                                             units = "days")/365.25), 
                         digits = 0)
# calculating htTKV
mriRel.m.d2$httkv = mriRel.m.d2$tkv/mriRel.m.d2$heightavm

# | # | # Preparing MAYO class information ----
# calculating MAYO class
mriRel.m.d2$mayo = apply(mriRel.m.d2, 1, function(x){
  findMayo(mayoChart = newmayo, agem = "Age",
           agep = x["ageR"],
           httkv = x["httkv"])
})

# which patients have double mayo class information?
# because several MRI measurements were performed for some patients at 
# different time points leading to more than one mayo class
doublemayosG = mriRel.m.d2 %>%  group_by(stdnr) %>%
  reframe(unique = paste(unique(mayo), collapse = " /////"),
          uniqcount = length(na.omit(unique(mayo))))

# extracting the first calculated mayo class per patient
wbcdoublemayoG = mriRel.m.d2 %>% group_by(stdnr) %>% 
  filter(!is.na(httkv)) %>% slice_min(date_mri, with_ties = F)

# saving the original mayo classes
mriRel.m.d2$mayo_org2filled_duprem = mriRel.m.d2$mayo

# replacing all mayo classes with the first calculated mayo class per patient
invisible(sapply(unique(mriRel.m.d2$stdnr), function(x){
  if (!all(is.na(wbcdoublemayoG[grep(x, wbcdoublemayoG$stdnr),"mayo"]))) {
    selctmayo = rep(pull(wbcdoublemayoG[grep(x, wbcdoublemayoG$stdnr),],"mayo"),
                    times = length(grep(x, mriRel.m.d2$stdnr)))
    mriRel.m.d2[grep(x, mriRel.m.d2$stdnr),"mayo"] <<- selctmayo 
  }
}))

mayoclassesReady = unique(mriRel.m.d2[,c("stdnr","mayo")])

# converting tibble to df
cannot.splt3 = as.data.frame(cannot.splt5) # 1481 35 # 228 patients

# merging the new mayo classes with annotation data
cannot.long = left_join(cannot.splt3, mayoclassesReady, by = "stdnr") 
# 1481 36 # 228 unique patient

# adjusting th levels of mayo class
cannot.long$mayo = as.factor(cannot.long$mayo)
cannot.long$mayoClass = cannot.long$mayo
levels(cannot.long$mayoClass) = c("1A" = "Less", "1B" = "Less", 
                                  "1C" = "Mid", "1D" = "More", "1E" = "More")

# | # | # Preparing Mutation information ----
# adjusting the levels of PKD1 and PKD2 genotype
cannot.long$first_mutation_category[cannot.long$first_mutation_category == "not performed"] = NA
cannot.long$first_mutation_category = droplevels(cannot.long$first_mutation_category)
levels(cannot.long$first_mutation_category) = list(unknown = "no mutation detected", 
                                                   unknown = "other",
                                                   PKD2 = "PKD2 truncating", 
                                                   PKD2 = "PKD2 non-truncating",
                                                   PKD1vus = "PKD1, not possible to decide truncating/non-truncating",
                                                   PKD1NT = "PKD1 non-truncating",
                                                   PKD1T = "PKD1 truncating")

cannot.long$second_mutation_category = droplevels(cannot.long$second_mutation_category)
levels(cannot.long$second_mutation_category) = list(PKHD1 = "PKHD1",
                                                    PKD2 = "PKD2 non-truncating",
                                                    PKD1vus = "PKD1, not possible to decide truncating/non-truncating",
                                                    PKD1NT = "PKD1 non-truncating",
                                                    PKD1T = "PKD1 truncating")
# 1481 37 # 228 unique patients
