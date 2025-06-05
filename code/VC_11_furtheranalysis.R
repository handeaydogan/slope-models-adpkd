load("~/Data/slope-models-adpkd/data/results/data_Validation.RData")
library(dplyr)
# Threshold for Slope: -3 ----
# # Creating binary variable for slope ≥ -3 or < 3 ml/min per year
# | # Proteome Model ----
cmbnData$predBinary = NA
cmbnData$predBinary[which(cmbnData$predSLOPEs >= -3)] = "EG"
cmbnData$predBinary[which(cmbnData$predSLOPEs < -3)] = "L"
cmbnData$predBinary = as.factor(cmbnData$predBinary)

cmbnData$beforeBinary = NA
cmbnData$beforeBinary[which(cmbnData$before >= -3)] = "EG"
cmbnData$beforeBinary[which(cmbnData$before < -3)] = "L"
cmbnData$beforeBinary = as.factor(cmbnData$beforeBinary)

cmbnData$age = coalesce(cannot.before.tolv[rownames(cmbnData),"merged_age"], 
                        metaclin2_slopes[rownames(cmbnData),"age"])

# | # Combined Model ----
cmbnDataExtra$predBinary = NA
cmbnDataExtra$predBinary[which(cmbnDataExtra$predSLOPEs >= -3)] = "EG"
cmbnDataExtra$predBinary[which(cmbnDataExtra$predSLOPEs < -3)] = "L"
cmbnDataExtra$predBinary = as.factor(cmbnDataExtra$predBinary)

cmbnDataExtra$beforeBinary = NA
cmbnDataExtra$beforeBinary[which(cmbnDataExtra$before >= -3)] = "EG"
cmbnDataExtra$beforeBinary[which(cmbnDataExtra$before < -3)] = "L"
cmbnDataExtra$beforeBinary = as.factor(cmbnDataExtra$beforeBinary)

cmbnDataExtra$age = coalesce(cannot.before.tolv[rownames(cmbnDataExtra),"merged_age"], 
                             metaclin2_slopes[rownames(cmbnDataExtra),"age"])

# | # MIC Model ----
metaclin2_slopes$identifier = metaclin2_slopes$identifiers

# future eGFR calculation 
timediffcp = list()
lapply(metaclin2_slopes$Name, function(x) timediffcp[[paste(x)]] <<- 1:10)

mayoFuturecp = mapply(function(x,y){
  calcFuteGFR(timedif = x, id = y, idcol = "Name" , 
              data = metaclin2_slopes, which = "mayo slopes")
}, x = timediffcp, y = names(timediffcp))
mayoFuturecp2 = as.data.frame(mayoFuturecp)
mayoFuturecp2 = as.list(mayoFuturecp2)

# converting them to a slope
mayoSlopescp = mapply(function(eGFR, year){
  if (all(is.na(eGFR))) {
    NA
  }
  else if (length(eGFR) > 1) {
    lm(eGFR~year)$coefficients[2]
  }
}, eGFR = mayoFuturecp2, year = timediffcp[names(mayoFuturecp2)])

# generating the df for MIC slopes
mayoSlopescp.df = cbind.data.frame("mayoSlope" = as.numeric(unlist(mayoSlopescp)), 
                                   "ID" = gsub("\\..*", "", names(mayoSlopescp)),
                                   metaclin2_slopes[match(gsub("\\..*", "",names(mayoSlopescp)), 
                                                          metaclin2_slopes$Name),
                                                    c("before", "identifier", "age", 
                                                      "eGFR", "gender", "mayo", 
                                                      "mayo5")],
                                   "date" = metaclin2_slopes[match(gsub("\\..*", "", names(mayoSlopescp)), 
                                                                   metaclin2_slopes$Name),"merged_dates"])
mayoSlopescp.df$ProjID = gsub("_.*", "", rownames(mayoSlopescp.df))

# merging it with SC's MIC slopes
mayoSlopes.df$ProjID = "sc"
rownames(mayoSlopes.df) = gsub("\\..*", "", rownames(mayoSlopes.df))

ggplotSlopes2 = ggplotSlopes
ggplotSlopes2$batchnames = rownames(ggplotSlopes2)
rownames(ggplotSlopes2) = ggplotSlopes2$identifier

cmsmpsmayos = intersect(mayoSlopes.df$identifier, ggplotSlopes$identifier)
mayoSlopes.df$date = as.Date(NA)
mayoSlopes.df[cmsmpsmayos,"date"] = ggplotSlopes2[cmsmpsmayos,"date"]
mayoSlopes.df[cmsmpsmayos,"ID"] = ggplotSlopes2[cmsmpsmayos,"batchnames"]
rownames(mayoSlopes.df) = mayoSlopes.df$ID

clnms = c("ProjID", "mayoSlope", "ID", "before", "identifier", 
          "age", "gender", "mayo", "date")
mayoslopesall = rbind.data.frame(mayoSlopes.df[,clnms],mayoSlopescp.df[,clnms])

mayoslopesall$select1p1 = paste0(mayoslopesall$ProjID, "&",mayoslopesall$identifier, 
                                 "&", mayoslopesall$date)

mayoslopesall$predBinary = NA
mayoslopesall$predBinary[which(mayoslopesall$mayoSlope >= -3)] = "EG"
mayoslopesall$predBinary[which(mayoslopesall$mayoSlope < -3)] = "L"
mayoslopesall$predBinary = as.factor(mayoslopesall$predBinary)

mayoslopesall$beforeBinary = NA
mayoslopesall$beforeBinary[which(mayoslopesall$before >= -3)] = "EG"
mayoslopesall$beforeBinary[which(mayoslopesall$before < -3)] = "L"
mayoslopesall$beforeBinary = as.factor(mayoslopesall$beforeBinary)
mayoslopesall$eGFR = coalesce(cannot.before.tolv[match(rownames(mayoslopesall), cannot.before.tolv$identifier), "merged_eGFR"],
                              mayoSlopescp.df[match(rownames(mayoslopesall), rownames(mayoSlopescp.df)),"eGFR"])
#grep("identifier", rownames(mayoslopesall), invert = T, value = T)==rownames(mayoSlopescp.df)

# # 1 sample per patient --> per project
library(dplyr)
myslpsnarm = mayoslopesall[-which(is.na(mayoslopesall$mayoSlope)),]
matchindates = myslpsnarm %>% group_by(ProjID, identifier) %>%
  reframe(mindate = min(date), maxdate = max(date),
          matchmin = paste0(ProjID,"&",identifier, "&", mindate),
          matchmax = paste0(ProjID,"&",identifier, "&", maxdate))
minsel = pull(matchindates, "matchmin")
maxsel = pull(matchindates, "matchmax")
alldatesmyo = cmbnData$select1p1


library(pROC)
dfs <- list(cmbnData, cmbnDataExtra, mayoslopesall)
names(dfs) = c("cmbnData", "cmbnDataExtra", "mayoslopesall")
dfs <- lapply(dfs, function(x){
  x$batch = rownames(x)
  return(x)
})

m1 = left_join(dfs[[1]], dfs[[2]], by = "batch", suffix = c("_c", "_ce"))
m2 = left_join(m1, dfs[[3]], by = "batch")

for_roc = m2[,c("batch", "ProjID_c","before_c", "beforeBinary_c", 
                "predSLOPEs", "predBinary_c", "predSLOPEsexted", 
                "predBinary_ce", "select1p1_c", "mayoSlope")]

rownames(for_roc) <- for_roc$batch
for_roc$batch <- NULL
for_roc = na.omit(for_roc)

for_roc = for_roc[grep(paste0("^",minsel,"$", collapse = "|"), for_roc$select1p1_c),]
for_roc_s = split(x = for_roc, f = for_roc$ProjID_c)

mod_combn = roc(as.numeric(for_roc_s$"sc"$beforeBinary_c),
                for_roc_s$"sc"$predSLOPEsexted)
mod_mayo = roc(as.numeric(for_roc_s$"sc"$beforeBinary_c),
               for_roc_s$"sc"$mayoSlope)

aucmod_combn <- round(as.numeric(mod_combn$auc), 4)
aucmod_mayo <- round(as.numeric(mod_mayo$auc), 4)

library(ggplot2)
rocCurve = ggroc(list(mod_combn, mod_mayo), size = 1.1, legacy.axes = T) +
  scale_colour_manual(name = "Models", values = c("#165079", "#FE9D52"),
                      labels = c(paste0("Combined Model - AUC: ", 
                                        round(aucmod_combn,2)),
                                 paste0("MIC Model - AUC: ", 
                                        round(aucmod_mayo,2)))) +
  theme_bw() + 
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
               color = "darkgrey", linetype = "dashed") +
  theme(legend.key.height = unit(10, "mm"))

rocCurve
ggsave(filename = paste0(file, "FigureS9.png"), width = 7.48, height = 6.55)

# Followup times ----
# TableS1
# preparing the data
library(dplyr)
merged.creas.final3$psc = "no"
sdonesc = na.omit(merged.creas.final3[match(cannot.before.tolv$identifier, 
                                              merged.creas.final3$patient_id), 
                                        "patient_id"])
merged.creas.final3$psc[grep(paste(sdonesc, collapse = "|"), 
                               merged.creas.final3$patient_id)] = "yes"

merged.creas.final3$pitc = "no"
sdoneitc = na.omit(merged.creas.final3[match(grep("xxxxxx", cmbnData$identifier, # IDs are masked
                                                   invert = T, value = T), 
                                              merged.creas.final3$patient_id), 
                                        "patient_id"])
merged.creas.final3$pitc[grep(paste(sdoneitc, collapse = "|"), 
                               merged.creas.final3$patient_id)] = "yes"

summarySlopes = merged.creas.final3 %>%  group_by(patient_id) %>%
  reframe(totalMeas = n(), # total measurments
          beforeTolv = length(na.omit(which(medication__tolvaptan == 1 & forwd == 1 & revs == 1 | medication__tolvaptan == 1 & forwd == 1 & revs == 2))),
          psc = psc, pitc = pitc, # SC and ITC
          usedForBS = length(na.omit(which(before == slope))),
          exam_date = exam_date, bTF = bTF)

# For Screening Cohort
sum_datsc = summarySlopes %>% group_by(patient_id) %>% 
  filter(bTF == T) %>% filter(psc == "yes") 
avcalcsc = sum_datsc %>%
  reframe(difdate = (max(exam_date) - min(exam_date)),
          bsu = length(which(bTF)))
timedifmeassc = sum_datsc %>%
  reframe(lags = exam_date - lag(exam_date))

# For Internal/Temporal Cohort
sum_datitc = summarySlopes %>% group_by(patient_id) %>% 
  filter(bTF == T) %>% filter(pitc == "yes") 
avcalcitc = sum_datitc %>%
  reframe(difdate = (max(exam_date) - min(exam_date)),
          bsu = length(which(bTF)))
timedifmeasitc = sum_datitc %>%
  reframe(lags = exam_date - lag(exam_date))

# For External Cohort
gonids = metaclin2_slopes[metaclin2_slopes$ProjID == "ec","identifiers"]
frtbs = cannot.long[grep(paste(gonids, collapse = "|"), cannot.long$protmatch),]
sum_datext = frtbs %>% group_by(stdnr) %>% filter(bTF == T)
avcalcec = sum_datext %>%
  reframe(difdate = (max(date_lab_serum) - min(date_lab_serum)),
          bsu = length(which(bTF)))
timedifmeasec = sum_datext %>%
  reframe(lags = date_lab_serum - lag(date_lab_serum))


# merge all 
vector_vers = c(sc, itc, ec,
                median(as.numeric(avcalcsc$difdate)), # median followup time for SC
                median(as.numeric(avcalcitc$difdate)), # median followup time for ITC
                median(as.numeric(avcalcec$difdate)), # median followup time for EC
                round(mean(as.numeric(avcalcsc$difdate)), digits = 1), # mean followup time for SC
                round(mean(as.numeric(avcalcitc$difdate)), digits = 1), # mean followup time for ITC
                round(mean(as.numeric(avcalcec$difdate)), digits = 1), # mean followup time for EC
                median(as.numeric(avcalcsc$bsu)), # median number of eGFR measurement for SC
                median(as.numeric(avcalcitc$bsu)), # median number of eGFR measurement for ITC
                median(as.numeric(avcalcec$bsu)), # median number of eGFR measurement for EC
                median(na.omit(as.numeric(timedifmeassc$lags))), # median time btw these measurements for SC
                median(na.omit(as.numeric(timedifmeasitc$lags))), # median time btw these measurements for ITC
                median(na.omit(as.numeric(timedifmeasec$lags)))) # median time btw these measurements for EC

med_fu_df = as.data.frame(matrix(nrow = 3, vector_vers))
rownames(med_fu_df) = med_fu_df$V1
med_fu_df$V1 = NULL
colnames(med_fu_df) = c("med_fu_time", "mean_fu_time", "med_#_eGFR_meas", "med_time_btw_meas")
med_fu_df$med_fu_time_y = round(med_fu_df$med_fu_time/365.25, digits = 1)
med_fu_df$mean_fu_time_y = round(med_fu_df$mean_fu_time/365.25, digits = 1)
med_fu_df$med_time_btw_meas_y = round(med_fu_df$med_time_btw_meas/365.25, 
                                      digits = 1)

med_fu_df2 = med_fu_df
med_fu_df2$med_fu_time = paste0(med_fu_df$med_fu_time, " [", 
                                med_fu_df$med_fu_time_y, "]")
med_fu_df2$mean_fu_time = paste0(med_fu_df$mean_fu_time, " [", 
                                med_fu_df$mean_fu_time_y, "]")
med_fu_df2$med_time_btw_meas = paste0(med_fu_df$med_time_btw_meas, 
                                      " [", med_fu_df$med_time_btw_meas_y, "]")
med_fu_df2$med_time_btw_meas_y = NULL
med_fu_df2$med_fu_time_y = NULL
med_fu_df2$mean_fu_time_y = NULL


med_fu_df3 = t(med_fu_df2)
colnames(med_fu_df3) = c("Screening Cohort (SC)", 
                         "Internal/Temporal Cohort (ITC)", 
                         "External Cohort (EC)")

rownames(med_fu_df3) = c("Median follow-up time, days [years]", 
                         "Mean follow-up time, days [years]", 
                         "Median number of eGFR measurements", 
                         "Median time between measurements, days [years]")

write.csv(med_fu_df3, file = paste0(file, "TableS1.csv"))
