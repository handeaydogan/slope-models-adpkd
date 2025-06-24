# VALIDATION ----
# model created in SC and we are now testing it on ITC and EC
source("~/Data/slope-models-adpkd/code/VC/VC_02_mergedClinical.R", echo = TRUE)
source("~/Data/slope-models-adpkd/code/VC/VC_03_proteome_load.R", echo = TRUE)
source("~/Data/slope-models-adpkd/code/VC/VC_03_proteome_open.R", echo = TRUE)

# | # Preparing the Data ----
# selecting the ones that have before slopes
metaclin2_slopes = metaclin2[which(metaclin2$bTF == T & !is.na(metaclin2$before)),] # 631  37

# removing the ones that have slope bigger than 5 and smaller than -10
rem.metaclinslopes = which(metaclin2_slopes$before > 5 | metaclin2_slopes$before < -10) # 30
metaclin2_slopes = metaclin2_slopes[-rem.metaclinslopes,] # 601  37

# the ones that have only 2 measurement for before slopes were removed
bsused = merged.creas.final3 %>%  group_by(patient_id) %>% 
  reframe(bu = length(na.omit(which(before == slope & bTF == T))))
wbstayin = setdiff(metaclin2_slopes$identifiers,bsused$patient_id[which(bsused$bu < 3)])
metaclin2_slopes = metaclin2_slopes[grep(paste(wbstayin, collapse = "|"),
                                         metaclin2_slopes$identifiers),] #589  37

metaclin2_slopes = metaclin2_slopes[-which(metaclin2_slopes$count < 3),] 
# 581  37 #8 removed

# selecting common samples between proteome and annotation file
instsampprotannot = intersect(colnames(data.slopesF),rownames(metaclin2_slopes))
data.slopesF = data.slopesF[,instsampprotannot]
metaclin2_slopes = metaclin2_slopes[instsampprotannot,]

#n patient
paste0("npatientitc=", length(unique(grep("identifier", metaclin2_slopes$identifiers, value = T)))) # 305
paste0("npatientec=", length(unique(grep("xxx", metaclin2_slopes$identifiers, value = T)))) # 173
#n observation
paste0("nobservationitc=", length((grep("identifier", metaclin2_slopes$identifiers, value = T)))) # 408
paste0("nobservationec=", length((grep("xxx", metaclin2_slopes$identifiers, value = T)))) # 173
#n followup range
paste0("nfollowupRangeitc=[", max(table(grep("identifier", metaclin2_slopes$identifiers, value = T))),"-", # 3
       min(table(grep("identifier", metaclin2_slopes$identifiers, value = T))), "]") # 1

# | # Loading the models ----
# | # | # Proteome Model - Data preparation ----
# loading the model
load("~/Data/slope-models-adpkd/data/results/mod_pm.RData")

# extracting clinical annotation which has the features used in this model
validation_onlyprot = createData(data.slopesF, metaclin2_slopes, selected.model)
onlyprot = validation_onlyprot$gendat #581   7

# predicting slope
onlyprot$predSLOPEs = predict(validation_onlyprot$model$`LR model`, 
                              newdata = onlyprot)
# including important information
onlyprot$date = as.Date(validation_onlyprot$cannot[rownames(onlyprot),"merged_dates"], 
                        format = "%Y-%m-%d")
onlyprot$identifier = validation_onlyprot$cannot[rownames(onlyprot),"identifiers"]
onlyprot$ProjID = validation_onlyprot$cannot[rownames(onlyprot),"ProjID"]
onlyprot$eGFR = validation_onlyprot$cannot[rownames(onlyprot),"eGFR"] #581  12 #no NA

# | # | # Combined Model - Data preparation ----
# loading the model
load("~/Data/slope-models-adpkd/data/results/mod_cm.RData")

# extracting clinical annotation which has the features used in this model
validation_onlyprot_extra = createData(data.slopesF, metaclin2_slopes, 
                                       selected.model.ext, which = "extra")
onlyprotExtra = validation_onlyprot_extra$gendat # 581  11

# predicting slope
onlyprotExtra$predSLOPEsexted = predict(validation_onlyprot_extra$model, 
                                        newdata = onlyprotExtra)

# including important information
onlyprotExtra$date = as.Date(validation_onlyprot_extra$cannot[rownames(onlyprotExtra),"merged_dates"],
                             format = "%Y-%m-%d")
onlyprotExtra$identifier = validation_onlyprot_extra$cannot[rownames(onlyprotExtra),"identifiers"]
onlyprotExtra$ProjID = validation_onlyprot_extra$cannot[rownames(onlyprotExtra),"ProjID"]
onlyprotExtra$eGFR = validation_onlyprot_extra$cannot[rownames(onlyprotExtra),"eGFR"] 
# 581  15 #40NAs - 20 in mayo and 20 in predSLOPEsexted


# | # | # Combined Genotype Model - Data preparation ----
# loading the model
load("~/Data/slope-models-adpkd/data/results/mod_cgm.RData")

# extracting clinical annotation which has the features used in this model
validation_onlyprot_extraGEN = createData(data.slopesF, 
                                          subset(metaclin2_slopes, !grepl("vus|unknown",
                                                                          metaclin2_slopes$selectedmut, 
                                                                          ignore.case = T)), 
                                          selected.model.extGen, which = "extra")
onlyprotExtraGEN = validation_onlyprot_extraGEN$gendat

# predicting slope
onlyprotExtraGEN$predSLOPEsextedGEN = predict(validation_onlyprot_extraGEN$model,
                                              newdata = subset(onlyprotExtraGEN, !grepl("vus|unknown", onlyprotExtraGEN$selectedmut, ignore.case = T)))

# including important information
onlyprotExtraGEN$date = as.Date(validation_onlyprot_extraGEN$cannot[rownames(onlyprotExtraGEN),"merged_dates"], 
                                format = "%Y-%m-%d")
onlyprotExtraGEN$identifier = validation_onlyprot_extraGEN$cannot[rownames(onlyprotExtraGEN),"identifiers"]
onlyprotExtraGEN$ProjID = validation_onlyprot_extraGEN$cannot[rownames(onlyprotExtraGEN),"ProjID"]
onlyprotExtraGEN$eGFR = validation_onlyprot_extraGEN$cannot[rownames(onlyprotExtraGEN),"eGFR"] 
#343  16 #32 NAs - 10 in mayo and 6 in selectedmut and 16 in predSLOPEsextedGEN

# | # | # Proteome4 Model - Data preparation ----
# loading the model
load("~/Data/slope-models-adpkd/data/results/mod_p4m.RData")

# extracting clinical annotation which has the features used in this model
validation_onlyprot44 = createData(data.slopesF, metaclin2_slopes, 
                                   selected.model.min, which = "extra")
onlyprot44 = validation_onlyprot44$gendat

# predicting slope
onlyprot44$predSLOPEs44 = predict(validation_onlyprot44$model,
                                  newdata = onlyprot44)

# including important information
onlyprot44$date = as.Date(validation_onlyprot44$cannot[rownames(onlyprot44),"merged_dates"], 
                          format = "%Y-%m-%d")
onlyprot44$identifier = validation_onlyprot44$cannot[rownames(onlyprot44),"identifiers"]
onlyprot44$ProjID = validation_onlyprot44$cannot[rownames(onlyprot44),"ProjID"]
onlyprot44$eGFR = validation_onlyprot44$cannot[rownames(onlyprot44),"eGFR"]

# | # Loading Screening Cohort's Data ----
load("~/Data/slope-models-adpkd/data/results/data_screening.RData")

# | # Hemoglobin Plots ----
# after slope calculation
hbb_as = c(validation_onlyprot$proteome[,"HBB"], data.before.tolv["P68871",])
hbb_as = data.frame("rows" = names(hbb_as), "hbb" = as.numeric(hbb_as))
hbb_as$pjID = gsub("_.*", "", hbb_as$rows)
hbb_as[grep("batch", hbb_as$pjID), "pjID"] = "sc"

histhemo = ggplot(data = hbb_as, aes( x = hbb, fill = pjID, color = pjID)) + 
  geom_histogram(alpha = 0.5, position = "identity", binwidth = 0.5) + 
  theme_bw(base_size = 15) +
  scale_fill_manual(name = "Cohorts", 
                    labels = c("Screening (SC)",
                               "Internal/Temporal (ITC)",
                               "External (EC)"), 
                    values = c("#D81B60", "#FFB000", "#648FFF")) +
  scale_color_manual(name = "Cohorts", 
                     labels = c("Screening (SC)",
                                "Internal/Temporal (ITC)",
                                "External (EC)"),
                     values = c("#D81B60", "#FFB000", "#648FFF")) +
  xlab("") + ylab("Frequency")


denshemo = ggplot(data = hbb_as, aes(x = hbb, fill = pjID, color = pjID)) + 
  geom_density(alpha = 0.2, position = "identity") +
  theme_bw(base_size = 15) +
  scale_fill_manual(name = "Cohorts", 
                    labels = c("Screening (SC)",
                               "Internal/Temporal (ITC)",
                               "External (EC)"), 
                    values = c("#D81B60", "#FFB000", "#648FFF")) +
  scale_color_manual(name = "Cohorts", 
                     labels = c("Screening (SC)",
                                "Internal/Temporal (ITC)",
                                "External (EC)"),
                     values = c("#D81B60", "#FFB000", "#648FFF")) +
  xlab("Hemoglobin Subunit Beta (HBB) Intensity") + ylab("Density")

library(ggpubr)
ggarrange(histhemo, denshemo, common.legend = T, legend = "bottom", 
          nrow = 2, labels = c("A", "B"))
ggsave(paste0(file, "FigureS6.png"),width = 12, height = 10 , units = "in")

# comparing HBB distributions
# by checking the distribution of hbb in cohorts are normally distributed
hbb_as %>% group_by(pjID) %>% reframe(st = shapiro.test(hbb)$p.value)

# since it is not, we use non-parametric version of anova to test the difference
kruskal.test(hbb ~ pjID, data = hbb_as)
boxplot(hbb~pjID, hbb_as)

# | # Plot Observed-Predicted Slopes ----
# preparing screening cohort
ggplotSlopes$ProjID = "sc"
ggplotSlopes$predSLOPEs44 = selected.model.min.pred
ggplotSlopes$eGFR = cannot.before.tolv[rownames(ggplotSlopes),"merged_eGFR"]
ggplotSlopes$date = cannot.before.tolv[rownames(ggplotSlopes),"exam_date"] 
# 214  11

# calculating x and y axes limits for the scatter plots
xyrange_cmb = range(na.omit(c(onlyprot$before, # observed slopes
                              onlyprot$predSLOPEs, # prediction - Proteome Model
                              onlyprot44$predSLOPEs44, # prediction - Proteome4 Model
                              onlyprotExtra$predSLOPEsexted, # prediction - Combined Model
                              onlyprotExtraGEN$predSLOPEsextedGEN))) # prediction - Combined Genotype Model
xy.limits_cmb <- c(xyrange_cmb[1] - 0.2, xyrange_cmb[2] + 0.2) 

# | # | # Proteome Model ----
# adding octreotid usage as a column for indicating them as stars in the plot
onlyprot$octreotid = "no"
onlyprot$octreotid[match(intersect(c(octreotid, octreotidext), 
                                   onlyprot$identifier),onlyprot$identifier)] = "yes"

# combining clinical data of screening and validation cohort
cmbnData = rbind.data.frame(ggplotSlopes[,c("predSLOPEs", "before","ProjID",
                                            "identifier","eGFR", "date", "octreotid"), ],
                            onlyprot[,c("predSLOPEs", "before", "ProjID", "identifier",
                                        "eGFR", "date", "octreotid")]) # 795 7
# identification variable is generated
cmbnData$select1p1 = paste0(cmbnData$ProjID, "&", cmbnData$identifier, 
                            "&", cmbnData$date)

# generating "newpatients" variable indicates whether the patients included in 
# the proteome of Screening Cohort (iSCP) but were sampled at different time 
# points in ITC or were newly recruited (nSCP). 
cmbnData$newpatients = "iSCP"
cmbnData[grep(paste(grep("identifier", setdiff(unique(onlyprot$identifier), 
                                       unique(ggplotSlopes$identifier)),
                         value = T),collapse = "|"), 
              cmbnData$identifier),"newpatients"] = "nSCP" # 795   9

# sample sizes of cohorts 
# slope predictions were performed by Proteome Model
mytextobs = paste0("cmbnData-Sample Sizes\nsc:", 
                   length(which(cmbnData$ProjID == "sc")),
                   "\nitc:", length(which(cmbnData$ProjID == "itc")),
                   "\nec:", length(which(cmbnData$ProjID == "ec")))

# patient sizes of cohorts
# slope predictions were performed by Proteome Model
mytext = paste0("cmbnData-Patient Sizes\nsc:", 
                length(unique(cmbnData[cmbnData$ProjID == "sc", "identifier"])),
                "\nitc:", length(unique(cmbnData[cmbnData$ProjID == "itc","identifier"])),
                "\nec:", length(unique(cmbnData[cmbnData$ProjID == "ec","identifier"])))

# scatter plot
library(ggplot2)
onlyProtVal = ggplot(aes(x = before, y = predSLOPEs), data = cmbnData) +
  geom_point(aes(color = ProjID, shape = newpatients), alpha = 0.5, size = 3) +
  geom_point(data = subset(cmbnData, octreotid == "yes"), # octreotid info
             shape = 8, size = 4, aes(color = ProjID),show.legend = F) +
  scale_x_continuous(limits = xy.limits_cmb) + 
  scale_y_continuous(limits = xy.limits_cmb) +
  scale_color_manual(name = "Cohorts", 
                     labels = c("Screening (SC)", "Internal/Temporal (ITC)",
                                "External (EC)"),
                     values = c("#D81B60", "#FFB000", "#648FFF")) +
  scale_fill_manual(name = "Cohorts", labels = c("Screening (SC)", 
                                                 "Internal/Temporal (ITC)", 
                                                 "External (EC)"),
                    values = c("#D81B60", "#FFB000", "#648FFF")) +
  scale_shape_manual(name = "Patient\nCategory", labels = c("iSCP", "nSCP"), 
                     values = c(16, 15)) +
  xlab("Observed Slope") + ylab("Predicted Slope") +
  geom_abline(intercept = 0, slope = 1, linewidth = 1, linetype = "dashed", 
              colour = "gray50") +
  stat_smooth(method = "rlm", se = T, 
              aes(color = ProjID, fill = ProjID), alpha = 0.1) +
  theme_bw(base_size = 15) +
  guides(colour = guide_legend(order = 1),
         fill = guide_legend(order = 1),
         shape = guide_legend(order = 2))
onlyProtVal

# | # | # Combined Model ----
# adding octreotid usage as a column for indicating them as stars in the plot
onlyprotExtra$octreotid = "no"
onlyprotExtra$octreotid[match(intersect(c(octreotid, octreotidext), 
                                        onlyprotExtra$identifier),
                              onlyprotExtra$identifier)] = "yes"

# combining clinical data of screening and validation cohort
cmbnDataExtra = rbind.data.frame(ggplotSlopes[,c("predSLOPEsexted", "before",
                                                 "ProjID","identifier","eGFR", "date",
                                                 "octreotid"), ],
                                 onlyprotExtra[,c("predSLOPEsexted","before",
                                                  "ProjID","identifier","eGFR", "date",
                                                  "octreotid")]) # 795   7
# identification variable is generated
cmbnDataExtra$select1p1 = paste0(cmbnDataExtra$ProjID, "&", cmbnDataExtra$identifier, 
                                 "&", cmbnDataExtra$date)

# generating "newpatients" variable indicates whether the patients included in 
# the proteome of Screening Cohort (iSCP) but were sampled at different time 
# points in ITC or were newly recruited (nSCP). 
cmbnDataExtra$newpatients = "iSCP"
cmbnDataExtra[grep(paste(grep("identifier", setdiff(unique(onlyprotExtra$identifier), 
                                            unique(ggplotSlopes$identifier)), 
                              value = T), collapse = "|"), 
                   cmbnDataExtra$identifier),"newpatients"] = "nSCP"
# removing NAs
cmbnDataExtra = cmbnDataExtra[-which(is.na(cmbnDataExtra$predSLOPEsexted)),] 
# 773   9

# sample sizes of cohorts 
# slope predictions were performed by Combined Model
mytextextobs = paste0("cmbnDataExtra-Sample Sizes\nsc:", 
                      length(which(cmbnDataExtra$ProjID == "sc")),
                      "\nitc:", length(which(cmbnDataExtra$ProjID == "itc")),
                      "\nec:", length(which(cmbnDataExtra$ProjID == "ec")))

# patient sizes of cohorts 
# slope predictions were performed by Proteome Model
mytextext = paste0("cmbnDataExtra-Patient Sizes\nsc:",
                   length(unique(cmbnDataExtra[cmbnDataExtra$ProjID == "sc","identifier"])),
                   "\nitc:", length(unique(cmbnDataExtra[cmbnDataExtra$ProjID == "itc","identifier"])),
                   "\nec:", length(unique(cmbnDataExtra[cmbnDataExtra$ProjID == "ec","identifier"])))

# scatter plots
onlyProtValExt = ggplot(aes(x = before, y = predSLOPEsexted), 
                        data = cmbnDataExtra) +
  geom_point(aes(color = ProjID), alpha = 0.5, size = 3) +
  scale_x_continuous(limits = xy.limits_cmb) + 
  scale_y_continuous(limits = xy.limits_cmb) +
  scale_color_manual(name = "Cohorts", 
                     labels = c("Screening (SC)", "Internal/Temporal (ITC)",
                                "External (EC)"),
                     values = c("#D81B60", "#FFB000", "#648FFF")) +
  scale_fill_manual(name = "Cohorts", 
                    labels = c("Screening (SC)", "Internal/Temporal (ITC)",
                               "External (EC)"),
                    values = c("#D81B60", "#FFB000", "#648FFF")) +
  xlab("Observed Slope") + ylab("Predicted Slope") +
  geom_abline(intercept = 0, slope = 1, linewidth = 1, 
              linetype = "dashed", colour = "gray50") +
  stat_smooth(method = "rlm", se = T, 
              aes(color = ProjID, fill = ProjID), alpha = 0.1) +
  theme_bw(base_size = 15)
onlyProtValExt
ggsave(paste0(file, "Figure3.pdf"),width = 12, height = 10 , units = "in")

onlyProtValExtstr = ggplot(aes(x = before, y = predSLOPEsexted), 
                           data = cmbnDataExtra) +
  geom_point(aes(color = ProjID,shape = newpatients), alpha = 0.5, size = 3) +
  geom_point(data = subset(cmbnDataExtra, octreotid == "yes"), shape = 8, 
             size = 4, aes(color = ProjID), show.legend = F) +
  scale_x_continuous(limits = xy.limits_cmb) + 
  scale_y_continuous(limits = xy.limits_cmb) +
  scale_color_manual(name = "Cohorts", 
                     labels = c("Screening (SC)", "Internal/Temporal (ITC)",
                                "External (EC)"),
                     values = c("#D81B60", "#FFB000", "#648FFF")) +
  scale_fill_manual(name = "Cohorts", 
                    labels = c("Screening (SC)", "Internal/Temporal (ITC)", 
                               "External (EC)"),
                    values = c("#D81B60", "#FFB000", "#648FFF")) +
  scale_shape_manual(name = "Patient\nCategory", labels = c("iSCP","nSCP"), 
                     values = c(16, 15)) +
  xlab("Observed Slope") + ylab("Predicted Slope") +
  geom_abline(intercept = 0, slope = 1, linewidth = 1, linetype = "dashed", 
              colour = "gray50") +
  stat_smooth(method = "rlm", se = T, aes(color = ProjID, fill = ProjID),
              alpha = 0.1) + theme_bw(base_size = 15) +
  guides(colour = guide_legend(order = 1),
         fill = guide_legend(order = 1),
         shape = guide_legend(order = 2))
onlyProtValExtstr

# | # | # Combined Genotype Model ----
# adding octreotid usage as a column for putting stars in the plot
onlyprotExtraGEN$octreotid = "no"
onlyprotExtraGEN$octreotid[match(intersect(c(octreotid, octreotidext), 
                                           onlyprotExtraGEN$identifier),
                                 onlyprotExtraGEN$identifier)] = "yes"

# combining clinical data of screening and validation cohort
cmbnDataGEN = rbind.data.frame(ggplotSlopes[,c("predSLOPEsextedGEN", "before", 
                                               "ProjID","identifier","eGFR","date",
                                               "octreotid"), ], # 557   7
                               onlyprotExtraGEN[,c("predSLOPEsextedGEN", 
                                                   "before", "ProjID", "identifier", 
                                                   "eGFR", "date", "octreotid")]) 

# identification variable is generated
cmbnDataGEN$select1p1 = paste0(cmbnDataGEN$ProjID, "&", cmbnDataGEN$identifier, 
                               "&", cmbnDataGEN$date)

# generating "newpatients" variable indicates whether the patients included in 
# the proteome of Screening Cohort (iSCP) but were sampled at different time 
# points in ITC or were newly recruited (nSCP). 
cmbnDataGEN$newpatients = "iSCP"
cmbnDataGEN[grep(paste(grep("identifier", setdiff(unique(onlyprot$identifier), 
                                          unique(ggplotSlopes$identifier)),
                            value = T),collapse = "|"), 
                 cmbnDataGEN$identifier),"newpatients"] = "nSCP"

# removing NAs
cmbnDataGEN = cmbnDataGEN[-which(is.na(cmbnDataGEN$predSLOPEsextedGEN)),] #441 9

# sample sizes of cohorts 
# slope predictions were performed by Combined Genotype Model
mytextgenobs = paste0("cmbnDataGEN-Sample Sizes\nsc:",
                      length(which(cmbnDataGEN$ProjID == "sc")),
                      "\nitc:", length(which(cmbnDataGEN$ProjID == "itc")),
                      "\nec:", length(which(cmbnDataGEN$ProjID == "ec")))

# patient sizes of cohorts 
# slope predictions were performed by Combined Genotype Model
mytextgen = paste0("cmbnDataGEN-Patient Sizes\nsc:",
                   length(unique(cmbnDataGEN[cmbnDataGEN$ProjID == "sc","identifier"])),
                   "\nitc:", length(unique(cmbnDataGEN[cmbnDataGEN$ProjID == "itc","identifier"])),
                   "\nec:", length(unique(cmbnDataGEN[cmbnDataGEN$ProjID == "ec","identifier"])))

# scatter plot
library(ggplot2)
onlyProtValExtGEN = ggplot(aes(x = before, y = predSLOPEsextedGEN), 
                           data = cmbnDataGEN) +
  geom_point(aes(color = ProjID, shape = newpatients), alpha = 0.5, size = 3) +
  geom_point(data = subset(cmbnDataGEN, octreotid == "yes"), shape = 8, 
             size = 4, aes(color = ProjID), show.legend = F) +
  scale_x_continuous(limits = xy.limits_cmb) + 
  scale_y_continuous(limits = xy.limits_cmb) +
  scale_shape_manual(name = "Patient\nCategory", labels = c("iSCP","nSCP"), 
                     values = c(16, 15)) +
  scale_color_manual(name = "Cohorts", labels = c("Screening (SC)", 
                                                  "Internal/Temporal (ITC)",
                                                  "External (EC)"),
                       values = c("#D81B60", "#FFB000", "#648FFF")) +
  scale_fill_manual(name = "Cohorts", labels = c("Screening (SC)",
                                                 "Internal/Temporal (ITC)",
                                                 "External (EC)"),
                    values = c("#D81B60", "#FFB000", "#648FFF")) +
  xlab("Observed Slope") + ylab("Predicted Slope") +
  geom_abline(intercept = 0, slope = 1, linewidth = 1, 
              linetype = "dashed", colour = "gray50") +
  stat_smooth(method = "rlm", se = T, 
              aes(color = ProjID, fill = ProjID), alpha = 0.1) +
  theme_bw(base_size = 15) + 
  guides(colour = guide_legend(order = 1),
         fill = guide_legend(order = 1),
         shape = guide_legend(order = 2))
onlyProtValExtGEN

# | # | # Proteome4 Model ----
# adding octreotid usage as a column for putting stars in the plot
onlyprot44$octreotid = "no"
onlyprot44$octreotid[match(intersect(c(octreotid, octreotidext), 
                                     onlyprot44$identifier),onlyprot44$identifier)] = "yes"

# combining clinical data of screening and validation cohort
cmbnData44 = rbind.data.frame(ggplotSlopes[,c("predSLOPEs44", "before", 
                                              "ProjID", "identifier", "eGFR", 
                                              "date","octreotid")], # 795   7
                              onlyprot44[,c("predSLOPEs44", "before", "ProjID",
                                            "identifier", "eGFR", "date", "octreotid")])

# identification variable is generated
cmbnData44$select1p1 = paste0(cmbnData44$ProjID, "&", cmbnData44$identifier, 
                               "&", cmbnData44$date)

# generating "newpatients" variable indicates whether the patients included in 
# the proteome of Screening Cohort (iSCP) but were sampled at different time 
# points in ITC or were newly recruited (nSCP). 
cmbnData44$newpatients = "iSCP"
cmbnData44[grep(paste(grep("identifier", setdiff(unique(onlyprot44$identifier), 
                                         unique(ggplotSlopes$identifier)),value = T),
                      collapse = "|"), cmbnData44$identifier),"newpatients"] = "nSCP"

# sample sizes of cohorts 
# slope predictions were performed by Proteome4 Model
mytext44obs = paste0("cmbnData44-Sample Sizes\nsc:",
                     length(which(cmbnData44$ProjID == "sc")),
                     "\nitc:", length(which(cmbnData44$ProjID == "itc")),
                     "\nec:", length(which(cmbnData44$ProjID == "ec")))

# patient sizes of cohorts 
# slope predictions were performed by Proteome4 Model
mytext44 = paste0("cmbnData44-Patient Sizes\nsc:",
                  length(unique(cmbnData44[cmbnData44$ProjID == "sc","identifier"])),
                  "\nitc:", length(unique(cmbnData44[cmbnData44$ProjID == "itc","identifier"])),
                  "\nec:", length(unique(cmbnData44[cmbnData44$ProjID == "ec","identifier"])))

# scatter plot
library(ggplot2)
onlyProtVal44 = ggplot(aes(x = before, y = predSLOPEs44), data = cmbnData44) +
  geom_point(aes(color = ProjID, shape = newpatients), alpha = 0.5, size = 3) +
  geom_point(data = subset(cmbnData44, octreotid == "yes"), 
             shape = 8, size = 4, aes(color = ProjID), show.legend = F) +
  scale_x_continuous(limits = xy.limits_cmb) + 
  scale_y_continuous(limits = xy.limits_cmb) +
  scale_color_manual(name = "Cohorts", labels = c("Screening (SC)", 
                                                  "Internal/Temporal (ITC)",
                                                  "External (EC)"),
                     values = c("#D81B60", "#FFB000", "#648FFF")) +
  scale_fill_manual(name = "Cohorts", labels = c("Screening (SC)",
                                                 "Internal/Temporal (ITC)",
                                                 "External (EC)"),
                    values = c("#D81B60", "#FFB000", "#648FFF")) +
  scale_shape_manual(name = "Patient\nCategory", labels = c("iSCP","nSCP"), 
                     values = c(16, 15)) +
  xlab("Observed Slope") + ylab("Predicted Slope") +
  geom_abline(intercept = 0, slope = 1, linewidth = 1, 
              linetype = "dashed", colour = "gray50") +
  stat_smooth(method = "rlm", se = T, aes(color = ProjID, 
                                          fill = ProjID), alpha = 0.1) +
  theme_bw(base_size = 15) +
  guides(colour = guide_legend(order = 1),
         fill = guide_legend(order = 1),
         shape = guide_legend(order = 2))
onlyProtVal44

# TableS8
samplesizesValidation = paste0(mytext, "\n\n", mytextext, "\n\n", mytextgen, 
                               "\n\n", mytext44)
writeLines(samplesizesValidation, paste0(file,"patientsizeValidation", 
                                         gsub("-", "", Sys.Date()), ".txt"))

patientsizesValidation = paste0(mytextobs, "\n\n", mytextextobs, "\n\n", 
                                mytextgenobs, "\n\n", mytext44obs)
writeLines(patientsizesValidation, paste0(file,"sampsizeValidation", 
                                          gsub("-", "", Sys.Date()), ".txt"))

# | # Plot Same Patients Slope Prediction at Different Time Points ----
# | # | # Proteome Model ----
cmbnData.by <- by(cmbnData,cmbnData$identifier,function(x) x)
cmbnData.by.s <- cmbnData.by[names(which(table(cmbnData$identifier) > 1))]

mapx <- do.call(rbind,lapply(cmbnData.by.s,function(x) {
  x <- x[order(x$date,decreasing = F),]
  diffd <- as.numeric(x$date[-1] - x$date[1])/365.25
  dfx <- data.frame(id = x$identifier[1], baselinePred = x$predSLOPEs[1], 
                    followPred = x$predSLOPEs[-1],
                    diffDate = diffd, before = x$before[1])
  return(dfx)
}))

cmbnData.byp <- ggplot(mapx, aes(x = baselinePred, y = followPred, 
                                 colour = factor(round(diffDate)))) +
  geom_point(size = 3, alpha = 0.8) +
  geom_line(aes(group = id), linewidth = 1.7, alpha = 0.8) + 
  geom_abline(slope = c(1,1,1), intercept = c(-1,0,1), lty = 2) + 
  theme_bw(base_size = 15) +
  scale_color_manual(name = "Time Difference (in Years)", 
                     values = c("#290AD8", "#3FA0FF", "#FFAD50", 
                                "#F76D70", "#D82632"), labels = c(0,1,2,3,4)) +
  xlab("Slope Prediction at Baseline") + ylab("Slope Prediction at Follow-ups")

# | # | # Combined Model ----
cmbnDataExtra.by <- by(cmbnDataExtra,cmbnDataExtra$identifier,function(x) x)
cmbnDataExtra.by.s <- cmbnDataExtra.by[names(which(table(cmbnDataExtra$identifier) > 1))]

mapx_ex <- do.call(rbind,lapply(cmbnDataExtra.by.s,function(x) {
  x <- x[order(x$date,decreasing = F),]
  diffd <- as.numeric(x$date[-1] - x$date[1])/365.25
  dfx <- data.frame(id = x$identifier[1], baselinePred = x$predSLOPEsexted[1],
                    followPred = x$predSLOPEsexted[-1], diffDate = diffd, 
                    before = x$before[1])
  return(dfx)
}))


cmbnDataExtra.byp <- ggplot(mapx_ex,aes(x = baselinePred, y = followPred,
                                        colour = factor(round(diffDate)))) +
  geom_point(size = 3, alpha = 0.8) +
  geom_line(aes(group = id), linewidth = 1.7, alpha = 0.8) +
  geom_abline(slope = c(1,1,1), intercept = c(-1,0,1), lty = 2) +
  theme_bw(base_size = 15) +
  scale_color_manual(name = "Time Difference (in Years)", 
                     values = c("#290AD8", "#3FA0FF", "#FFAD50", 
                                "#F76D70", "#D82632"), labels = c(0,1,2,3,4)) +
  xlab("Slope Prediction at Baseline") + ylab("Slope Prediction at Follow-ups")


# | # | # Combined Genotype Model ----
cmbnDataGEN.by <- by(cmbnDataGEN,cmbnDataGEN$identifier,function(x) x)
cmbnDataGEN.by.s <- cmbnDataGEN.by[names(which(table(cmbnDataGEN$identifier) > 1))]

mapx_gen <- do.call(rbind,lapply(cmbnDataGEN.by.s,function(x) {
  x <- x[order(x$date,decreasing = F),]
  diffd <- as.numeric(x$date[-1] - x$date[1])/365.25
  dfx <- data.frame(id = x$identifier[1], baselinePred = x$predSLOPEsextedGEN[1],
                    followPred = x$predSLOPEsextedGEN[-1],
                    diffDate = diffd,before = x$before[1])
  return(dfx)
}))

cmbnDataGEN.byp <- ggplot(mapx_gen,aes(x = baselinePred, y = followPred,
                                       colour = factor(round(diffDate)))) +
  geom_point(size = 3, alpha = 0.8) +
  geom_line(aes(group = id), linewidth = 1.7, alpha = 0.8) +
  geom_abline(slope = c(1,1,1),intercept = c(-1,0,1),lty = 2) +
  theme_bw(base_size = 15) +
  scale_color_manual(name = "Time Difference (in Years)", 
                     values = c("#290AD8", "#3FA0FF", "#FFAD50", 
                                "#F76D70", "#D82632"), labels = c(0,1,2,3,4)) +
  xlab("Slope Prediction at Baseline") + ylab("Slope Prediction at Follow-ups")

# | # | # Proteome4 Model ----
cmbnData44.by <- by(cmbnData44,cmbnData44$identifier,function(x) x)
cmbnData44.by.s <- cmbnData44.by[names(which(table(cmbnData44$identifier) > 1))]

mapx44 <- do.call(rbind,lapply(cmbnData44.by.s,function(x) {
  x <- x[order(x$date,decreasing = F),]
  diffd <- as.numeric(x$date[-1] - x$date[1])/365.25
  dfx <- data.frame(id = x$identifier[1], baselinePred = x$predSLOPEs[1],
                    followPred = x$predSLOPEs[-1],
                    diffDate = diffd, before = x$before[1])
  return(dfx)
}))

cmbnData44.byp <- ggplot(mapx44,aes(x = baselinePred, y = followPred,
                                    colour = factor(round(diffDate)))) +
  geom_point(size = 3, alpha = 0.8) + 
  geom_line(aes(group = id), linewidth = 1.7, alpha = 0.8) + 
  geom_abline(slope = c(1,1,1), intercept = c(-1,0,1), lty = 2) + 
  theme_bw(base_size = 15) +
  scale_color_manual(name = "Time Difference\n(in Years)", 
                     values = c("#290AD8", "#3FA0FF", "#FFAD50", 
                                "#F76D70", "#D82632"), labels = c(0,1,2,3,4)) +
  xlab("Slope Prediction at Baseline") + ylab("Slope Prediction at Follow-ups")

save.image(paste0(file,"data_Validation.RData"))
writeLines(capture.output(sessionInfo()),
           paste0(file,"sessionInfo03valcomb",
                  gsub("-", "", Sys.Date()),".txt"))

library(ggpubr)
library(MASS)
require(grid)
leftsd = ggarrange(onlyProtVal,  onlyProtValExtstr, onlyProtValExtGEN,
                   onlyProtVal44, nrow = 4, common.legend = T, legend = "bottom",
                   labels = c("A", "B", "C", "D"))

rightsd = ggarrange(cmbnData.byp , cmbnDataExtra.byp, 
                    cmbnDataGEN.byp, cmbnData44.byp,
                    nrow = 4, common.legend = T, legend = "bottom",
                    labels = c("E", "F", "G", "H"))
ggarrange(leftsd, rightsd, widths = c(1,1.3))
ggsave(paste0(file, "FigureS7.png"),width = 22 , height = 33 , units = "in")

# | # Future eGFRs ----
# | # | # Internal/Temporal Cohort ----
# preparing the data
cmbnData[intersect(rownames(cmbnData), rownames(metaclin2_slopes)), "age"] = metaclin2_slopes[intersect(rownames(cmbnData), rownames(metaclin2_slopes)), "age"]
cmbnData[intersect(rownames(cmbnData), rownames(metaclin2_slopes)), "gender"] = metaclin2_slopes[intersect(rownames(cmbnData), rownames(metaclin2_slopes)), "gender"]
cmbnData[intersect(rownames(cmbnData), rownames(metaclin2_slopes)), "mayo5"] = metaclin2_slopes[intersect(rownames(cmbnData), rownames(metaclin2_slopes)), "mayo5"]

selectiondf = data.frame("rn" = grep("^itc_", rownames(cmbnData), value = T),
                         "identifiers" = cmbnData$identifier[grep("^itc_", rownames(cmbnData))],
                         "age" = cmbnData$age[grep("^itc_", rownames(cmbnData))])

# selecting the first measurement for the scatter plot
library(tidyverse)
min_selectiondf = selectiondf %>% group_by(identifiers) %>% slice_min(age)
min_selectiondf = min_selectiondf[-which(duplicated(min_selectiondf[,c("age","identifiers")])),]
cmbnData22222 = cmbnData[min_selectiondf$rn,]

#for calculating future eGFR --> you need to get positive time difference
calctimedif_itc = sapply(min_selectiondf$identifiers, function(x){
  calcTimeDiff(x, merged.creas.final3, cmbnData22222,
               "patient_id", "identifier", "exam_date", 
               "date",eGFR = "patient__eGFR")
},simplify = FALSE, USE.NAMES = TRUE)

# to select indexes after proteome sampling
indxx_itc = sapply(sapply(calctimedif_itc, "[[", "timedif"), function(x) x > 0)

# to select time differences after proteome sampling
timediffall_itc = sapply(sapply(calctimedif_itc, "[[", "timedif"), 
                         function(x) x[which(x > 0)]) 

#to select eGFRs after proteome sampling
eGFRsclc_itc = mapply(function(x,y){
  x[y]},x = sapply(calctimedif_itc, "[[", "eGFRs"), y = indxx_itc)

# Future eGFR calculation with Proteome Model
onlyProtFuture_itc = mapply(function(x,y){
  calcFuteGFR(timedif = x, id = y, idcol = "identifier" , data = cmbnData22222, 
              which = "prot slopes", predtype = "predSLOPEs")
}, x = timediffall_itc, y = names(timediffall_itc))

# Future eGFR calculation with Combined Model
cmbnDataExtra22222 = cmbnDataExtra[intersect(min_selectiondf$rn, 
                                             rownames(cmbnDataExtra)),]
extendedProtFuture_itc = mapply(function(x,y){
  calcFuteGFR(timedif = x, id = y, idcol = "identifier" , data = cmbnDataExtra22222, 
              which = "prot slopes", predtype = "predSLOPEsexted")
}, x = timediffall_itc, y = names(timediffall_itc))

# Future eGFR calculation with MIC Model
mayoFuture_itc = mapply(function(x,y){
  calcFuteGFR(timedif = x, id = y, idcol = "identifier" , 
              data = cmbnData22222, which = "mayo slopes")
}, x = timediffall_itc, y = names(timediffall_itc))

remthese_itc = names(which(sapply(mayoFuture_itc, is.null) | unlist(lapply(sapply(mayoFuture_itc, is.na), any))))
mayoFuture_itc[remthese_itc] <- NA
extendedProtFuture_itc[remthese_itc] <- NA 

scatterFuture_itc = do.call(rbind.data.frame, lapply(names(timediffall_itc), function(x){
  # if a patient & their time difference data & predicted future eGFR from
  # predicted slopes (Proteome Model and Combined Model) & actual eGFR 
  # measurement exist 
  if (!is.null(x) & !is.null(as.numeric(timediffall_itc[[x]])) & !is.null(onlyProtFuture_itc[[x]]) & !is.null(extendedProtFuture_itc[[x]]) & !is.null(eGFRsclc_itc[[x]])) {
    # generate a vector for that patient containing all these variables
    # and future eGFR from mayo slope
    return(cbind.data.frame(id = x, 
                            timeidffyear = as.numeric(timediffall_itc[[x]]), 
                            onlyProtFuteGFR = onlyProtFuture_itc[[x]], 
                            extFuteGFR = extendedProtFuture_itc[[x]], 
                            realFuteGFR = eGFRsclc_itc[[x]], 
                            mayoFut = mayoFuture_itc[[x]]))
  }
})) # and then combine them rowwise

scatterFuture_itcMaxd = scatterFuture_itc %>% group_by(id) %>% 
  filter(timeidffyear == max(timeidffyear)) # selecting the maximum timediffyear
scatterFuture_itcMaxd2 = scatterFuture_itcMaxd
scatterFuture_itcMaxd2$deltaSonlyP = scatterFuture_itcMaxd2$realFuteGFR - scatterFuture_itcMaxd2$onlyProtFuteGFR # delta eGFR Proteome Model
scatterFuture_itcMaxd2$deltaSPext = scatterFuture_itcMaxd2$realFuteGFR - scatterFuture_itcMaxd2$extFuteGFR # delta eGFR Combined Model
scatterFuture_itcMaxd2$deltaSPmayo = scatterFuture_itcMaxd2$realFuteGFR - scatterFuture_itcMaxd2$mayoFut # delta eGFR MIC Model

library(reshape2)
scatterFuture_itcMaxd2[,c("extFuteGFR", "onlyProtFuteGFR", 
                          "realFuteGFR", "mayoFut")] = NULL

# melting the data for ggplot
scatterFuture_itcMaxd2 = melt(scatterFuture_itcMaxd2, 
                              id.vars = c("timeidffyear", "id"))
# scatter plot
scatterOnlyProtFig_itc = ggplot(scatterFuture_itcMaxd2, 
                                aes(x = timeidffyear, y = value)) + 
  geom_point(aes(color = variable), size = 3, alpha = 0.7) + 
  scale_y_continuous(limits = ylimits.ed, breaks = seq(-50, 50, by = 10)) + 
  theme_bw(base_size = 17) + scale_x_continuous(breaks = seq(0, 8, 1)) + 
  geom_hline(yintercept = 0, linewidth = 2, colour = "darkred", linetype = 3) +
  scale_colour_manual(name = "Models", 
                      values = c("#3C97DA", "#165079", "#FE9D52"), 
                      labels = c("Proteome Model", "Combined Model",
                                 "MIC Model")) +
  scale_fill_manual(name = "Models", values = c("#3C97DA", "#165079", "#FE9D52"), 
                    labels = c("Proteome Model", "Combined Model", "MIC Model")) + 
  xlab(expression("Latest eGFR for each patient (in years)")) + 
  ylab(expression(~Delta~"eGFR"["Observed-Predicted"])) +
  stat_smooth(method = "lm", se = T, aes(fill = variable, colour = variable), 
              alpha = 0.3, level = 0.95) +
  stat_smooth(method = "lm", se = F, aes(fill = variable, colour = variable)) +
  theme(legend.position = 'bottom', 
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(size = 16))
scatterOnlyProtFig_itc

# | # | # External Cohort ----
#for calculating future eGFR --> you need to get positive time difference
calctimedif_ec = sapply(grep("xxx", cmbnData$identifier, value = T), function(x){
  calcTimeDiff(x, cannot.long, cmbnData, 
               "protmatch", "identifier", "exam_date", 
               "date",eGFR = "patient__eGFR")
}, simplify = FALSE, USE.NAMES = TRUE)

# to select indexes after proteome sampling
indxx_ec = sapply(sapply(calctimedif_ec, "[[", "timedif"), function(x) x > 0)

# to select time differences after proteome sampling
timediffall_ec = sapply(sapply(calctimedif_ec, "[[", "timedif"), 
                        function(x) x[which(x > 0)])

#to select eGFRs after proteome sampling
eGFRsclc_ec = mapply(function(x,y){ 
  x[y]}, x = sapply(calctimedif_ec, "[[", "eGFRs"), y = indxx_ec)

# Future eGFR calculation with Proteome Model
onlyProtFuture_ec = mapply(function(x,y){
  calcFuteGFR(timedif = x, id = y, idcol = "identifier" , data = cmbnData, 
              which = "prot slopes", predtype = "predSLOPEs")
}, x = timediffall_ec, y = names(timediffall_ec))

# Future eGFR calculation with Combined Model
extendedProtFuture_ec = mapply(function(x,y){
  calcFuteGFR(timedif = x, id = y, idcol = "identifier" , data = cmbnDataExtra, 
              which = "prot slopes", predtype = "predSLOPEsexted")
}, x = timediffall_ec, y = names(timediffall_ec))

# Future eGFR calculation with MIC Model
mayoFuture_ec = mapply(function(x,y){
  calcFuteGFR(timedif = x, id = y, idcol = "identifier" , 
              data = cmbnData, which = "mayo slopes")
}, x = timediffall_ec, y = names(timediffall_ec))

remthese_ec = names(which(sapply(mayoFuture_ec, is.null) | unlist(lapply(sapply(mayoFuture_ec, is.na), any))))
mayoFuture_ec[remthese_ec] <- NA
extendedProtFuture_ec[remthese_ec] <- NA 

# if a patient & their time difference data & predicted future eGFR from
# predicted slopes (Proteome Model and Combined Model) & actual eGFR 
# measurement exist 
scatterFuture_ec = do.call(rbind.data.frame, lapply(names(timediffall_ec), function(x){
  if (!is.null(x) & !is.null(as.numeric(timediffall_ec[[x]])) & !is.null(onlyProtFuture_ec[[x]]) & !is.null(extendedProtFuture_ec[[x]]) & !is.null(eGFRsclc_ec[[x]])) {
    # generate a vector for that patient containing all these variables
    # and future eGFR from mayo slope
    return(cbind.data.frame(id = x, 
                            timeidffyear = as.numeric(timediffall_ec[[x]]), 
                            onlyProtFuteGFR = onlyProtFuture_ec[[x]],
                            extFuteGFR = extendedProtFuture_ec[[x]], 
                            realFuteGFR = eGFRsclc_ec[[x]], 
                            mayoFut = mayoFuture_ec[[x]]))
  }
}))# and then combine them rowwise

scatterFuture_ecMaxd = scatterFuture_ec %>% group_by(id) %>%
  filter(timeidffyear == max(timeidffyear)) # selecting the maximum timediffyear
scatterFuture_ecMaxd2 = scatterFuture_ecMaxd[-which(duplicated(scatterFuture_ecMaxd)),]
scatterFuture_ecMaxd2$deltaSonlyP = scatterFuture_ecMaxd2$realFuteGFR - scatterFuture_ecMaxd2$onlyProtFuteGFR # delta eGFR Proteome Model
scatterFuture_ecMaxd2$deltaSPext = scatterFuture_ecMaxd2$realFuteGFR - scatterFuture_ecMaxd2$extFuteGFR # delta eGFR Combined Model
scatterFuture_ecMaxd2$deltaSPmayo = scatterFuture_ecMaxd2$realFuteGFR - scatterFuture_ecMaxd2$mayoFut # delta eGFR MIC Model

scatterFuture_ecMaxd2[,c("extFuteGFR", "onlyProtFuteGFR",
                         "realFuteGFR", "mayoFut")] = NULL

# melting the data for ggplot
scatterFuture_ecMaxd2 = melt(scatterFuture_ecMaxd2, 
                             id.vars = c("timeidffyear","id"))
# scatter plot
scatterOnlyProtFig_ec = ggplot(scatterFuture_ecMaxd2, aes(x = timeidffyear, y = value)) +
  geom_point(aes(color = variable), size = 3, alpha = 0.7) + 
  scale_y_continuous(limits = ylimits.ed, breaks = seq(-50, 50, by = 10)) + 
  theme_bw(base_size = 17) + scale_x_continuous(breaks = seq(0, 8, 1)) + 
  geom_hline(yintercept = 0, linewidth = 2, colour = "darkred", linetype = 3) +
  scale_colour_manual(name = "Models", 
                      values = c("#3C97DA", "#165079", "#FE9D52"), 
                      labels = c("Proteome Model", "Combined Model", 
                                 "MIC Model")) +
  scale_fill_manual(name = "Models", 
                    values = c("#3C97DA", "#165079", "#FE9D52"), 
                    labels = c("Proteome Model", "Combined Model", 
                               "MIC Model")) + 
  xlab(expression("Latest eGFR for each patient (in years)")) + 
  ylab(expression(~Delta~"eGFR"["Observed-Predicted"])) +
  stat_smooth(method = "lm", se = T, aes(fill = variable, colour = variable), 
              alpha = 0.3, level = 0.95) +
  stat_smooth(method = "lm", se = F, aes(fill = variable, colour = variable)) +
  theme(legend.position = 'bottom', 
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(size = 16))
scatterOnlyProtFig_ec

# Plotting both 
library(ggpubr)
ggarrange(scatterOnlyProtFig_itc, scatterOnlyProtFig_ec, common.legend = T, 
          legend = "bottom", ncol = 2, labels = c("A", "B"))
ggsave(paste0(file, "FigureS8.png"),width = 15, height = 10 , units = "in")

# | # Table1 and TableS2 preparation ----
# Table 1 - cont'd
metaclin2_slopes$CKD = eGFRtoStageAB(metaclin2_slopes$eGFR)
library(tibble)
firstvalidation = metaclin2_slopes %>% group_by(ProjID) %>% 
  reframe(n = n(),
          sexpct = length(which(gender == "female"))*100/length(na.omit(gender)),
          age = paste0(round(median(age),1)," [",
                       round(IQR(age),1), "]"),
          egfr = paste0(round(median(eGFR),1)," [",
                        round(IQR(eGFR),1), "]"),
          tkv = paste0(round(median(tkvmerged,na.rm = T), 1)," [",
                       round(IQR(tkvmerged, na.rm = T),1), "]"))  %>%
  t() %>% as.data.frame() %>% rownames_to_column(var = "RowNames") %>% unname

mayotb = metaclin2_slopes %>% group_by(ProjID, mayo5) %>% 
  reframe(n = n()) %>% as.data.frame() %>% unname()
mayotb2 = do.call(cbind.data.frame, split(mayotb, f = mayotb[,1]))[,c(2,3,6)]

ckdtb = metaclin2_slopes %>% group_by(ProjID, CKD) %>% 
  reframe(n = n()) %>% as.data.frame() %>% unname()
ckdtb2 = split(ckdtb, f = ckdtb[,1])
ckdtb2$`ec` = rbind(ckdtb2$`ec`, c(ec, 5, 0))
ckdtb2 = do.call(cbind.data.frame,ckdtb2)[,c(2,3,6)]

colnames(mayotb2) = c("rowns", "vars", NA)
colnames(ckdtb2) = c("rowns", "vars", NA)

# Table S2 - cont'd
firstvalidation2 = metaclin2_slopes %>%  group_by(ProjID) %>% 
  reframe(hypt_n = length(na.omit(hypertension)),
          hypt_none = length(na.omit(which(hypertension == "0"))),
          hypt_yes = length(na.omit(grep("1|2|3",hypertension))),
          mut_n = length(na.omit(selectedmut)),
          mut_pkd1t = length(na.omit(which(selectedmut == "PKD1T"))),
          mut_pkd1nt = length(na.omit(which(selectedmut == "PKD1NT"))),
          mut_pkd2 = length(na.omit(which(selectedmut == "PKD2")))) %>% 
  t() %>% as.data.frame() %>% 
  rownames_to_column(var = "RowNames") %>% unname 

table1validationch = do.call(rbind.data.frame,
                             lapply(list(firstvalidation,mayotb2,
                                         ckdtb2,firstvalidation2),function(x){
                                           x[,2] = as.character(x[,2])
                                           colnames(x) = c("rowns","vars")
                                           return(x) 
                                         }))

table1subsetCoh = rbind.data.frame(table1subsetCoh[c(1:which(table1subsetCoh$rowns == 4)),], c(5, 0), 
                                   table1subsetCoh[c((which(table1subsetCoh$rowns == 4) + 1):nrow(table1subsetCoh)),])
colnames(table1validationch) = c("rowns","itc", "ec")
table1validationch = table1validationch[-grep("ProjID", table1validationch$rowns),]
table1complete = cbind.data.frame(table1subsetCoh[-grep("us_", table1subsetCoh$rowns),], 
                                  table1validationch)
table1complete = table1complete[,c(1,2,4,5)]
colnames(table1complete) = c("rowns", "sc", "itc", "ec")
library(xlsx)
write.xlsx(table1complete,
           file = paste0(file,"Table1.xlsx"),
           sheetName = "Table1")

supplmtitc = metaclin2_slopes %>% filter(ProjID == "itc") %>% 
  reframe(n = n(),
          posfam = length(which(postfam == "true"))*100/length((postfam)),
          hypt_n = length(na.omit(hypertension)),
          hypt_none = length(na.omit(which(hypertension == "0"))),
          hypt_l35 = length(na.omit(which(hypertension == "1"))),
          hypt_g35 = length(na.omit(which(hypertension == "2"))),
          hypt_ukage = length(na.omit(which(hypertension == "3"))),
          us_n = length(na.omit(usyeni)),
          us_none = length(na.omit(which(usyeni == "0"))),
          us_l35 = length(na.omit(which(usyeni == "1"))),
          us_g35 = length(na.omit(which(usyeni == "2"))),
          us_uk = length(na.omit(which(usyeni == "-1")))) %>% 
  t() %>% as.data.frame() %>% rownames_to_column(var = "RowNames") %>% unname 

supplCohCharact = cbind.data.frame(subsebtADPKDsupplmt, supplmtitc)
write.xlsx(supplCohCharact,
           file = paste0(file,"Table1.xlsx"),
           sheetName = "TableS2", append = TRUE)
