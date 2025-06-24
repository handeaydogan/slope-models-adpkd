# SLOPE PREDICTIONS ----
source("~/Data/slope-models-adpkd/code/00functions.R", echo = TRUE)
source("~/Data/slope-models-adpkd/code/SC/SC_01_annot.half.R", echo = TRUE)
source("~/Data/slope-models-adpkd/code/SC/SC_02_interventions.R", echo = TRUE)
source("~/Data/slope-models-adpkd/code/SC/SC_03_mergedCreasUpd.R", echo = TRUE)
source("~/Data/slope-models-adpkd/code/SC/SC_04_annots_merging.R", echo = TRUE)
source("~/Data/slope-models-adpkd/code/SC/SC_05_prots.R", echo = TRUE)

# | # Preparing the Data ----
# selecting the ones that have before slopes
cannot.before.tolv = cannot.slopesF[which(cannot.slopesF$bTF == T & !is.na(cannot.slopesF$before)),] # 234 306

# removing the ones that have slope bigger than 5 and smaller than -10
rem.slopes = which(cannot.before.tolv$before > 5 | cannot.before.tolv$before < -10)
cannot.before.tolv = cannot.before.tolv[-rem.slopes,] # 218 307

# the ones that have only 2 measurement for before slopes were removed
bsused = merged.creas.final3 %>%  group_by(patient_id) %>% 
  reframe(bu = length(na.omit(which(before == slope))))
keepthem = rownames(cannot.before.tolv)[match(setdiff(cannot.before.tolv$patient_id,
                                                      bsused$patient_id[which(bsused$bu < 3)]),
                                              cannot.before.tolv$patient_id)]
cannot.before.tolv = cannot.before.tolv[keepthem,] # 214 307
data.before.tolv = data.slopesF[,rownames(cannot.before.tolv)]

#adding clinical features to annot 
cannot.before.tolv$CKD.merged.eGFR = eGFRtoStageAB(cannot.before.tolv$merged_eGFR)
cannot.before.tolv$ageGroup = ageGroupFact(cannot.before.tolv$merged_age)

# adjusting the levels of PKD1 and PKD2 genotype
cannot.before.tolv$Class_ADPKD[which(is.na(cannot.before.tolv$Class_ADPKD))] = "unknown"
cannot.before.tolv$selectedmut = as.factor(cannot.before.tolv$Class_ADPKD)
levels(cannot.before.tolv$selectedmut) = list(unknown = "unknown", 
                                              PKD1vus = "PKD1 VUS",
                                              PKD2 = "PKD2",PKD2vus = "PKD2 VUS",
                                              PKD1NT = "PKD1 non-truncating",
                                              PKD1T = "PKD1 truncating")
# adjusting th levels of mayo class
cannot.before.tolv$merged_mayo = droplevels(cannot.before.tolv$merged_mayo)
cannot.before.tolv$mayoClass = cannot.before.tolv$merged_mayo
levels(cannot.before.tolv$mayoClass) = c("1A" = "Less", "1B" = "Less", 
                                         "1C" = "Mid", "1D" = "More", 
                                         "1E" = "More")
cannot.before.tolv$mayoClass = factor(cannot.before.tolv$mayoClass,
                                      c("Less", "Mid", "More")) # 214 311

# | # Feature Selection ----
lassodat = cbind.data.frame("before" = cannot.before.tolv$before,
                            t(data.before.tolv))
lassodat.mm2 <- modmat(lassodat) # 214 283

# | # | # LIMMA  ----
slope.limma1 = limma.numeric.resultsUp(~before, cannot.before.tolv,
                                       t(na.omit(lassodat.mm2)),pval = 0.05)

#this is for heatmap
limma.sig = getNamesPro(gsub("`","",rownames(slope.limma1[[3]])))
limma.sig[which(is.na(limma.sig))] = rownames(slope.limma1[[3]])[which(is.na(limma.sig))]

#this is for intersecting it with limma
intlimma = unique(gsub("&.*", "", gsub("`","",rownames(slope.limma1[[3]])))) # 27

# # | # | # LASSO ----
library(caret)
library(glmnet)
lassodat.mm2 = as.data.frame(lassodat.mm2)
colnames(lassodat.mm2) = gsub("`","",colnames(lassodat.mm2))

# creating partitions, train and test objects
set.seed(12345)
flsNew = createDataPartition(lassodat.mm2$before, p = 2/3, times = 100)
trCont = trainControl(method = "cv", number = 10, savePredictions = "all")
testCont = trainControl(method = "none")
lambdas = 10^seq(5, -5, length = 100)

myScaledClassoPred = list()
myScaledClassoPred = lapply(flsNew, function(x){
  
  # test and training data
  train = lassodat.mm2[x,]
  test = lassodat.mm2[-x,]
  
  # training model
  modeltrain = train(before~., data = train, weights = abs(train[,"before"]),
                     method = "glmnet", na.action = na.omit, trControl = trCont,
                     tuneGrid = expand.grid(alpha = 1, lambda = lambdas))
  allcoef = coef(modeltrain$finalModel, modeltrain$bestTune$lambda)
  
  # testing model
  modeltest = train(before~., data = test,
                    method = "glmnet", na.action = na.omit, trControl = testCont,
                    tuneGrid = expand.grid(alpha = 1, 
                                           lambda = modeltrain$bestTune$lambda))
  
  # performance related measures
  predictions = predict(modeltrain, newdata = test)
  modelperf = data.frame(RMSE = RMSE(predictions, test$before), 
                         Rsq = R2(predictions, test$before))
  
  return(list("coefs" = allcoef, "modelperf" = modelperf, 
              "predictions" = predictions, "model" = modeltrain))
})

# checking the errors of the test models
hist(unlist(lapply(lapply(myScaledClassoPred, "[[", "modelperf"), function(y) y$RMSE)))
hist(unlist(lapply(lapply(myScaledClassoPred, "[[", "modelperf"), function(y) y$Rsq)))

# selected proteins in each model
selectedOnes = lapply(lapply(myScaledClassoPred, "[[", "coefs"), function(x){
  sel = gsub("`","", x@Dimnames[[1]][x@i + 1][-1])
  return(sel)
})
sort(table(unlist(selectedOnes)))
lassoresults = cbind.data.frame(sort(table(unlist(selectedOnes))), 
                                getNamesPro(names(sort(table(unlist(selectedOnes))))))

# selecting the proteins that included in 75% of the models generated
selected.lasso = gsub("&.*","",names(which(table(unlist(selectedOnes)) > 100*0.75)))

# | # Model Generation ----
# | # | # First Models ----

# extracting the proteins that selected by lasso
slope.prots.int = c(selected.lasso, "before")
beforetolvall = cbind.data.frame(t(data.before.tolv), cannot.before.tolv)
model.slopes = beforetolvall[,slope.prots.int] # 214 7

#get the gene names
colnamesnew = getNamesPro(colnames(model.slopes))
colnamesnew[which(is.na(colnamesnew))] = colnames(model.slopes)[which(is.na(colnamesnew))]
colnames(model.slopes) = colnamesnew

#to select the most stable ones, backward selection will be used
model.slopes = na.omit(model.slopes) # 214   7
model.all = lm(before~., data = model.slopes)
summary(model.all)
# Residual standard error: 2.029 on 207 degrees of freedom
# Multiple R-squared:  0.3176,	Adjusted R-squared:  0.2979 
# F-statistic: 16.06 on 6 and 207 DF,  p-value: 3.834e-15

library(MASS)
stepwise.all.rem = stepAIC(model.all, direction = "both")
summary(stepwise.all.rem) 
# Residual standard error: 2.029 on 207 degrees of freedom
# Multiple R-squared:  0.3176,	Adjusted R-squared:  0.2979 
# F-statistic: 16.06 on 6 and 207 DF,  p-value: 3.834e-15

# selected proteins --> pca
selected_fs = c(gsub(" ","",strsplit2(stepwise.all.rem$call$formula, "\\+")[3,]), "before")
selprots = na.omit(getIDspro(setdiff(selected_fs, "before"),
                             proteinIDtoGenes_v2$ProteinID,
                             proteinIDtoGenes_v2$Genes))

# | # | # Cross Validating the Models ----
# setting CVs
library(caret)
set.seed(95200295)
model.slopes.new = na.omit(model.slopes[,selected_fs]) #214   7
summary(lm(before~. ,model.slopes.new))
# Residual standard error: 2.029 on 207 degrees of freedom
# Multiple R-squared:  0.3176,	Adjusted R-squared:  0.2979 
# F-statistic: 16.06 on 6 and 207 DF,  p-value: 3.834e-15

# creating folds
flsH = createDataPartition(model.slopes.new$before, times = 100, p = 2/3)
flsH = lapply(flsH, function(x) rownames(model.slopes.new)[x])

# | # | # | # Linear Regression ----
LR.proteins = lapply(flsH, function(x){
  modelgen.proteome(dat = model.slopes.new[,selected_fs], 
                    dependent = "before", stepwise = F,
                    trainmi = T, preview = F, folds = x, model = "lr")
})
LR.proteins.rmse = mean(unlist(lapply(LR.proteins, function(x) x$rmse)))
hist(sapply(LR.proteins,"[[", "Accuracy.tests"))
plot(density(sapply(LR.proteins,"[[", "Accuracy.tests")))

# | #  Prediction ----
# predict with Proteome Model
selected.model = modelgen.proteome(dat = model.slopes.new[,selected_fs],
                                   dependent = "before", stepwise = F,
                                   trainmi = F, preview = F, model = "lr")
file = "~/Data/slope-models-adpkd/data/results/"
save(selected.model, file = paste0(file,"mod_pm.RData"))
selected.model.pred = predict(selected.model$`LR model`, model.slopes.new)

model.slopes.new$age = cannot.before.tolv$merged_age[match(rownames(model.slopes.new),rownames(cannot.before.tolv))]
model.slopes.new$gender = cannot.before.tolv$merged_gender[match(rownames(model.slopes.new),rownames(cannot.before.tolv))]

model.slopes.new$mayo = cannot.before.tolv$mayoClass[match(rownames(model.slopes.new),rownames(cannot.before.tolv))]
model.slopes.new$eGFR = cannot.before.tolv$merged_eGFR[match(rownames(model.slopes.new),rownames(cannot.before.tolv))]

model.slopes.new$selectedmut = cannot.before.tolv$selectedmut[match(rownames(model.slopes.new),rownames(cannot.before.tolv))]
levels(model.slopes.new$selectedmut)

# predict with Combined Model
selected.model.ext = lm(before~SERPINF1+GPX3+AFM+FERMT3+CFHR1+RARRES2+eGFR+gender+age+mayo, 
                        data = model.slopes.new)
summary(selected.model.ext)
save(selected.model.ext, file = paste0(file,"mod_cm.RData"))
selected.model.ext.pred = predict(selected.model.ext, model.slopes.new)

# predict with Combined Genotype Model
selected.model.extGen = lm(before~SERPINF1+GPX3+AFM+FERMT3+CFHR1+RARRES2+eGFR+gender+age+mayo+selectedmut,
                           data = subset(model.slopes.new, !grepl("vus|unknown", model.slopes.new$selectedmut, ignore.case = T)))
summary(selected.model.extGen)
save(selected.model.extGen, file = paste0(file,"mod_cgm.RData"))
selected.model.extGEN.pred = predict(selected.model.extGen, subset(model.slopes.new, !grepl("vus|unknown", model.slopes.new$selectedmut, ignore.case = T)))

# predict with Proteome4 Model
selected.model.min = lm(before~SERPINF1+GPX3+AFM+CFHR1, data = model.slopes.new)
summary(selected.model.min) #adjR2 0.2506 
save(selected.model.min, file = paste0(file,"mod_p4m.RData"))
selected.model.min.pred = predict(selected.model.min, model.slopes.new)
save(selected.model.min.pred, file = paste0(file,"pred_p4m.RData"))

# | # Accuracy ----
lr.fold.acc = cbind.data.frame("LRaccuracy" = sapply(LR.proteins,"[[", "Accuracy.tests"))

# Plot Accuracies of Test Folds
xy.limits <- range(c(lr.fold.acc$LRaccuracy))
xy.limits[1] = xy.limits[1] - 0.05
xy.limits[2] = xy.limits[2] + 0.05

ggplot(lr.fold.acc) +
  geom_density(aes(x = LRaccuracy, colour = "#8c50a0"), fill = "#8c50a0", 
               alpha = 0.6, linewidth = 1) +
  scale_x_continuous(limits = xy.limits) + theme_bw(base_size = 20) +
  scale_colour_manual('Models', limits = c("Linear Regression"), values = c('#8c50a0')) +
  scale_colour_manual("Models", values = c('#8c50a0'), labels = c("Linear Regression")) +
  theme_bw() + xlab(expression(""~R^2)) + ylab("Density") +
  geom_vline(xintercept = density(lr.fold.acc$LRaccuracy)$x[which.max(density(lr.fold.acc$LRaccuracy)$y)],
             color = "#8c50a0", size = 1, linetype = 2) +
  theme(legend.position = 'none')
ggsave(paste0(file,"FigureS5.png"),width = 10, height = 7 , units = "in")

# | # Future eGFRs ----
# preparinf the data
model.slopes.new$predSLOPEs = selected.model.pred
model.slopes.new$predSLOPEsexted = selected.model.ext.pred
model.slopes.new[names(selected.model.extGEN.pred),"predSLOPEsextedGEN"] = selected.model.extGEN.pred
model.slopes.new$IDD = rownames(model.slopes.new)
model.slopes.new$diff = model.slopes.new$before - model.slopes.new$predSLOPEs

model.slopes.new$merged_CKD = cannot.before.tolv$CKD.merged.eGFR[match(rownames(model.slopes.new),rownames(cannot.before.tolv))]
model.slopes.new$merged_CKD2 = model.slopes.new$merged_CKD
model.slopes.new$merged_CKD2[grep("2|3|4", model.slopes.new$merged_CKD2)] = "2-4"
levels(model.slopes.new$merged_CKD2) = c("1" = "1", "2-4" = "2", "2-4" = "3a",
                                         "2-4" = "3b", "2-4" = "4")

model.slopes.new$identifier = cannot.before.tolv$identifier[match(model.slopes.new$IDD, 
                                                  cannot.before.tolv$ID.sc)]
model.slopes.new$mayo5 = cannot.before.tolv[match(model.slopes.new$IDD, 
                                                  cannot.before.tolv$ID.sc), "merged_mayo"]

# for calculating future eGFR --> you need to get positive time difference
calctimedif = sapply(model.slopes.new$identifier, function(x){
  calcTimeDiff(x, merged.creas.final3, cannot.slopesF, 
               "patient_id", "identifier", "exam_date", 
               "protSamp",eGFR = "patient__eGFR")
},simplify = FALSE,USE.NAMES = TRUE)

# to select indexes after proteome sampling
indxx = sapply(sapply(calctimedif, "[[", "timedif"), function(x) x > 0) 

# to select time differences after proteome sampling
timediffall = sapply(sapply(calctimedif, "[[", "timedif"), function(x) x[which(x > 0)]) 

# to select eGFRs after proteome sampling
eGFRsclc = mapply(function(x,y){ 
  x[y]
},x = sapply(calctimedif, "[[", "eGFRs"),y = indxx)

# Future eGFR calculation with Proteome Model
onlyProtFuture = mapply(function(x,y){
  calcFuteGFR(timedif = x, id = y, idcol = "identifier" , data = model.slopes.new, 
              which = "prot slopes", predtype = "predSLOPEs")
},x = timediffall, y = names(timediffall))

# Future eGFR calculation with Combined Model
extendedProtFuture = mapply(function(x,y){
  calcFuteGFR(timedif = x, id = y, idcol = "identifier" , data = model.slopes.new, 
              which = "prot slopes", predtype = "predSLOPEsexted")
},x = timediffall, y = names(timediffall))

# Future eGFR extraction from the annotation file
actualFuture = mapply(function(x,y){
  merged.creas.final3[y[x],"patient__eGFR"]
},x = indxx, y = sapply(calctimedif, "[[", "index"))

# Future eGFR calculation with MIC Model
mayoFuture = mapply(function(x,y){
  calcFuteGFR(timedif = x, id = y, idcol = "identifier" ,
              data = model.slopes.new, which = "mayo slopes")
},x = timediffall, y = names(timediffall))

remthese = names(which(sapply(mayoFuture, is.null) | unlist(lapply(sapply(mayoFuture, is.na), any))))
mayoFuture[remthese] <- NA

scatterFuture = do.call(rbind.data.frame, lapply(names(timediffall), function(x){
  # if a patient & their time difference data & predicted future eGFR from
  # predicted slopes (Proteome Model and Combined Model) & actual eGFR 
  # measurement exist 
  if (!is.null(x) & !is.null(as.numeric(timediffall[[x]])) & !is.null(onlyProtFuture[[x]]) & !is.null(extendedProtFuture[[x]]) & !is.null(actualFuture[[x]])) {
    # generate a vector for that patient containing all these variables
    # and future eGFR from mayo slope
    return(cbind.data.frame(id = x, timeidffyear = as.numeric(timediffall[[x]]), 
                            onlyProtFuteGFR = onlyProtFuture[[x]], 
                            extFuteGFR = extendedProtFuture[[x]], 
                            realFuteGFR = actualFuture[[x]], 
                            mayoFut = mayoFuture[[x]]))
  }
})) # and then combine them row-wise

scatterFutureMaxd = scatterFuture %>% group_by(id) %>% 
  filter(timeidffyear == max(timeidffyear)) #selecting the maximum time difference
scatterFutureMaxd2 = scatterFutureMaxd
scatterFutureMaxd2$deltaSonlyP = scatterFutureMaxd2$realFuteGFR - scatterFutureMaxd2$onlyProtFuteGFR #delta eGFR Proteome Model
scatterFutureMaxd2$deltaSPext = scatterFutureMaxd2$realFuteGFR - scatterFutureMaxd2$extFuteGFR #delta eGFR Combined Model
scatterFutureMaxd2$deltaSPmayo = scatterFutureMaxd2$realFuteGFR - scatterFutureMaxd2$mayoFut #delta eGFR MIC Model

scatterFutureMaxd2[,c("extFuteGFR","onlyProtFuteGFR","realFuteGFR","mayoFut")] = NULL

# melting the data for ggplot
scatterFutureMaxd2 = melt(scatterFutureMaxd2, id.vars = c("timeidffyear","id"))

# assigning limits
ylimits.ed = range(na.omit(c(scatterFutureMaxd$realFuteGFR - scatterFutureMaxd$onlyProtFuteGFR,
                             scatterFutureMaxd$realFuteGFR - scatterFutureMaxd$extFuteGFR,
                             scatterFutureMaxd$realFuteGFR - scatterFutureMaxd$mayoFut))) + c(-2,2)

library(rcartocolor)
# 4 samples from 2 patients were removed --> they don't have mayo class 
scatterFOnlyProtFig = ggplot(scatterFutureMaxd2, aes(x = timeidffyear, y = value)) + 
  geom_point(aes(color = variable), size = 3, alpha = 0.7) + 
  scale_y_continuous(limits = ylimits.ed, breaks = seq(-50, 50, by = 10)) + 
  theme_bw(base_size = 17) + scale_x_continuous(breaks = seq(0,8,1)) + 
  geom_hline(yintercept = 0, linewidth = 2, colour = "darkred", linetype = 3) +
  scale_colour_manual(name = "Models", 
                      values = c("#3C97DA", "#165079", "#FE9D52"), 
                      labels = c("Proteome Model", "Combined Model", "MIC Model")) +
  scale_fill_manual(name = "Models", 
                    values = c("#3C97DA", "#165079", "#FE9D52"), 
                    labels = c("Proteome Model", "Combined Model", "MIC Model")) + 
  xlab(expression("Latest eGFR for each patient (in years)")) + 
  ylab(expression(~Delta~"eGFR"["Observed-Predicted"])) +
  stat_smooth(method = "lm", se = T, aes(fill = variable, colour = variable), 
              alpha = 0.3, level = 0.95) +
  stat_smooth(method = "lm", se = F, aes(fill = variable, colour = variable)) +
  theme(legend.position = 'bottom', 
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(size = 16))
scatterFOnlyProtFig

# MIC Model - Slope calculation 
mayoIds = model.slopes.new[which(!is.na(model.slopes.new$mayo)), "identifier"]
time1to7 = rep_len(list(c(1:7)), length(mayoIds))
names(time1to7) = mayoIds
mayoFtime1to7 = mapply(function(x,y){
  calcFuteGFR(timedif = x, id = y, idcol = "identifier" , 
              data = model.slopes.new, which = "mayo slopes")
},x = time1to7, y = names(time1to7), SIMPLIFY = F)

mayoSlopes = mapply(function(eGFR, year){
    lm(eGFR~year)$coefficients[2] #this is already years 
}, eGFR = mayoFtime1to7, year = time1to7[names(mayoFtime1to7)])

mayoSlopes.df = cbind.data.frame("mayoSlope" = unlist(mayoSlopes), 
                                 "ID" = names(mayoSlopes),
                                 model.slopes.new[match(gsub("\\..*", "",names(mayoSlopes)), 
                                                        model.slopes.new$identifier),
                                                c("before","predSLOPEs","predSLOPEsexted",
                                                  "identifier", "age", "merged_CKD2", "gender", "mayo")])
library(caret)
mayoSlopes.df$diff = mayoSlopes.df$before - mayoSlopes.df$mayoSlope
countmayo = mayoSlopes.df %>% group_by(merged_CKD2) %>% #generating this for boxplot
  reframe(rmse = paste0("RMSE=", round(RMSE(pred = mayoSlope, obs = before),2)), 
          n = paste0("n=",n()))

#adding octreotid usage as a column for putting stars
ggplotSlopes = cbind.data.frame(model.slopes.new[,c("predSLOPEs","predSLOPEsexted",
                                                    "predSLOPEsextedGEN","before",
                                                    "merged_CKD2","identifier")], "octreotid" = "no")
ggplotSlopes$octreotid[match(intersect(octreotid, ggplotSlopes$identifier),ggplotSlopes$identifier)] = "yes"

# | # Proteome Model ----
# | # | # Plot Observed-Predicted Slopes ----
xy.limits <- range(na.omit(c(ggplotSlopes$before, ggplotSlopes$predSLOPEs, 
                             ggplotSlopes$predSLOPEsexted)))
library(grid)
scatterPredvsObs = ggplot(data = ggplotSlopes, aes(x = before, y = predSLOPEs, 
                                                 colour = merged_CKD2)) +
  geom_point(size = 3) + #, aes(shape=octreotid))+
  theme_bw(base_size = 17) +
  #geom_point(size=2)+
  scale_colour_manual(name = "CKD Stages", values = c("#e5c779", "#d43f96")) +
  #scale_shape_manual(name="Octreotid Usage", values=c(16,8)) +
  scale_x_continuous(limits = xy.limits) + scale_y_continuous(limits = xy.limits) +
  xlab("Observed Slope") + ylab("Predicted Slope") +
  geom_abline(intercept = 0, slope = 1, linewidth = 1, 
              linetype = "dashed", colour = "gray50") +
  stat_smooth(method = "rlm", se = T, fill = 4, colour = 4) +
  theme(legend.position = 'bottom', 
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(size = 16))

# | # | # CKD1 va CKD2-4 Plot Observed-Predicted Slopes ----
library(caret)
distance.list = lapply(c("1", "2-4"), function(x){
  preds = predict(selected.model$`LR model`, model.slopes.new[model.slopes.new$merged_CKD2 == x,])
  rmse = RMSE(preds, model.slopes.new[model.slopes.new$merged_CKD2 == x,"before"])
  distance = model.slopes.new[model.slopes.new$merged_CKD2 == as.character(x),"before"] - preds
  quant = quantile(distance, probs = c(0, 0.33, 0.66, 1))
  data.frame("distance" = distance, "p0" = rep(quant[1], length(distance)), 
             "p33" = rep(quant[2], length(distance)), 
             "p66" = rep(quant[3], length(distance)),
             "p100" = rep(quant[4], length(distance)),
             "RMSE" = rep(round(rmse, digits = 2), 
                          length(which(model.slopes.new$merged_CKD2 == x))), 
             row.names = rownames(model.slopes.new[model.slopes.new$merged_CKD2 == as.character(x),]),
             "observed" = model.slopes.new[model.slopes.new$merged_CKD2 == x,"before"], 
             "predicted" = preds)
})

names(distance.list) = c("1", "2-4")
distance.df = do.call(rbind.data.frame, distance.list)
distance.df = cbind.data.frame(distance.df, ckd = gsub("\\..*","",rownames(distance.df)))
rownames(distance.df) = gsub(".*\\.","",rownames(distance.df))

library(dplyr)
countss = dplyr::count(distance.df, ckd, RMSE, p0, p33, p66, p100)
countss$n = paste0("n=", countss$n)
countss$RMSE = paste0("RMSE=", countss$RMSE)

# | # Combined Model ----
# | # | # Plot Observed-Predicted Slopes ----
scatterPredvsObsExt = ggplot(data = ggplotSlopes, aes(x = before,
                                                      y = predSLOPEsexted,
                                                      colour = merged_CKD2)) +
  geom_point(size = 3) +#, aes(shape=octreotid))+
  theme_bw(base_size = 17) +
  #geom_point(size=2) +
  scale_colour_manual(name = "CKD Stages", values = c("#e5c779", "#d43f96")) +
  #scale_shape_manual(name="Octreotid Usage", values=c(16,8)) +
  scale_x_continuous(limits = xy.limits) + scale_y_continuous(limits = xy.limits) +
  xlab("Observed Slope") + ylab("Predicted Slope") +
  geom_abline(intercept = 0, slope = 1, linewidth = 1, 
              linetype = "dashed", colour = "gray50") +
  stat_smooth(method = "rlm", se = T, fill = 4, colour = 4) +
  theme(legend.position = 'bottom', legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(size = 16))

# | # | # CKD1 va CKD2-4 Plot Observed-Predicted Slopes ----
model.slopes.new.ext = model.slopes.new
model.slopes.new.ext = model.slopes.new.ext[-which(is.na(model.slopes.new.ext$mayo)),]

distance.list.ext = lapply(c("1", "2-4"), function(x){
  preds = predict(selected.model.ext, model.slopes.new.ext[model.slopes.new.ext$merged_CKD2 == x,])
  rmse = RMSE(preds, model.slopes.new.ext[model.slopes.new.ext$merged_CKD2 == x,"before"])
  distance = model.slopes.new.ext[model.slopes.new.ext$merged_CKD2 == as.character(x),"before"] - preds
  quant = quantile(distance, probs = c(0, 0.33, 0.66, 1))
  data.frame("distance" = distance, "p0" = rep(quant[1], length(distance)),
             "p33" = rep(quant[2], length(distance)), 
             "p66" = rep(quant[3], length(distance)),
             "p100" = rep(quant[4], length(distance)),
             "RMSE" = rep(round(rmse, digits = 2), 
                          length(which(model.slopes.new.ext$merged_CKD2 == x))),
             row.names = rownames(model.slopes.new.ext[model.slopes.new.ext$merged_CKD2 == as.character(x),]),
             "observed" = model.slopes.new.ext[model.slopes.new.ext$merged_CKD2 == x,"before"], 
             "predicted" = preds)
})

names(distance.list.ext) = c("1", "2-4")
distance.df.ext = do.call(rbind.data.frame, distance.list.ext)
distance.df.ext = cbind.data.frame(distance.df.ext, ckd = gsub("\\..*","",rownames(distance.df.ext)))
rownames(distance.df.ext) = gsub(".*\\.","",rownames(distance.df.ext))

countss.ext = dplyr::count(distance.df.ext, ckd, RMSE, p0, p33, p66, p100)
countss.ext$n = paste0("n=", countss.ext$n)
countss.ext$RMSE = paste0("RMSE=", countss.ext$RMSE)

####

aonlyprot.bp = left_join(distance.df[,c("ckd","distance")],
                         countss[,c("ckd","RMSE","n")], by = "ckd")
bextenprot.bp = left_join(distance.df.ext[,c("ckd","distance")],
                          countss.ext[,c("ckd","RMSE","n")], by = "ckd")
cmayoslopes.bp = left_join(mayoSlopes.df[,c("merged_CKD2","diff")],
                           countmayo[,c("merged_CKD2","rmse","n")], 
                           by = "merged_CKD2")
colnames(cmayoslopes.bp) = colnames(aonlyprot.bp)

df_names = c("aonlyprot.bp", "bextenprot.bp","cmayoslopes.bp")
combinedBoxPlots = do.call(rbind.data.frame, lapply(df_names, function(x){
  data.frame(id = x, eval(parse(text = x)))}))

ylimitsbp = range(na.omit(combinedBoxPlots$distance))
ylimitsbp = ylimitsbp + c(-2,0.5) #changed to have more space for text 

combinedBoxPlotsFigC = ggplot(data = combinedBoxPlots, aes(y = distance)) +
  geom_boxplot(aes(x = id, colour = factor(ckd)), lwd = 1) + 
  theme_bw(base_size = 17) +
  scale_x_discrete(labels = c("Proteome Model", "Combined Model", "MIC Model")) +
  geom_text(aes(y = min(ylimitsbp) + 0.8, label = n, x = id, 
                colour = factor(ckd)), fontface = "bold", 
            position = position_dodge(0.75)) + #changed to have more space for text 
  geom_text(aes(y = min(ylimitsbp) + 0.2, label = RMSE, x = id, 
                colour = factor(ckd)), fontface = "bold", 
            position = position_dodge(0.75)) +
  scale_y_continuous(limits = ylimitsbp) +
  theme(strip.background = element_rect(colour = "black",fill = "white", linewidth = 0)) +
  scale_colour_manual(name = "CKD Stages", values = c("#e5c779", "#d43f96")) +
  xlab("Models") + ylab(expression(~Delta~"Slope"["Observed-Predicted"])) +
  theme(legend.position = 'bottom', 
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(size = 16))

library(MASS)
library(ggpubr)
ggarrange(ggarrange(scatterPredvsObs, scatterPredvsObsExt,
                    ncol = 2,  nrow = 1, labels = c("A", "B"), legend = "none"), 
          ggarrange(combinedBoxPlotsFigC, scatterFOnlyProtFig,
                    ncol = 2,  nrow = 1, labels = c("C", "D")),
          common.legend = T, nrow = 2, legend = "bottom")
ggsave(paste0(file,"Figure2.png"), width = 17, height = 17 , units = "in", dpi = "retina",)

# | # Model Tables ----
library(sjPlot)
proteome.model = lm(before~SERPINF1+GPX3+AFM+FERMT3+CFHR1+RARRES2, data = model.slopes.new)

# TableS6
table3 = tab_model(proteome.model, dv.labels = c("Proteome Model"))

proteome.ext.model = lm(before~SERPINF1+GPX3+AFM+FERMT3+CFHR1+RARRES2+age+gender+eGFR+mayo, 
                        data = model.slopes.new)
tab_model(proteome.ext.model, dv.labels = c("Combined Model"))

clinmodel = lm(before~age+gender+eGFR+mayo, data = model.slopes.new)
tab_model(clinmodel, dv.labels = c("Clinical Model"))

# to make the Proteome and Combined Models comparable, 
# 4 samples missing mayo class information are removed 
obs212 = intersect(names(proteome.ext.model$fitted.values), names(clinmodel$fitted.values))

proteome.model212 = lm(before~SERPINF1+GPX3+AFM+FERMT3+CFHR1+RARRES2, 
                       data = model.slopes.new[obs212,])
clinmodel212 = lm(before~age+gender+eGFR+mayo, data = model.slopes.new[obs212,])
proteome.ext.model212 = lm(before~SERPINF1+GPX3+AFM+FERMT3+CFHR1+RARRES2+age+gender+eGFR+mayo, 
                           data = model.slopes.new[obs212,])

# Table2
table4 = tab_model(proteome.model212,clinmodel212, proteome.ext.model212,
                   dv.labels = c("Proteome Model", "Clinical Model", 
                               "Combined Model"),
                   pred.labels = c("(Intercept)","SERPINF1", "GPX3", "AFM", 
                                   "FERMT3", "CFHR1", "RARRES2", "Age", 
                                   "Sex [Male]", "eGFR", "MAYO [Mid]", 
                                   "MAYO [More]"))

combGenmodel = lm(before~SERPINF1+GPX3+AFM+FERMT3+CFHR1+RARRES2+age+gender+eGFR+mayo+selectedmut,
                  data = subset(model.slopes.new, !grepl("vus|unknown", model.slopes.new$selectedmut,
                                                         ignore.case = T)))

# to make the Proteome, Combined and Combined Genotype Models comparable, 
# 100 samples missing either mayo class and genotype information are removed 
obs114 = intersect(names(combGenmodel$fitted.values), names(clinmodel$fitted.values))

# subsetting the models according to samples which have 
# mayo class and genotype information
proteome.model114 = lm(before~SERPINF1+GPX3+AFM+FERMT3+CFHR1+RARRES2, 
                       data = model.slopes.new[obs114,])
clinmodel114 = lm(before~age+gender+eGFR+mayo, data = model.slopes.new[obs114,])
clinmodel1143 = lm(before~age+gender+eGFR+mayo+selectedmut, 
                   data = model.slopes.new[obs114,])
proteome.ext.model114 = lm(before~SERPINF1+GPX3+AFM+FERMT3+CFHR1+RARRES2+age+gender+eGFR+mayo, 
                           data = model.slopes.new[obs114,])
combGenmodel114 = lm(before~SERPINF1+GPX3+AFM+FERMT3+CFHR1+RARRES2+age+gender+eGFR+mayo+selectedmut, 
                   data = model.slopes.new[obs114,])

# TableS7
tab_model(proteome.model114,clinmodel114,clinmodel1143, 
          proteome.ext.model114, combGenmodel114,
          dv.labels = c("Proteome Model", "Clinical Model", 
                        "Clinical Genotype Model", "Combined Model", 
                        "Combined Genotype Model"),
          pred.labels = c("(Intercept)","SERPINF1", "GPX3", "AFM", "FERMT3", 
                          "CFHR1", "RARRES2", "Age", "Sex [Male]", "eGFR", 
                          "MAYO [Mid]", "MAYO [More]", "PKD1 [NT]", "PKD1 [T]"))

# TableS9
tab_model(selected.model.min, dv.labels = c("Proteome4 Model"))

# | # Heatmaps ----
library(ComplexHeatmap)
library(circlize)
library(makeunique)

# preparing the sata for heatmap
hm.slopes = t(data.before.tolv[unique(c(selprots,intlimma)),])
colnames(hm.slopes) = make_unique(getNamesPro(colnames(hm.slopes)), sep = "#",
                                  wrap_in_brackets = F)
# generating colors to be used for column names in the heatmap
color.cofs = rep("black", length(colnames(hm.slopes)))
color.cofs[match(c("SERPINF1", "GPX3", "AFM", "FERMT3", "CFHR1", "RARRES2"), 
                 colnames(hm.slopes))] = "red"
names(color.cofs) <- paste(colnames(hm.slopes))

# calculating the mean slope for the heatmap legend
midSlope = mean(cannot.before.tolv[rownames(hm.slopes),"before"], na.rm = T)

# preparing the row annotation for heatmap
ha = rowAnnotation(
  Slopes = cannot.before.tolv[rownames(hm.slopes),"before"],
  eGFR = cannot.before.tolv[rownames(hm.slopes), "merged_eGFR"],
  MAYO = cannot.before.tolv[rownames(hm.slopes), "mayoClass"],
  Age = cannot.before.tolv[rownames(hm.slopes), "merged_age"],
  Sex = cannot.before.tolv[rownames(hm.slopes),"merged_gender"],
  
  col = list(Slopes = colorRamp2(c(-8.5,midSlope, 3.7), c("#004D66", "#0AC2FF", "#C2F0FF")),
             eGFR = colorRamp2(c(12, 72.5, 133), c("#4b2991", "#f7667c", "#edd9a3")),
             MAYO = c("Less" = "#d9ed92", "Mid" = "#52b69a", "More" = "#1e6091"),
             Age = colorRamp2(c(18.7, 40, 79.1), c("#e4c7f1", "#9999EA", "#4545D9")),
             Sex = c("female" = "pink", "male" = "cornflowerblue")
  ))
# | # | # Original ----
library(ComplexHeatmap)
library(circlize)
set.seed(1234567)
hmclustOrig = draw(Heatmap(hm.slopes,
                           column_names_gp = gpar(cex = 1, col = color.cofs, 
                                                  fontface = "bold"),
                           clustering_distance_rows = function(x) as.dist(1 - cor(t(x))),
                           clustering_method_rows = "average",
                           clustering_distance_columns = function(x) as.dist(1 - cor(t(x))),
                           clustering_method_columns = "average",
                           row_km = 3, column_km = 3,
                           row_dend_reorder = T, column_dend_reorder = T,
                           show_row_dend = F,
                           #col = colorRamp2(c(-3,0,3), 
                                            #c("green", "black", "red")),
                           heatmap_legend_param = list( title = "Intensity"),
                           show_row_names = F, right_annotation = ha))

# | # | # Sorted according to eGFR ----
hmclustOrigSorted = hmclustOrig

#obtaining row and column orders
clustrows = row_order(hmclustOrig)
clustcol = column_order(hmclustOrig)

rowmylcustdata = data.frame(cluster = rep(names(clustrows),
                                          lengths(clustrows)),
                            id = unlist(clustrows))
rowmylcustdata$batchID = rownames(hm.slopes)[rowmylcustdata$id]
rowmylcustdata$eGFR = cannot.before.tolv[rowmylcustdata$batchID,"merged_eGFR"]
rowmylcustdata$m = cannot.before.tolv[rowmylcustdata$batchID,"mayoClass"]
rownames(rowmylcustdata) = rowmylcustdata$batchID

ord_rmcd = rowmylcustdata[order(rowmylcustdata$cluster, rowmylcustdata$eGFR),]
ord_rmcdspt = split(ord_rmcd, ord_rmcd$cluster)
hmclustOrigSorted@ht_list[[names(hmclustOrigSorted@ht_list)]]@row_order_list = lapply(ord_rmcdspt, function(x) x$id)[c(2,1,3)]

pdf(paste0(file, "Figure1.pdf"), width = 10, height = 8,
    bg = "white", paper = "A4")
hmclustOrigSorted
dev.off()

# check correlation between CST3 and proteins
cst3cordn = model.slopes.new
cst3cordn2 = cbind.data.frame(cst3cordn, "CST3" = data.before.tolv["P01034", rownames(cst3cordn)])
corCST3 = cor(cst3cordn2[c(selected_fs[-7], "CST3", "eGFR", "before")], 
              use = "complete.obs") # max 0.39568677 - SERPINF1
write.csv(corCST3, file = paste0(file,"TableS10.csv"))

#all 29 markers - cor with before and eGFR
allcompcor = cbind.data.frame(cst3cordn, hm.slopes[rownames(cst3cordn),])

library(Hmisc)
allcompcorSC <- rcorr(as.matrix(allcompcor[c(colnames(hm.slopes), "eGFR")]))
cor_all <- allcompcorSC$r
p_all <- allcompcorSC$P

cortabsSC = data.frame("corGFR" = cor_all[,"eGFR"], "pGFR" = p_all[,"eGFR"],
                       "sigstarGFR" = sigstar(p_all[,"eGFR"]))
cortabsSC
write.csv(cortabsSC, file = paste0(file,"TableS4.csv"))

# | # | # Enrichment of Proteins ----
# | # | # | # 29 Proteins against all detected proteins ----
library(gprofiler2)
# setting archived url for gprofiler
set_base_url("https://biit.cs.ut.ee/gprofiler_archive3/e111_eg58_p18")

# creating a key for protein IDs and names in heatmap
hm.slopesKEY = data.frame(protnam = unique(c(selprots,intlimma)),
                          genenam = getNamesPro(unique(c(selprots,intlimma)))) 
# 29  2

# splitting the protein IDs containing ;
split_pro_hm <- strsplit(hm.slopesKEY$protnam, ";")
split_pro_hm_unls <- unlist(split_pro_hm)
hm.slopesKEY <- data.frame(proteinID = split_pro_hm_unls, 
                           geneName = rep(hm.slopesKEY$genenam,  
                                          lengths(split_pro_hm))) # 35  2

#converting protein IDs to gene ENSG IDs 
hm_ensgID = gconvert(hm.slopesKEY$proteinID, 
                     organism = "hsapiens", 
                     target = "ENSG") # 38 7

#some protein IDs match with more than one protein (due to ambiguity?) 
hm_ensgID$matched = hm.slopesKEY[match(hm_ensgID$input,hm.slopesKEY$proteinID),
                                 "geneName"] 
hm_ensgID = hm_ensgID[-which(hm_ensgID$matched != hm_ensgID$name),] # 34 8
hm_ensgID_unique = unique(hm_ensgID$target) # 28

# creating a key for all proteins that are detected, 
# even though they are detected in less than 20% of the samples
detectedprots = data.frame("ProteinID" = rownames(data.slopesFold),
                           "Genes" = getNamesPro(rownames(data.slopesFold)))
# 398   2

# splitting the protein IDs containing ;
split_pro_all <- strsplit(detectedprots$ProteinID, ";") 
split_pro_all_unls <- unlist(split_pro_all)
detectedprots <- data.frame(proteinID = split_pro_all_unls, geneName = rep(detectedprots$Genes, lengths(split_pro_all)))
#566   2

#converting protein IDs to gene ENSG IDs 
all_ensgID = gconvert(detectedprots$proteinID, 
                      organism = "hsapiens", 
                      target = "ENSG") # 603   7
all_ensgID$matched = sapply(all_ensgID$input, function(x) detectedprots$geneName[grep(paste0("\\b",x,"\\b"), detectedprots$proteinID)])
all_ensgID$tf = NA
all_ensgID$tf = mapply(function(name, matched){# creating a true false variable - tf
  # find a match between the name came up from converter and 
  # the name from the proteome data containing ";"
  grpfnc = grep(paste0("\\b",name,"\\b"), matched) # this check for matching whole words
  ifelse(length(grpfnc) > 0, T, F) # if there is a match, assign tf to true or false
},name = all_ensgID$name, all_ensgID$matched) # 603   9

# only novels contain first occurrence and they are only occurrence
none_novels = which(grepl("^None$", all_ensgID$name, ignore.case = T) & grepl("novel transcript|novel gene|novel protein",all_ensgID$description))
linc_orf_as = which(grepl("^LINC", all_ensgID$name) | grepl("open reading frame|antisense",all_ensgID$description))

all_ensgIDrefined = all_ensgID[-c(unique(none_novels, linc_orf_as)),] # 575   9

# the ones that match with our gene names
changereq_df = all_ensgIDrefined[all_ensgIDrefined$tf == T,] # 496   9
changereq = changereq_df[which(gsub(".*\\.", "", changereq_df$target_number) != "1"),"input"] # 3

wbchange_enr = grep(paste(changereq, collapse = "|"), all_ensgIDrefined$input) # 6
wbmerged_enr = all_ensgIDrefined[wbchange_enr,] # getting those 6
all_ensgIDrefined = all_ensgIDrefined[-wbchange_enr,]  # 569   9 
# and removing them from data

wbmerged_enr$tf = c(F,T,T,T,F,T) # changing those 3 proteins to get additional matches 
wbmerged_enr = wbmerged_enr[wbmerged_enr$tf == T,] # and remove the false matches 
# 4 matches stayed

#523   9 # getting the first hits
all_ensgIDrefined = all_ensgIDrefined[grep("\\.1", all_ensgIDrefined$target_number),] 

#527   9 # merging it with the 3 protein hits
all_ensgIDrefined = rbind.data.frame(all_ensgIDrefined, wbmerged_enr) 

# performing enrichment
set.seed(12345)
hm_gonew = gost(query = hm_ensgID_unique, organism = "hsapiens", 
                ordered_query = F,  custom_bg = unique(all_ensgIDrefined$target), 
                user_threshold = 0.1, correction_method = "fdr", evcodes = T)
hm_gonew.res = as.data.frame(hm_gonew$result)
hm_gonew.res = apply(hm_gonew.res,2,as.character) # 15 16
write.csv(x = hm_gonew.res, file = paste0(file, "DataS3.csv"))

hm_gonew.res.sim = getrrvgo(hm_gonew)
sort(table(hm_gonew.res.sim$reducedTerms$parentTerm),decreasing = T)  # 3
hm_gtp = unique(cbind.data.frame(hm_gonew.res.sim$reducedTerms$parent, 
                                 hm_gonew.res.sim$reducedTerms$parentTerm))
# write.csv(x = hm_gtp, file = paste0(file, "parent_enriched_GOBP_DataS3.csv"))

# | # | # | # 3 cluster of 29 Proteins against all detected proteins ----
# clustered proteins - enrichment
clusteredprots = lapply(clustcol, function(x) colnames(hm.slopes)[x])
finalclspr = lapply(clusteredprots, function(x){
  x = gsub("#.*", "", x)
  hm.slopesKEY$proteinID[match(x,hm.slopesKEY$geneName)]
}) #getting protein IDs

#getting the ENSG IDs
finalclspr2 = lapply(finalclspr, function(x) gconvert(x, organism = "hsapiens", target = "ENSG"))
finalclspr2 = lapply(finalclspr2, function(x) x[grep("\\.1", x$target_number),])
finalclspr2 = lapply(finalclspr2, function(x) x$target)

# performing enrichment
set.seed(123)
finalclsprGO = lapply(finalclspr2,function(x){
  gost(query = x, organism = "hsapiens", 
       custom_bg = unique(all_ensgIDrefined$target),
       ordered_query = F, user_threshold = 0.1, 
       correction_method = "fdr", evcodes = T)
})

# extracting GO:BP and KEGG terms
enrichmentClusters = do.call(rbind.data.frame, lapply(names(finalclsprGO), function(x){
  y = finalclsprGO[[x]]$result
  y$cluster = as.character(x)
  y = y[grep("GO:BP|KEGG", y$source),]
  return(y)}))
enrichmentClusters = as.data.frame(enrichmentClusters)
enrichmentClusters = apply(enrichmentClusters,2,as.character)
write.csv(x = enrichmentClusters, file = paste0(file, "DataS4.csv"))
# 217  17

# getting the parent enriched GO:BP terms
enrichmentClusters.sim = lapply(finalclsprGO, getrrvgo)
lapply(enrichmentClusters.sim, function(x) sort(table(x$reducedTerms$parentTerm),
                                                decreasing = T))

hmclust_gtp = lapply(enrichmentClusters.sim, function(x){
  unique(cbind.data.frame(x$reducedTerms$parent,x$reducedTerms$parentTerm))})
hmclust_gtp = do.call(rbind.data.frame, hmclust_gtp)
write.csv(x = hmclust_gtp, file = paste0(file, "DataS5.csv"))

# | # | # Enrichment of Samples ----
clusteredRowsData = do.call(rbind.data.frame, lapply(names(clustrows), function(x){
  ids = rownames(hm.slopes)[clustrows[[x]]]
  y = cannot.before.tolv[ids,]
  y$cluster = x
  return(y)
}))

# Demographics of Clusters
library(dplyr)
supplTable1 = clusteredRowsData %>% group_by(cluster) %>% 
  reframe(female = length(which(merged_gender == "female")),
          male = length(which(merged_gender == "male")),
          mayoLess = length(which(mayoClass == "Less")), 
          mayoMid = length(which(mayoClass == "Mid")),
          mayoMore = length(which(mayoClass == "More")), 
          ageRange = paste0(round(median(merged_age),1)," [", 
                            round(IQR(merged_age),1), "]"),
          eGFRrange = paste0(round(median(merged_eGFR),1)," [",
                             round(IQR(merged_eGFR),1), "]"),
          slopeRange = paste0(round(median(before),1)," [", 
                              round(IQR(before),1), "]"))

# Are age or sex equally distributed within clusters
library(tidyverse)
library(rstatix)
egfrsupp <- clusteredRowsData %>%
  t_test(merged_eGFR ~ cluster, p.adjust.method = "bonferroni")

slopesupp <- clusteredRowsData %>%
  t_test(before ~ cluster, p.adjust.method = "bonferroni")

agesuppp <- clusteredRowsData %>%
  t_test(merged_age ~ cluster, p.adjust.method = "bonferroni")

library(xlsx)
write.xlsx(supplTable1, file = paste0(file,"TableS3.xlsx"), sheetName = "sTable1")
write.xlsx(egfrsupp, file = paste0(file,"TableS3.xlsx"), sheetName = "egfrsupp", append = TRUE)
write.xlsx(slopesupp, file = paste0(file,"TableS3.xlsx"), sheetName = "slopesupp", append = TRUE)
write.xlsx(agesuppp, file = paste0(file,"TableS3.xlsx"), sheetName = "agesuppp", append = TRUE)

divrows = model.matrix(~0+clusteredRowsData$cluster)
propTestGender = lapply(colnames(divrows), function(x) {
  tables = rbind(table(divrows[,x], clusteredRowsData$gender)[2,],
                 table(clusteredRowsData$merged_gender))
  greater = prop.test(tables, alternative = "greater")
  greaterp = greater$p.value
  greaterods = greater$estimate
  
  less = prop.test(tables, alternative = "less")
  lessp = less$p.value
  lessods = less$estimate
  
  tests = list(greater, less)
  results = data.frame("pvalg" = greaterp, "pvall" = lessp, 
                       "oddsg" = greaterods, "oddsl" = lessods)
  return(list(results, tables, tests))
})

names(propTestGender) = colnames(divrows)

gendersuppp = lapply(propTestGender, "[[", 1)
gendersupppdf = do.call(rbind.data.frame, lapply(names(gendersuppp),function(x){
  cbind.data.frame("name" = c(as.character(gsub(".*\\$","",x)), "whole cohort"),
                   "props" = rownames(gendersuppp[[x]]),gendersuppp[[x]])}))
rownames(gendersupppdf) = NULL
write.xlsx(gendersupppdf, file = paste0(file,"TableS3.xlsx"), 
           sheetName = "gendersuppp", append = TRUE)

# | # Table1 and TableS2 preparation ----
cannot.slopesF$selctmt = as.factor(cannot.slopesF$Class_ADPKD)
levels(cannot.slopesF$selctmt) = list(unknown = "unknown", 
                                      PKD1vus = "PKD1 VUS",
                                      PKD2 = "PKD2",
                                      PKD2vus = "PKD2 VUS",
                                      PKD1NT = "PKD1 non-truncating",
                                      PKD1T = "PKD1 truncating")
#subset ADPKD COHORT 
subsADHPKD = cannot.before.tolv %>% 
  reframe(n = n(),
          sexpct = round(length(which(merged_gender == "female"))*100/length(na.omit(merged_gender)),1),
          age = paste0(round(median(merged_age),1)," [",
                       round(IQR(merged_age),1), "]"),
          egfr = paste0(round(median(merged_eGFR),1)," [",
                        round(IQR(merged_eGFR),1), "]"),
          tkv = paste0(round(median(tomography__kidney_volume,na.rm = T), 1), " [",
                      round(IQR(tomography__kidney_volume,na.rm = T),1), "]")) %>% 
  t() %>%  as.data.frame() %>% 
  tibble::rownames_to_column(var = "RowNames") %>% unname()

subsmayotb = cannot.before.tolv %>% group_by(merged_mayo) %>% 
  reframe(n = n()) %>% as.data.frame() %>% unname()
subsckdtb = cannot.before.tolv %>% group_by(merged_CKD) %>% 
  reframe(n = n()) %>% as.data.frame() %>% unname()
subsADHPKD2 = cannot.before.tolv %>% 
  reframe(
    hypt_n = length(na.omit(arterial_hypertension)),
    hypt_none = length(na.omit(which(arterial_hypertension == "0"))),
    hypt_yes = length(na.omit(grep("1|2|3",arterial_hypertension))),
    
    us_n = length(na.omit(urological_symptoms)),
    us_none = length(na.omit(which(urological_symptoms == "0"))),
    us_yes = length(na.omit(grep("1|2|-1",urological_symptoms))),
    
    mut_n = length(na.omit(selectedmut)),
    mut_pkd1t = length(na.omit(which(selectedmut == "PKD1T"))),
    mut_pkd1nt = length(na.omit(which(selectedmut == "PKD1NT"))),
    mut_pkd2 = length(na.omit(which(selectedmut == "PKD2")))) %>% t() %>%
  as.data.frame() %>% tibble::rownames_to_column(var = "RowNames") %>% unname() 

table1subsetCoh = do.call(rbind.data.frame,lapply(list(subsADHPKD,subsmayotb,
                                                       subsckdtb,subsADHPKD2),
                                                function(x) {
                                                  x[,2] = as.character(x[,2])
                                                  colnames(x) = c("rowns","vars")
                                                  return(x)}))

subsebtADPKDsupplmt = cannot.before.tolv %>%
  reframe(n = n(),
          posfam = round(length(which(family__family_members_with_adpkd == "true"))*100/length((family__family_members_with_adpkd)),1),
          hypt_n = length(na.omit(arterial_hypertension)),
          hypt_none = length(na.omit(which(arterial_hypertension == "0"))),
          hypt_l35 = length(na.omit(which(arterial_hypertension == "1"))),
          hypt_g35 = length(na.omit(which(arterial_hypertension == "2"))),
          hypt_ukage = length(na.omit(which(arterial_hypertension == "3"))),
          
          us_n = length(na.omit(urological_symptoms)),
          us_none = length(na.omit(which(urological_symptoms == "0"))),
          us_l35 = length(na.omit(which(urological_symptoms == "1"))),
          us_g35 = length(na.omit(which(urological_symptoms == "2"))),
          us_uk = length(na.omit(which(urological_symptoms == "-1")))) %>% 
  t() %>% as.data.frame() %>% 
  tibble::rownames_to_column(var = "RowNames") %>% unname() 

save.image(paste0(file,"data_screening.RData"))
writeLines(capture.output(sessionInfo()), paste0(file,"sessionInfo_screening",
                                                 gsub("-", "", Sys.Date()),
                                                 ".txt"))
# savehistory(file = paste0(file,"history_screening.Rhistory"))

library(sjPlot)
model.slopes$CST3 = as.numeric(data.before.tolv["P01034",rownames(model.slopes)])

# TableS11
tab_model(lm(before~SERPINF1+GPX3+AFM+FERMT3+CFHR1+RARRES2, model.slopes),
          lm(before~SERPINF1+GPX3+AFM+FERMT3+CST3+RARRES2, model.slopes),
          lm(before~CST3+GPX3+AFM+FERMT3+CFHR1+RARRES2, model.slopes),
          lm(before~SERPINF1+CST3+AFM+FERMT3+RARRES2+CFHR1, model.slopes),
          lm(before~SERPINF1+GPX3+AFM+FERMT3+CFHR1+RARRES2+CST3, model.slopes))
