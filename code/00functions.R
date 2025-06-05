# generating star notation for p-value 
sigstar = function(pval){
  # pval: numeric, p-value
  cc = ifelse(pval < 0.001, "***",
              ifelse(pval < 0.01, "**",
                     ifelse(pval < 0.05, "*", "")))
  return(cc)
} 

# dataframe containing age and the limits of htTKV
# this will be used within the function below (findMayo)
# to assign mayo classes
library(readr)
newmayo <- as.data.frame(read_csv("~/Data/slope-models-adpkd/data/newmayo.csv"))

# calculating the mayo class
# # usage # findMayo(mayoChart = newmayo, agem = "Age", agep = 20, httkv = 200)
findMayo = function(mayoChart, agem, agep, httkv){
  # mayoChart: dataframe, newmayo variable
  # agem: string, name of the column containing age in mayoChart
  # agep: numeric, age of the patient
  # httkv: numeric, htTKV of the patient
  if (is.na(agep) | is.na(httkv) | agep < 15 | agep > 80) {
    mayo = NA
    return(mayo)}
  else{
    httkv = round(10*as.numeric(httkv),digits = 0)/10
    agep = round(as.numeric(agep),digits = 0)
    for (i in 1:length(mayoChart[,agem])) {
      if (mayoChart[i,agem] == agep) {
        if (httkv <= mayoChart[i,"X1A"]) {
          mayo = "1A"}else if (httkv <= mayoChart[i,"X1B"]) {
            mayo = "1B"}else if (httkv <= mayoChart[i,"X1C"]) {
              mayo = "1C"}else if (httkv <= mayoChart[i,"X1D"]) {
                mayo = "1D"}else{mayo = "1E"}
        return(mayo)
      }}}}



# computing eGFR from gender & age & serum creatinine
computeGFR_upd <- function(x, dob="age", cdate=NULL, 
                           creatinine="crea", gender="gender"){
  
  # x: dataframe, clinical annotation 
  # dob: string, variable name of either date of birth or age (if cdate is null)
  # cdate: string, variable name of the measurement date
  # creatinine: string, variable name of the creatinine measurement
  # gender: string, variable name of the sex 
  
  #if age variable doesnt exists 
  # there is measurement date, dob will be date of birth
  if (!is.null(cdate)) {
    # calculate age
    age <- as.numeric(difftime(x[,cdate], x[,dob], units = "weeks"))/52.2
  }else{
    # if age variable exists
    age <- x[,dob]
  }
  # calculating eGFR
  gfr <- mapply(function(crea,age,sex){
    # if either serum creatinine, age or sex is missing, do not calculate eGFR
    if (any(is.na(crea),is.na(age),is.na(sex))) 
      return(NA)
    # if all of them exist, calculate eGFR with thresholds below
    if (length(grep("f|w", sex)) > 0 & crea <= 0.7) { # female & Scr <= 0.7
      egfr <- 144*((crea/0.7)^(-0.329)) * 0.993^age
    }else if (length(grep("f|w", sex)) > 0 & crea > 0.7) { # female & Scr > 0.7
      egfr <- 144*((crea/0.7)^(-1.209)) * 0.993^age
    }else if (length(grep("f|w", sex)) == 0 & crea <= 0.9) { # male & Scr <= 0.9
      egfr <- 141*((crea/0.9)^(-0.411)) * 0.993^age
    }else if (length(grep("f|w", sex)) == 0 & crea > 0.9) { # male & Scr > 0.9
      egfr <- 141*((crea/0.9)^(-1.209)) * 0.993^age
    }
    return(egfr)
  },x[,creatinine],age,x[,gender])
  .gfr <<- data.frame(crea = x[,creatinine], age = age, gender = x[,gender])
  return(gfr) 
}

# calculate CKD stages
eGFRtoStage <- function(egfrF)
  #egfrF: numeric, eGFR
{
  calcStage <- rep(5,length(egfrF))
  calcStage[which(is.na(egfrF))] = NA
  calcStage[which(egfrF >= 15)] <- 4
  calcStage[which(egfrF >= 30)] <- 3
  calcStage[which(egfrF >= 60)] <- 2
  calcStage[which(egfrF >= 90)] <- 1
  
  return(calcStage)
}

# calculate the stages with 3A and 3B
eGFRtoStageAB <- function(egfrF)
  #egfrF: numeric, eGFR
{
  calcStage <- rep(5,length(egfrF))
  calcStage[which(is.na(egfrF))] = NA
  calcStage[which(egfrF >= 15)] <- "4"
  calcStage[which(egfrF >= 30)] <- "3b"
  calcStage[which(egfrF >= 44)] <- "3a"
  calcStage[which(egfrF >= 60)] <- "2"
  calcStage[which(egfrF >= 90)] <- "1"
  
  return(calcStage)
}

# it divides the samples into groups according to tolvaptan usage
# this function is important for patients using tolv on/off
# slope calculation will be done separately and mean will be used for 
# "during" tolvaptan usage period 
# "medication__tolvaptan" variable contains all intervention, 
# including tolvaptan, dialysis, nephrectomy and kidney transplant
tuvele = function(x){
  # x: dataframe, clinical data
  
  # extract intervention information
  tolv = x[,"medication__tolvaptan"]
  tolv[which(is.na(tolv))] = -99
  
  # assigning the first group as 0
  first = tolv[1] 
  groups = vector()
  group.no = 0
  counter = 1
  for (i in tolv) { # go through each sample from a patient
    if (i != first) { # if it does not match with the first group
      group.no = group.no + 1 # assign it to second group
    }
    first = i # change the first group to second group
    groups[counter] = group.no
    counter = counter + 1
    # and repeat until the last sample
  }
  return(groups) # return assigned groups
}

# calculating rolling mean to identify assigned groups (before, on or after)
# for calculating slope
rolling.mean = function(vect){
  # vect: numeric vector, to be calculated 
  forwd = vector()
  vect[which(is.na(vect))] = 0
  max = vect[1]
  for (i in 1:length(vect)) {
    next.ind = vect[i]
    max = max(max, next.ind)
    forwd[i] = max
  }
  return(forwd)
}

# calculating the slope
get_coef_new = function(yy){
  # yy: dataframe, clinical annotation
  require(MASS)
  
  # if any of the samples has eGFR measurement which is NA
  if (any(is.na(yy$patient__eGFR))) { 
    indna = which(is.na(yy$patient__eGFR))
    yy = yy[-indna,] # remove that sample
  }
  
  # if after the removal, the sample size of that patient is 1 or 
  # all of the intervention information (medication__tolvaptan) is NA
  if (nrow(yy) == 1 | all(is.na(yy$medication__tolvaptan))) {
    yy$slope = NA # don't calculate the slope
    yy$rmse = NA # don't calculate the rmse
  }
  else{ # if not
    
    # check if the length of unique eGFR measurement or exam date is 1 
    if (length(unique(yy$patient__eGFR)) == 1 | length(unique(yy$exam_date)) == 1) {  
      yy$slope = NA # don't calculate the slope
      yy$rmse = NA # don't calculate the rmse
    }
    else{ # if not
      coef = rlm(patient__eGFR~exam_date, 
                 data = yy, maxit = 10000) # generate rlm model
      yy$slope = coef$coefficients[2]*365.25 # calculate the yearly slope 
      yy$rmse = summary(coef)$sigma # calculate rmse
    }
  }
  return(yy)
}

# reordering function for heatmap
reorderHeatUPP = function(scaleDat2, annot2, by, decreasing=F){
  # scaleDat2: dataframe, omics data
  # annot2: dataframe, clinical annotation
  # by: string or vector of strings, name of the column containing
  #     variable(s) to be sorted 
  if (decreasing) {
    if (length(by) == 3)
    {
      scaleDat2 = scaleDat2[order(annot2[,by[1]], annot2[,by[2]],
                                  annot2[,by[3]], decreasing = T),]
    }
    
    else if (length(by) == 2) {
      scaleDat2 = scaleDat2[order(annot2[,by[1]], annot2[,by[2]], 
                                  decreasing = T),]
    }
    
    else if (length(by) == 1) {
      scaleDat2 = scaleDat2[order(annot2[,by[1]], decreasing = T),]
    }
  }
  else if (!decreasing) {
    if (length(by) == 3)
    {
      scaleDat2 = scaleDat2[order(annot2[,by[1]], annot2[,by[2]], 
                                  annot2[,by[3]]),]
    }
    
    else if (length(by) == 2) {
      scaleDat2 = scaleDat2[order(annot2[,by[1]], annot2[,by[2]]),]
    }
    
    else if (length(by) == 1) {
      scaleDat2 = scaleDat2[order(annot2[,by[1]]),]
    }
  }
  else{
    print("no")
  }
  return(scaleDat2)
}

# heatmap function
heatmapSlopes <- function(scaledDat, annot, eGFR, kmR, kmC, dist,
                          valsOrdR=NULL, disGroup=NULL, ageGroup=NULL, ...){
  # scaledDat: dataframe, omics data
  # annot: dataframe, clinical annotation
  # eGFR: string, name of the column containing eGFR 
  # kmR: numeric, the number of clusters for k-means clustering - samples
  # kmC: numeric, the number of clusters for k-means clustering - features
  # dist: string or function, for calculating distance matrix
  # valsOrdR: string or vector of strings, name of the column containing
  #           variable(s) to be sorted
  # disGroup: string, disease group to be included in the heatmap
  # ageGroup: string, age group to be included in the heatmap
  require(ComplexHeatmap)
  ht_opt$message = FALSE
  require(hopach)
  require(ComplexHeatmap)
  require(circlize)
  set.seed(1234)
  
  scaledDat = scale(scaledDat) # scaling the data
  
  # finding common samples between scaledDat and annot
  if (!missing(disGroup) & !missing(ageGroup)) { 
    annotInd = which(annot[,"ageGroup"] != ageGroup & annot[,"Zuordnung"] != disGroup)
    annot = annot[annotInd,]
    intersect = intersect(rownames(annot), rownames(scaledDat))
  }
  else if (!missing(disGroup) & missing(ageGroup)) {
    annotInd = which(annot[,"Zuordnung"] != disGroup)
    annot = annot[annotInd,]
    intersect = intersect(rownames(annot), rownames(scaledDat))
  }
  else if (missing(disGroup) & !missing(ageGroup)) {
    annotInd = which(annot[,"ageGroup"] != ageGroup)
    annot = annot[annotInd,]
    intersect = intersect(rownames(annot), rownames(scaledDat))
  }
  else if (missing(disGroup) & missing(ageGroup)) {
    intersect = intersect(rownames(annot), rownames(scaledDat))
  }
  annot = annot[intersect,]
  scaledDat = scaledDat[intersect,]
  
  # priting dimentions 
  print(dim(annot))
  print(dim(scaledDat))
  
  # if dist argument is missing, assign it to function below
  if (missing(dist)) { 
    dist = function(x) as.dist(1 - cor(t(x)))
  }
  
  # generating row annotations - samples 
  ha = rowAnnotation(
    eGFR = annot[,eGFR],
    MAYO = annot[,"merged_mayo"],
    Age = annot[,"merged_age"],
    Sex = annot[,"merged_gender"],
    
    # assigning colors
    col = list(Age = colorRamp2(c(18.7, 40, 79.1), 
                                c("#e4c7f1", "#9999EA", "#4545D9")),
               Sex = c("female" = "pink", "male" = "cornflowerblue"),
               eGFR = colorRamp2(c(12, 72.5, 133), 
                                 c("#4b2991", "#f7667c", "#edd9a3")),
               MAYO = c("1A" = "#d9ed92", "1B" = "#99d98c", "1C" = "#52b69a", 
                        "1D" = "#168aad", "1E" = "#1e6091")
    ))
  
  # generating the heatmap
  if (missing(valsOrdR)) {
    set.seed(123123)
    Heatmap(scaledDat,
            column_names_gp = gpar(cex = 0.5),
            clustering_distance_rows = dist,
            clustering_method_rows = "average",
            row_km = kmR,
            column_km = kmC,
            row_dend_reorder = T,
            column_dend_reorder = T,
            clustering_distance_columns = dist,
            clustering_method_columns = "average",
            heatmap_legend_param = list( title = "Intensity"),
            show_row_names = F,
            right_annotation = ha, ...) 
  }
  else if (!missing(valsOrdR)) { # if variable(s) to be sorted is provided
    Heatmap(scaledDat,
            column_names_gp = gpar(cex = 0.5),
            row_order = match(rownames(reorderHeatUPP(scaledDat, annot,
                                                      valsOrdR)), 
                              rownames(scaledDat)), #this part is added
            cluster_rows = F,
            row_dend_reorder = F,
            column_dend_reorder = T,
            clustering_distance_columns = dist,
            clustering_method_columns = "average",
            heatmap_legend_param = list( title = "Intensity"),
            show_row_names = F,
            right_annotation = ha, ...)
  }
}

# imputation by sampling from the 5th percentile of the data
imputeProteomics <- function(data, lower=0.05){
  # data: dataframe, data to be imputed
  # lower: numeric, assigning the percentile of the data
  sel <- which(is.na(data)) # select which ones are NA
  data[sel] <- sample(data[which(data <= quantile(data, 
                                                  probs = lower, na.rm = T))],
                      length(sel),replace = T) # impute
  return(data)
}

# getting the names of the proteins according to 
# proteinIDtoGenes_v2 variable which is generated in later scripts, 
# contains protein ID and gene names
getNamesPro = function(protID, proteinID=proteinIDtoGenes_v2$ProteinID, 
                       Genes=proteinIDtoGenes_v2$Genes){
  # protID: string or vector of strings, protein IDs to be converted
  # proteinID: vector of strings, protein IDs
  # Genes: vector of strings, gene names
  ind = match(protID,proteinID)
  return(Genes[ind])
}
# matching the exact gene name with protein ID
getfullIDlist = function(geneName, proteinID, Genes){
  # geneName: string or vector of strings, gene names 
  # proteinID: vector of strings, protein IDs
  # Genes: vector of strings, gene names
  
  # pattern for obtaining exact gene name match
  geneName = paste0("\\b",geneName,"\\b")
  
  # if there is only one match 
  if (length(geneName) < 2) {
    indx = grep(geneName, Genes, perl = T) # extract the index
  }
  else if (length(geneName) > 1) { # if there are multiple matches
    indx = grep(paste(geneName, collapse = "|"), Genes, perl = T) # extract indexes
  }
  return(proteinID[indx]) # return protein IDs
}

# matching the gene name to a protein ID
getIDspro = function(geneName, proteinID, Genes){
  # geneName: string or vector of strings, gene names
  # proteinID: vector of strings, protein IDs
  # Genes: vector of strings, gene names
  ind = match(geneName,Genes)
  return(proteinID[ind])
}

# changing the column names
colnames_change = function(x) {
  colnames(x) = paste0("&", colnames(x))
  x
}

modmat = function(dat){
  contrast = lapply(dat[, sapply(dat, is.factor), drop = FALSE],
                  function(x) colnames_change(contrasts(x)))
  dat1 = model.matrix(~ ., data = dat, contrasts.arg = contrast)
  return(dat1)
}

# calculate age category
ageGroupFact = function(ageNum){
  # ageNum: numeric, age
  calcAgeGroup = rep("Young",length(ageNum))
  calcAgeGroup[which(is.na(ageNum))] <- NA
  calcAgeGroup[which(ageNum <= 39)] <- "Young"
  calcAgeGroup[which(ageNum >= 40)] <- "MiddleAge"
  calcAgeGroup[which(ageNum >= 60)] <- "Old"
  
  calcAgeGroup <- as.factor(calcAgeGroup)
  calcAgeGroup = factor(calcAgeGroup, # assigning levels
                        levels = c("Young", "MiddleAge", "Old"))
  return(calcAgeGroup)
}

# limma numeric
# # usage # limma.numeric.resultsUp(~before, cannot.beofre.tolv, 
#                                   t(na.omit(lassodat.mm2)), pval = 0.05)
limma.numeric.resultsUp <- function(indVar, desing.data, fit.data, 
                                    pval=0.05, adjMeth="fdr", 
                                    indiv=T, rem.int=F, compp=NULL){
  # indVar: formula, 
  # design.data: dataframe, annotation data (rows=samples, columns=features),
  # fit.data: dataframe, omics data (rows=proteins, columns=samples), 
  # indiv: T/F, returning either topTable and summary or only topTable
  # rem.int: T/F, whether to remove intercept
  # compp: matrix, a matrix containing the comparisons to be checked
  require(limma)
  myList = list()
  
  # finding common samples between desing.data and fit.data
  intsct = intersect(rownames(desing.data), colnames(fit.data))
  desing.data = desing.data[intsct,]
  fit.data = fit.data[,intsct]
  
  # checking if the variables entered in indVar formula exist in desing.data
  indvV1 = strsplit(gsub("~", "", paste(indVar, collapse = "")),
                    split = " + ", fixed = T) # splitting the formula
  
  # if every variable in the formula doesn't exist in desing.data
  if (!all(indvV1[[1]] %in% colnames(desing.data))) {
    message(paste(indvV1[[1]][which(!(indvV1[[1]] %in% colnames(desing.data)))],
                  "variable is missing in the design data", "\n", sep = " "))
    stop() # stop the function and print the message above
  }
  
  # checking if there are samples to remove
  removed = unique(c(do.call(c, lapply(indvV1[[1]], function(x) which(is.na(desing.data[,x])))), # samples that has NA in it in the corresponding features
                   match(setdiff(colnames(fit.data),intsct), colnames(fit.data)))) # samples that are not matching with fit.data
  
  # if there are samples to remove
  if (length(removed) > 0) {
    
    #remove them
    desing.data = desing.data[-removed,]
    fit.data = fit.data[,-removed]
    
    # print this message
    message(paste("Removing", length(removed) , "samples", sep = " "))
  }
  else if ( length(removed) == 0) { # if not
    message(paste("No sample was removed")) # print this message
  }
  
  # checking if any variable in the formula exists in fit.data
  if (any(indvV1[[1]] %in% rownames(fit.data))) { # if yes
    message(paste("Removing features from fit.data: ", 
                  indvV1[[1]][which(indvV1[[1]] %in% rownames(fit.data))], 
                  sep = " ")) # printing the name of the variable being removed
    
    # remove them from fit.data
    fit.data = fit.data[-na.omit(match(indvV1[[1]] ,rownames(fit.data))),]
  }
  
  # to remove the intercept
  if (rem.int == T) {
    design <- model.matrix(as.formula(paste(gsub("~", "~ 0 + ", indVar), 
                                            collapse = "")), data = desing.data)
  }
  else{
    design <- model.matrix(as.formula(indVar), data = desing.data)
  }
  
  # generate lmFit
  fit <- lmFit(fit.data, design)
  
  # to check if there are comparisons to be made
  if (!missing(compp)) {
    fit2 <- contrasts.fit(fit, compp)
    fit2 <- eBayes(fit2)
  }else{
    fit2 <- eBayes(fit)
  }
  
  if (indiv) { # both individual topTables and summary is given 
    result <- decideTests(fit2, p.value = pval, adjust.method = adjMeth)
    sum <- summary(result)
    
    topTab <- topTable(fit2, p.value = pval, adjust.method = adjMeth, 
                       number = 500)
    myResults <- list("summary" = sum, "decideTest" = result, 
                      "toptable" = topTab, "fit" = fit2)
    return(myResults)
  }else{ #only topTables are saved
    # if there are multiple independent variables
    if (grepl("+", paste(indVar, collapse = ""), fixed = T)) {
      indvV = strsplit(gsub("~", "", paste(indVar, collapse = "")),
                       split = " + ", fixed = T) # separate them 
      indvV = unlist(indvV)
      length = length(indvV)
      
      for (i in 1:length) {
        name = paste(indvV[i])
        # generate list with multiple elements containing topTable
        myList[[name]] = topTable(fit2, p.value = pval, adjust.method = adjMeth,
                                  coef =  i + 1, number = 500) 
      }
    }else{ # if only one independent variable
      indvV = gsub("~", "", paste(indVar, collapse = ""))
      name = paste(indvV)
      myList[[name]] = topTable(fit2, p.value = pval, adjust.method = adjMeth,
                                coef = indvV, number = 500)
    }
    return(myList) # return the generated list
  }
  
}


modelgen.proteome = function(dat, dependent, trainmi=F, preview=T, stepwise=T,
                             features=NULL, folds=NULL, model="lr", ...){
  # dat: dataframe, data containing features
  # dependent: string, name of the dependent variable
  # trainmi: T/F, whether training should be done
  # preview: T/F, whether to see the preview of the model
  # stepwise: T/F, whether to perform stepwise selection
  # features: string or vector of strings, which features to be included
  # folds: vector of strings, rownames of the samples selected for training fold
  # model: string, only linear regression model can be selected
  require(caret)
  # list is generated
  myList = list()
  # original version of the data is saved
  dat.orig = dat
  
  # removing NAs
  if (missing(features)) { # if no feature was selected previously
    dat.na.omit = na.omit(dat) # removing NAs
    myList[["dim"]] = dim(dat.na.omit) # saving the dimension
    
    # generating the formula
    formula = paste(dependent, "~", 
                    paste(setdiff(colnames(dat.na.omit), dependent), 
                          collapse = "+"))
    formula = as.formula(formula) # saving it as formula
  }
  else{ # if features were selected previously
    dat.na.omit = na.omit(dat[,c(features,dependent)]) # removing NAs
    myList[["dim"]] = dim(dat.na.omit) #saving the dimension
    
    # generating the formula
    formula = paste(dependent, "~", 
                    paste(setdiff(colnames(dat.na.omit), dependent), 
                          collapse = "+"))
    formula = as.formula(formula) # saving it as formula
  }
  
  # saving the final formula
  myList[["formula"]] = formula
  
  # assigning folds
  if (trainmi) { # if user wants training first
    # check the correct format
    if (missing(folds)) { 
      stop("Did you forget to put folds? as rownames ...")
    }
    else if (!is.character(folds)) {
      stop("as rownames (as.character) ...")
    }
    else{ # if the format was correct
      test <- dat.na.omit[setdiff(rownames(dat.na.omit), folds),]
      train <- dat.na.omit[folds,] # generating test and train data
    }
  }
  else{ # if user doesn't want training
    test = dat.na.omit
    train = test
  }
  
  # removing any feature if there is only 1 unique value
  removing = which(lapply(train, function(x) length(unique(x))) < 2)
  myList[["REMOVED"]] = names(removing) # saving removed features
  if (length(removing) > 0) { #if there are removed features
    for (i in 1:length(removing)) {
      train = train[,-removing[i]] # removing one by one
      
      # updating the model
      formula = update(formula, paste0(".~. - ", names(removing[i])))
    }
  }
  
  # saving the final formula
  myList[["formula"]] = formula
  
  # Model Generation
  if (model == "lr") {
    lgmod = lm(formula, data = train) #generating linear regression model
    myList[["LR model"]] = lgmod # saving the model
    message("LR is generated")
    
    accuracy = summary(lgmod)$adj.r.squared # calculating accuracy
    
    if (stepwise == T) { # if stepwise selection is prefered - BIC
      stepwise.all = stepAIC(lgmod, direction = "both", k = log(nrow(train)))
      myList[["STEPWISE LR model"]] = stepwise.all # saving the stepwise model
      message("stepwise LR is generated")
      
      accuracy = summary(stepwise.all)$adj.r.squared # updating the accuracy
    }
  }
  else{ # if any other string is provided instead of model = "lr"
    stop(paste0("only lr, not ", model))
  }
  
  # Preview
  if (preview) { # if user wants a preview
    message("model accuracy: ", accuracy) # showing the accuracy
    message("Showing only the preview. The code execution will stop here.")
    return(print(sorted)) # then stopping the execution
  }else{ # or skipping the preview
    message("preview is skipped")
  }
  
  # Testing
  if (trainmi) { # if user wants training
    # predicting 
    pred.test = if (stepwise) {
      predict(stepwise.all, newdata = test)
    }else{
      predict(lgmod, newdata = test)
    }
    
    message("prediction is done") 
    
    # calculating model performance measures
    acc.test = 1 - var(test[,dependent] - pred.test) / var(test[,dependent])
    acc.test2 = 1 - ((1 - acc.test)*(nrow(test) - 1)/(nrow(test) - ncol(test)))
    rmse.test = RMSE(pred.test,test[,dependent])
   
    # saving those measures
    myList[["Accuracy.tests"]] = acc.test
    myList[["adjr2"]] = acc.test2
    myList[["rmse"]] = rmse.test
    myList[["predictions test"]] = pred.test
    
    # printing those measures
    message("r2: ", format(acc.test, digits = 3))
    message("r2.adj: ", format(acc.test2, digits = 3))
    message("test rmse: ", format(rmse.test, digits = 3))
  }
  return(myList)
}

# calculating future eGFR by MIC model
futureeGFR <- function(age, eGFR, gender, mayo, future, na = T){
  
  # function to calculate future eGFR
  get <- function(age, eGFR, gender, future, mayo1, mayo2){
    return(21.18-1.26*gender-0.26*age+0.9*eGFR+mayo1-0.23*future+0.19*gender*future-0.02*age*future+0.001*eGFR*future+mayo2*future)
  }
  
  # assigning sex
  if (grepl("w|f",gender)) {
    gender <- 1
  }else if (grepl("m",gender)) {
    gender <- 0
  }else{
    gender <- NA
  }
  
  # calculating future eGFR according to mayo classes
  pred <- switch(mayo,
                 '1A' = get(age,eGFR,gender,future,0,0),
                 '1B' = get(age,eGFR,gender,future,0.58,-1.33),
                 '1C' = get(age,eGFR,gender,future,-1.14,-2.63),
                 '1D' = get(age,eGFR,gender,future,-1.93,-3.48),
                 '1E' = get(age,eGFR,gender,future,-6.26,-4.78),
                 ifelse(na, NA, NULL))
  return(pred)
}

# calculating the time difference between two samples
calcTimeDiff = function(id, data1, data2, d1id, d2id, d1date, d2date, which="before", eGFR=NULL){
  # id: string, id of the patient
  # data1: dataframe, clinical annotation containing each sample for patients
  # data2: dataframe, clinical annotation of patients in the proteome
  # d1id: string, column name of the IDs in data1
  # d2id: string, column name of the IDs in data2
  # d1date: string, column name of the dates in data1
  # d2date: string, column name of the dates in data2
  # which: string, column name of the slopes (before, on or after any intervention)
  # eGFR: string, column name of the eGFR measurement in data1
  if (which == "before") { # to select which slope to consider
    # matching the id to its samples 
    # before any intervention and has before slope
    indx = (grepl(id, data1[,d1id]) & data1$slope == data1$before & !is.na(data1$slope))
  }else{
    indx = (grepl(id, data1[,d1id]) & !is.na(data1$slope))
  }
  
  # calculating the time difference between each sample of a patient and 
  # sample of a patient which proteomic measurement was done 
  timedif = data1[indx,as.character(d1date)] - data2[grep(id, data2[,d2id]),as.character(d2date)]
  timedif = timedif/365.25 # converting it to year
  
  # generating a list to return
  if (!missing(eGFR)) {
    # provide eGFR value with corresponding time difference and indexes
    returnList = list("timedif" = timedif, "index" = which(indx), 
                      "eGFRs" = data1[which(indx),eGFR])
  }else{
    returnList = list("timedif" = timedif, "index" = which(indx))
  }
  return(returnList)
}

# calculating future eGFR based on the slope type (proteome? or MIC?)
calcFuteGFR = function(timedif, id, idcol, data, which, predtype, ...){
  # timedif: numeric vector, time difference in years
  # id: string, id of the patient
  # idcol: string, column name of the variable containing IDs
  # data: dataframe, containing calculated slopes
  # which: string, defining which slopes to use 
  # predtype: string, column name of the variable containing predicted slopes
  id = which(data[,idcol] == id) # selecting the patient
  if (missing(which)) { # if which predicted slope will be used is missing
    stop("mayo slopes OR prot slopes") # stopping the execution
  }else if (which == "mayo slopes") { # if it is MIC slopes
    futeGFR = unlist(lapply(as.numeric(timedif), function(x){ 
      # use futureeGFR function
      futureeGFR(age = data[id,"age"], 
                 eGFR = data[id,"eGFR"],
                 gender = data[id,"gender"],
                 mayo = as.character(data[id,"mayo5"]), 
                 future = x, ...)
    }))
  }else if (which == "prot slopes") {#if it is proteome related slopes
    futeGFR = unlist(lapply(as.numeric(timedif), function(x){
      # take the eGFR measurement from the proteome sampling date and
      # sum it with (slope*time)
      data[id,"eGFR"] + data[id,predtype]*x}))
  }else{
    stop("mayo slopes OR prot slopes") # stopping the execution
  }
  return(futeGFR)
}

getrrvgo <- function(x, type = "GO:BP", orgdb = "org.Hs.eg.db", 
                     method = "Rel", threshold=0.7 , ...){
  # x: list, result of the gost function from gprofiler2 package
  # type: string, type of enrichment term 
  # orgdb: string, package for annotation
  # method: string, distance method
  # threshold: numeric, similarity threshold 
  require(rrvgo) 
  # subseting the enrichment according to type 
  got <- subset(x$result, source == type)[,c("term_id","p_value")]
  
  # calculating similarity matrix
  simMatrix <- calculateSimMatrix(got$term_id, orgdb = orgdb, 
                                  ont = gsub("GO:","",type), method = method)
  scores <- setNames(-log10(got$p_value),got$term_id)
  
  # reducing the terms according to semantic similarity
  reducedTerms <- reduceSimMatrix(simMatrix,scores,
                                  threshold = threshold,
                                  orgdb = orgdb)
  return(list(simMatrix = simMatrix, reducedTerms = reducedTerms))
}

# plotting the similarity matrix 
scp = function(simMatrix, reducedTerms, algorithm = c("pca", "umap"),
               onlyParents = FALSE, size = "score", addLabel = TRUE,
               labelSize = 3, max.overlaps=10, flip=F, legends=F, ...){
  # simMatrix: matrix, simMatrix generated above
  # reducedTerms: dataframe, reducedTerms generated above
  # algorithm: string, algorithm to use
  # onlyParents: T/F, showing only parent terms or all terms
  # size: string, column name of reducedTerms fo deciding size of the points
  # addLabel: T/F, adding the text labels or not
  # labelSize: numeric, the size of the text labels
  # max.overlaps: numeric, max overlap of the text labels
  # flip: T/F, flipping the axes of the plot or not
  # legends: T/F, adding the legends of the plot or not
  
  # checking whether the packages are installed
  if (!all(sapply(c("ggplot2", "ggrepel", "umap"), requireNamespace, 
                  quietly = TRUE))) {
    stop("Packages ggplot2, ggrepel, umap and/or its dependencies not available. ", 
         "Consider installing them before using this function.", 
         call. = FALSE)
  }
  
  # if showing only parent terms
  if (onlyParents) {
    x <- as.data.frame(table(reducedTerms$parentTerm))
    reducedTerms <- reducedTerms[reducedTerms$term == reducedTerms$parentTerm, 
    ]
    simMatrix <- simMatrix[reducedTerms$go, reducedTerms$go]
    reducedTerms[, size] <- x$Freq[match(reducedTerms$term, 
                                         x$Var1)]
  }
  
  # calculating pca or umap
  x <- switch(match.arg(algorithm), 
              pca = cmdscale(as.matrix(as.dist(1 - simMatrix)), eig = TRUE, k = 2)$points, 
              umap = umap::umap(as.matrix(as.dist(1 - simMatrix)))$layout)
  
  # generating the data for plotting 
  df <- cbind(as.data.frame(x), reducedTerms[match(rownames(x), reducedTerms$go), 
                                             c("term", "parent", "parentTerm", size)])
  
  # plotting
  p <- ggplot2::ggplot(df, ggplot2::aes(x = V1, y = V2, color = parentTerm)) + 
    ggplot2::geom_point(ggplot2::aes_string(size = size), 
                        alpha = 0.5) +
    ggplot2::scale_size_continuous(range = c(0,25)) + 
    ggplot2::scale_x_continuous(name = "") + ggplot2::scale_y_continuous(name = "") + 
    ggplot2::theme_minimal() + ggplot2::theme(axis.text.x = ggplot2::element_blank(), 
                                              axis.text.y = ggplot2::element_blank())
  if (addLabel) { # adding labels 
    p <- p + ggrepel::geom_label_repel(ggplot2::aes(label = parentTerm), 
                                       data = subset(df, parent == rownames(df)), 
                                       box.padding = grid::unit(1, "lines"), 
                                       size = labelSize, ..., 
                                       max.overlaps = max.overlaps)
  }else {
    p
  }
  
  if (legends == F) { # adding legends
    p <- p +  ggplot2::scale_color_discrete(guide = "none") + 
      ggplot2::scale_size_continuous(guide = "none",range = c(0,25))
  }
  
  if (flip) { # flipping the axes
    p <- p + coord_flip()
  }
  p
}

createData = function(protData, annotData, model, which){
  require(limma)

  # protData: dataframe, omics data
  # annotData: dataframe, clinical annotation
  # model: object, model to be used for prediction
  # which: string, a way to extract the formula form the "model"  
  
  # extracting the formula fo the model
  if (missing(which)) {
    formula = model$formula
  }else{
    formula = model$call$formula
  }
  
  # extracting features used in the models
  sel = c(gsub(" ","",strsplit2(formula, "\\+")[3,]), "before")
  
  # subsetting the omics data and clinical annotation
  protData = t(protData[,rownames(annotData)])
  colnames(protData) = getNamesPro(colnames(protData), 
                                   proteinID = uniprotkb$Entry, 
                                   Genes = uniprotkb$gene)
  
  # merging them
  validation_s = cbind.data.frame(protData, annotData)
  validation_s$merged_eGFR = validation_s$patient__eGFR
  validation_s = validation_s[,sel]
  return(list("proteome" = protData, "cannot" = annotData, 
              "gendat" = validation_s, "model" = model))
}

# getting the mass of the proteins
getMass = function(selprot, protID, genes, view=T){
  # selprot: string, name of the proteins whose mass will be calculated
  # protID: vector of strings, protein IDs 
  # genes: vector of strings, gene names 
  # view: T/F, whether to view the table of proteins and molecular weights
  require(protr) # getUniProt function to obtain seq
  require(Peptides) # mw function to calculate molecular weights from seq
  
  options(warn = 2) # to catch errors and warnings
  
  # generating an ID list similar to enrichment
  df.mass = proteinIDtoGenes_v2[match(selprot,proteinIDtoGenes_v2$ProteinID),]
  
  # splitting the protein IDs containing ;
  prots <- strsplit(df.mass$ProteinID, ";")
  prots_unls <- unlist(prots)
  df.mass <- data.frame(protID = prots_unls, 
                        geneName = rep(df.mass$Genes, lengths(prots)))
  df.mass$mass = NA # to record calculated mass
  df.mass$warning = NA # to record problems
  
  
  sapply(df.mass$protID, function(x){
    seq = try(getUniProt(x)) # trying to get the sequence
    if (is.null(attr(seq,"class"))) { # if there is no error
      mwsq = try(mw(getUniProt(x))) # trying to calculate molecular weight
      if (is.null(attr(mwsq,"class"))) { # if there is no error
        
        # formatting the weight
        df.mass$mass[df.mass$protID == x] <<- format(round(mwsq/1000, 2), 
                                                     nsmall = 2)
      }else{ # if there is a warning
        options(warn = 0) # reverting this to default
        print(x) #  printing the protein ID
        mwsq = mw(getUniProt(x)) # forcing it to calculate
        
        # formatting the weight
        df.mass$mass[df.mass$protID == x] <<- format(round(mwsq/1000, 2), 
                                                     nsmall = 2)
        options(warn = 2) # reverting this back 
        
        # recording it as warning
        df.mass$warning[df.mass$protID == x] <<- "warning" 
      }
    }else if (attr(seq,"class") == "try-error") { # if there is an error
      
      # recording it as error
      df.mass$warning[df.mass$protID == x] <<- "error"
    }
  })
  
  # renaming the column names and removing the row names
  colnames(df.mass) = c("Protein ID", "Gene Name","Mass (kDa)", "WoE")
  rownames(df.mass) = NULL
  
  if (view == T) { # if user wants as a table
    
    # generate the table
    library(kableExtra)
    df.mass %>% kbl(align = "cc") %>%
      kable_styling(bootstrap_options = c("striped", "hover", "condensed"), 
                    full_width = F)  %>%
      column_spec(3, color = ifelse(is.na(df.mass$`Mass (kDa)`), "gray",
                                    ifelse(df.mass$`Mass (kDa)` > 60, "#d40000",
                                           ifelse(df.mass$`Mass (kDa)` > 45, "#ff7b00",
                                                  ifelse(df.mass$`Mass (kDa)` > 15,
                                                         "#ffd000", "mediumseagreen"))))) %>%
      column_spec(c(1,3), bold = T)
  }
  options(warn = 0) # reverting this to default
  return(df.mass)
}

# Plotting Related Functions ----
require(ggplot2)

plot.prcomp <- function(x, y, annot, show.legend = T, col = NULL,
                        legend.pos = "topright", ...){
  if (missing(annot) | missing(y)) {
    tmp <- factor(rep("black",nrow(x$x)))
    levels(tmp) <- "none"
  }else{
    tmp <- factor(apply(annot[rownames(x$x),y,drop = F],1,
                        paste, collapse = " - "))
  }
  if (!is.null(col) & length(col) == nrow(x$x)) {
    cold <- col
    warning("Turning off legend as no grouping available. (Check length of colours).")
    show.legend <- F
  }else if (!is.null(col) & length(col) < nrow(x$x)) {
    cold <- col[as.numeric(tmp)]
    fill.col <- col[1:length(levels(tmp))]
    leg <- levels(tmp)
  }else{
    cold <- as.numeric(tmp)
    fill.col <- 1:max(as.numeric(tmp))
    leg <- levels(tmp)
  }
  plot(x$x, col = cold, ...)
  
  if (show.legend)
    legend(legend.pos, fill = fill.col, legend = leg)
}

plot.prcomp <- function(pca, annot, pca.x = "PC1", pca.y = "PC2",
                        f.palette = "BrBG", slot="x", verbose = F, ...){
  varsExp = round(pca$sdev^2/sum(pca$sdev^2)*100, digits = 2)
  if (!missing(annot)) {
    if (!is.null(dim(annot))) {
      fd <- data.frame(annot[rownames(pca[[slot]]),],pca[[slot]])
    }else{
      fd <- data.frame(annot,pca[[slot]])
    }
  }else{
    fd <- data.frame(pca[[slot]])
  }
  if (any(colnames(fd) == "")) {
    sel <- which(colnames(fd) == "")
    colnames(fd)[sel] <- paste("PC", 1:length(sel), sep = "")
  }
  if (is.null(colnames(fd))) {
    colnames(fd) <- paste("PC", 1:ncol(fd), sep = "")
  }
  if (verbose)
    print("Setting ggplot")
  p <- ggplot(fd,aes_string(x = pca.x, y = pca.y))
  prms <<- list(...)
  prms_rem <- prms[setdiff(names(prms),c("fill","shape","colour","shape.name",
                                         "colour.name","fill.name"))]
  if (any(names(prms) == "size")) {
    if (!is.numeric(prms$size))
      prms_rem <- prms[setdiff(names(prms),c("fill","shape","colour",
                                             "shape.name","colour.name",
                                             "fill.name","size"))]
  }
  
  if (verbose) {
    print("Additional slots")
    print(prms_rem)
  }
  if (missing(annot)) {
    if (verbose)
      print("Adding points")
    p <- p + do.call(geom_point,c(prms_rem))
  }else if (is.vector(annot) & length(annot) == nrow(pca[[slot]])) {
    if (verbose)
      print("Changing shape")
    if (any(names(prms) == "shape")) {
      p <- p + do.call(geom_point,c(list(aes_string(shape = 'factor(annot)')),
                                    prms_rem))
    }else{
      p <- p + do.call(geom_point,c(list(aes_string(colour = 'annot')),
                                    prms_rem))
    }
  }else{
    if (verbose)
      print("Changing fill")
    if (any(names(prms) == "fill")) {
      if (any(names(prms) == "fill.name")) {
        f.name <- prms$fill.name
      }else{
        f.name <- prms$fill
      }
      print("Setting shape to not solids...")
      p <- p + scale_shape_manual(values = 21:25) + 
        scale_fill_brewer(name = f.name, palette = f.palette, 
                          guide = guide_legend(override.aes = aes(shape = 21)))
    }
    dl <- list()
    if (any(names(prms) == "shape")) {
      dl <- list(shape = prms$shape)
    }
    if (any(names(prms) == "colour")) {
      dl <- c(dl,list(colour = prms$colour))
    }
    if (any(names(prms) == "fill")) {
      dl <- c(dl,list(fill = prms$fill))
    }
    if (any(names(prms) == "size")) {
      dl <- c(dl,list(size = prms$size))
    }
    p <- p + do.call(geom_point,c(list(do.call(aes_string,dl)),prms_rem))
    if (any(names(prms) == "shape.name")) {
      p <- p + scale_shape_discrete(name = prms$shape.name)
    }
    if (any(names(prms) == "colour.name")) {
      if (any(names(p$labels) == "colour")) {
        p$labels$colour <- prms$colour.name
      }else{
        warning("No colour scale was identified.")
      }
    }
  }
  p + xlab(paste0("PC1: ",varsExp[1], "%")) + 
    ylab(paste0("PC2: ", varsExp[2], "%"))
}

bi = function(p, pca, pca.x = "X1", pca.y = "X2", slot.scores = "scores",
              slot.loadings = "loadings", slot.names = "center", 
              threshold = 0.75, threshold.sel = NULL, round = 2){
  data <- data.frame(obsnames = rownames(pca[[slot.scores]]), 
                     pca[[slot.scores]])
  varsExp = round(pca$sdev^2/sum(pca$sdev^2)*100, digits = round)
  
  p <- p + geom_hline(aes(yintercept = 0), size = .2) + 
    geom_vline(aes(xintercept = 0), size = .2) +
    xlab(paste0("PC1: ",varsExp[1], "%")) + 
    ylab(paste0("PC2: ", varsExp[2], "%"))
  datapc <- data.frame(varnames = names(pca[[slot.names]]), 
                       pca[[slot.loadings]])
  mult <- min(
    (max(data[,pca.y]) - min(data[,pca.y])/(max(datapc[,pca.y]) - min(datapc[,pca.y]))),
    (max(data[,pca.x]) - min(data[,pca.x])/(max(datapc[,pca.x]) - min(datapc[,pca.x])))
  )
  datapc <- transform(datapc, v1 = .7 * mult * (get(pca.x)), 
                      v2 = .7 * mult * (get(pca.y))
  )
  q1up <- quantile(datapc$v1[datapc$v1 > 0],probs = threshold)
  q1lo <- quantile(datapc$v1[datapc$v1 < 0],probs = threshold)
  q2up <- quantile(datapc$v2[datapc$v2 > 0],probs = threshold)
  q2lo <- quantile(datapc$v2[datapc$v2 < 0],probs = threshold)
  
  selv1 <- which(datapc$v1 >= q1up | datapc$v1 <= q1lo)
  selv2 <- which(datapc$v2 >= q2up | datapc$v2 <= q2lo)
  .dd <<- list(selv1, selv2, datapc, q1up, q1lo, q2up, q2lo)
  
  if (!is.null(threshold.sel)) {
    selv1 <- c(selv1[order(datapc$v1[selv1],decreasing = F)[1:floor(threshold.sel/2)]],selv1[order(datapc$v1[selv1],decreasing = T)[1:floor(threshold.sel/2)]])
    selv2 <- c(selv2[order(datapc$v2[selv2],decreasing = F)[1:floor(threshold.sel/2)]],selv2[order(datapc$v2[selv2],decreasing = T)[1:floor(threshold.sel/2)]])
  }
  sel <- unique(c(selv1,selv2))
  tmp <<- datapc
  p <- p +
    geom_text_repel(data = datapc[sel,], aes(x = v1, y = v2, label = varnames), 
                    size = 5, vjust = 1, color = "red", max.overlaps = 20) +
    geom_segment(data = datapc[sel,], aes(x = 0, y = 0, xend = v1, yend = v2), 
                 arrow = arrow(length = unit(0.2,"cm")), alpha = 0.75, 
                 color = "red")
  p
}

plot.spca = function(pca, annot, pca.x = "X1", pca.y = "X2", f.palette = "BrBG",
                     slot = "scores", verbose = F, biplot = F, 
                     slot.loadings = "loadings", slot.names = "center",
                     threshold = 0.75, threshold.sel = 20, round = 2, ...){
  p <- plot.prcomp(pca, annot, pca.x = pca.x, pca.y = pca.y, 
                   f.palette = f.palette, slot = slot, verbose = verbose, ...)
  if (biplot)
    p <- bi(p, pca, pca.x, pca.y, slot.scores = slot, 
            slot.loadings = slot.loadings, slot.names = slot.names, 
            threshold = threshold, threshold.sel = threshold.sel, round)
  p
}

runspca <- function(data,...){
  require(sparsepca)
  data.m <- data
  if (any(is.na(data.m)))
    data.m[which(is.na(data.m))] <- 0
  pca <- spca(t(data.m), ...)
  return(pca)
}

summary.spca <- function(object, ...){
  chkDots(...)
  vars <- object$sdev^2
  vars <- vars/sum(vars)
  importance <- rbind("Standard deviation" = object$sdev,
                      "Proportion of Variance" = round(vars, 5),
                      "Cumulative Proportion" = round(cumsum(vars), 5))
  colnames(importance) <- colnames(object$scores)
  object$importance <- importance
  class(object) <- "summary.spca"
  object
}

print.summary.spca <- function(x, digits = max(3L, getOption("digits") - 3L), 
                               ...){
    cat("Importance of components:\n")
    print(x$importance, digits = digits, ...)
    invisible(x)
  }