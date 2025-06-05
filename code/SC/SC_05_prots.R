# all non-ADPKD and Mayo Class 2 ADPKD are removed and 
# batch correction was performed only on ADPKD patients

options(warn = 0)
library(diann)
library(vsn)
library(openxlsx)
library(sva)
library(readxl)

# Loading proteome data ----
setwd("~/Data/slope-models-adpkd/data")
df.sc <- diann_load("sc_from_sc_24FracLib.tsv") # file is loaded

# selecting the significant ones and melting the df
protein.groups.sc <- diann_maxlfq(df.sc[df.sc$Q.Value <= 0.01 & df.sc$PG.Q.Value <= 0.01,],
                                    # significant ones are selected here
                                    group.header = "Protein.Group", #melt
                                    id.header = "Precursor.Id",
                                    quantity.header = "Precursor.Normalised")

colnames(protein.groups.sc) <- gsub(".*\\\\","",colnames(protein.groups.sc))
colnames(protein.groups.sc) <- gsub(".raw","",colnames(protein.groups.sc))

# adding pool samples
pools.comin = as.data.frame(matrix(nrow = length(grep("pool", colnames(protein.groups.sc), ignore.case = T)), ncol = ncol(merged.cannot)))
colnames(pools.comin) = colnames(merged.cannot)
rownames(pools.comin) = grep("pool", colnames(protein.groups.sc), ignore.case = T, value = T)
merged.cannot = rbind(merged.cannot, pools.comin)

merged.cannot$Zuordnung[grep("pool", rownames(merged.cannot), ignore.case = T)] = "pool"
merged.cannot$batch[grep("pool", rownames(merged.cannot), ignore.case = T)] = gsub("_.*", "", rownames(merged.cannot)[grep("pool", rownames(merged.cannot), ignore.case = T)])
# 387 309

# selecting only ADPKD
notADPKD = which(grepl("^IgA$|^AB$", merged.cannot$Zuordnung) | merged.cannot$batch == "xxx") # IDs are masked
notADPKD = c(notADPKD,which(is.na(merged.cannot$patient_id) & merged.cannot$Zuordnung == "ADPKD"),grep(paste(wbdel, collapse = "|"), merged.cannot$patient_id))# 97

merged.cannot.slopes = merged.cannot[-notADPKD,]
commonsmpy = intersect(rownames(merged.cannot.slopes),colnames(protein.groups.sc))

# selecting common samples between proteome and annotation file
protein.groups.sc2 = protein.groups.sc[,commonsmpy] # 689 288
merged.cannot.sc.slopes = merged.cannot.slopes[commonsmpy,] # 288 309

data.na.sc <- apply(protein.groups.sc2,1,function(x) length(which(is.na(x)))/length(x))
data.sc.slopes <- protein.groups.sc2[which(data.na.sc <= 0.8),] # if it has NA more than 80%, remove it

# n patient=264
paste0("npatient=", length(unique(grep("identifier", merged.cannot.sc.slopes$patient_id, value = T))))
# n observation=264
paste0("nobservation=", length((grep("identifier", merged.cannot.sc.slopes$patient_id, value = T))))
# n followup range
paste0("nfollowupRange=[", max(table(grep("identifier", merged.cannot.sc.slopes$patient_id, value = T))),"-", # 1
       min(table(grep("identifier", merged.cannot.sc.slopes$patient_id, value = T))), "]") # 1

# normalization & imputation
data.norm.slopes = justvsn(as.matrix(data.sc.slopes)) # 398 288 samples
set.seed(12345)
data.imp.slopes <- imputeProteomics(data.norm.slopes)# 295 samples

# PCA
data.imp.slopes.pca = prcomp(scale(t(data.imp.slopes)))
dataimpPCAplot = plot(data.imp.slopes.pca, merged.cannot.sc.slopes, size = 3,
                      colour = "batch", shape = "Zuordnung") + theme_bw() +
  scale_color_discrete(name = "Batch", 
                       type = c("#03045e", "#0077b6", 
                                "#00b4d8","#90e0ef")) +
  scale_shape_discrete(name = "Group")

# outlier removals & batch
data.imp.slopes = data.imp.slopes[,-which(data.imp.slopes.pca$x[,"PC2"] < -12.5)]
merged.cannot.sc.slopes = merged.cannot.sc.slopes[-which(data.imp.slopes.pca$x[,"PC2"] < -12.5),] # 294 samples
merged.cannot.sc.slopes$Zuordnung = droplevels(merged.cannot.sc.slopes$Zuordnung)
merged.cannot.sc.slopes$batch = droplevels(merged.cannot.sc.slopes$batch)

data.comb.slopes <- ComBat(data.imp.slopes,merged.cannot.sc.slopes[colnames(data.imp.slopes),"batch"])
data.comb.slopes.pca <- prcomp(scale(t(data.comb.slopes)))
datacombPCAplot = plot(data.comb.slopes.pca, merged.cannot.sc.slopes,
                       colour = "batch", shape = "Zuordnung") + theme_bw() +
  scale_color_discrete(name = "Batch") +
  scale_shape_discrete(name = "Group")

# outlier removals
data.imp.slopes2 = data.imp.slopes[,setdiff(colnames(data.imp.slopes), names(which(data.comb.slopes.pca$x[,"PC2"] > 13)))] # 291 samples
merged.cannot.sc.slopes2 = merged.cannot.sc.slopes[setdiff(rownames(merged.cannot.sc.slopes), names(which(data.comb.slopes.pca$x[,"PC2"] > 13))),]
data.comb.slopes2 <- ComBat(data.imp.slopes2,merged.cannot.sc.slopes2[colnames(data.imp.slopes2),"batch"])
data.comb.slopes.pca2 <- prcomp(scale(t(data.comb.slopes2)))
datacombPCAplot2 = plot(data.comb.slopes.pca2, merged.cannot.sc.slopes2,
                        colour = "batch", shape = "Zuordnung") + theme_bw() +
  scale_color_discrete(name = "Batch") +
  scale_shape_discrete(name = "Group")

# outlier removals
data.imp.slopes3 = data.imp.slopes2[,setdiff(colnames(data.imp.slopes2), names(which(data.comb.slopes.pca2$x[,"PC2"] < -12)))] #288 samples
merged.cannot.sc.slopes3 = merged.cannot.sc.slopes2[setdiff(rownames(merged.cannot.sc.slopes2), names(which(data.comb.slopes.pca2$x[,"PC2"] < -12))),]

# final batch effect correction
data.comb.slopes3 <- ComBat(data.imp.slopes3,merged.cannot.sc.slopes3[colnames(data.imp.slopes3),"batch"])
data.comb.slopes.pca3 <- prcomp(scale(t(data.comb.slopes3)))
datacombPCAplot3 = plot(data.comb.slopes.pca3, merged.cannot.sc.slopes3, size = 3,
                        colour = "batch", shape = "Zuordnung") + theme_bw() +
  scale_color_discrete(name = "Batch", 
                       type = c("#03045e", "#0077b6", 
                                "#00b4d8","#90e0ef")) +
  scale_shape_discrete(name = "Group")

# before and after batch effect correction - PCAs
library(ggpubr)
ggarrange(dataimpPCAplot, datacombPCAplot3, ncol = 2,  nrow = 1, 
          labels = c("A","B"), common.legend = T,  legend = "right")
file = "~/Data/slope-models-adpkd/data/results/"
ggsave(paste0(file,"FigureS10.png"),width = 15, height = 10 , units = "in")

# removing the pools
data.slopesF <- t(scale(t(data.comb.slopes3[,rownames(merged.cannot.sc.slopes3)[which(merged.cannot.sc.slopes3$Zuordnung != "pool")]])))

cannot.slopesF <- merged.cannot.sc.slopes3[colnames(data.slopesF),]
cannot.slopesF = droplevels(cannot.slopesF) # 264
cannot.slopesF[,c("CKD.Stage.A", "Protein")] = NULL

#n patient=257
paste0("npatient=", length(unique(grep("identifier", cannot.slopesF$patient_id, value = T))))
#n observation=257
paste0("nobservation=", length((grep("identifier", cannot.slopesF$patient_id, value = T))))
#n followup range
paste0("nfollowupRange=[", max(table(grep("identifier", cannot.slopesF$patient_id, value = T))),"-", # 1
       min(table(grep("identifier", cannot.slopesF$patient_id, value = T))), "]") # 1

# generating a df for protein IDs and their corresponding gene names
proteinIDtoGenes = cbind(df.sc$Protein.Ids,df.sc$Genes)
proteinIDtoGenes_v2 = as.data.frame(unique(proteinIDtoGenes))
colnames(proteinIDtoGenes_v2) <- c( "ProteinID","Genes")

# Heatmap generation - sorted according to eGFR
library(ComplexHeatmap)
library(circlize)
cannot.slopesF$merged_mayo = as.factor(cannot.slopesF$merged_mayo)
cannot.slopesF$merged_mayo = droplevels(cannot.slopesF$merged_mayo)
wholeptoreomeHM = heatmapSlopes(scaledDat = t(data.slopesF),
                                annot = cannot.slopesF,
                                eGFR = "merged_eGFR", kmR = 3, kmC = 4,
                                show_column_names = F,
                                valsOrdR = "merged_eGFR")

library(stringr)
allprotssc = unlist(str_split(getNamesPro(rownames(data.slopesF)), pattern = ";"))
allprotssc = allprotssc[-which(allprotssc == "")]

library(readxl)
uniprotkb <- read_excel("~/Data/slope-models-adpkd/data/uniprotkb_Human_AND_reviewed_true_AND_m_2023_11_16.xlsx")

# obtaining which proteins are common in screening and validation cohorts
library(limma)
uniprotkb$gene = as.character(sapply(uniprotkb$`Gene Names`, function(x) strsplit2(x, " ")[1]))
load("~/Data/slope-models-adpkd/data/data_s_ecDat.RData")
allprotsitc = sort(getNamesPro(rownames(data.s.ecDat), 
                                proteinID = uniprotkb$Entry, 
                                Genes = uniprotkb$gene))

commonprots = intersect(allprotsitc,allprotssc) # 271

commonprotsIDs = getfullIDlist(commonprots, 
                               proteinID = proteinIDtoGenes_v2$ProteinID, 
                               Genes = proteinIDtoGenes_v2$Genes) # 285
commonprotsIDsproteome = intersect(commonprotsIDs, rownames(data.slopesF)) 
# 281 (4 of them had other IDs as well which were included)

data.slopesFold = data.slopesF
data.slopesF = data.slopesF[commonprotsIDsproteome,]
