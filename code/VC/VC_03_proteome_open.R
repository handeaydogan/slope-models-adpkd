
#n patient
paste0("npatientitc=", length(unique(grep("identifier", metaclin2$identifiers, value = T)))) # 465
paste0("npatientec=", length(unique(grep("xxx", metaclin2$identifiers, value = T)))) # 221
#n observation
paste0("nobservationitc=", length((grep("identifier", metaclin2$identifiers, value = T)))) # 623
paste0("nobservationec=", length((grep("xxx", metaclin2$identifiers, value = T)))) # 221
#n followup range
paste0("nfollowupRangeitc=[", max(table(grep("identifier", metaclin2$identifiers, value = T))),"-", # 3
       min(table(grep("identifier", metaclin2$identifiers, value = T))), "]") # 1

# normalization & imputation
data.norm.slopes = justvsn(as.matrix(data.s.ecDat))
set.seed(12345)
data.imp.slopes <- imputeProteomics(data.norm.slopes) # 338 844

# PCA
data.imp.slopes.pca = prcomp(scale(t(data.imp.slopes)))
dataimpPCAplotCP = plot(data.imp.slopes.pca, metaclin2, size = 3,
                        colour = "column", shape = "factor(ProjID)") + 
  theme_bw() +
  scale_color_discrete(name = "Batch", 
                       type = c("#03045e", "#00b4d8","#90e0ef")) +
  scale_shape_discrete(name = "Cohorts", 
                       labels = c("itc" = "Internal/Temporal (ITC)", 
                                  "ec" = "External (EC)"))

# outlier removals
rem = which(data.imp.slopes.pca$x[,"PC2"] > 15)
data.imp.slopes = data.imp.slopes[,setdiff(colnames(data.imp.slopes),
                                           names(rem))] # 338 839
metaclin2 = metaclin2[setdiff(rownames(metaclin2),names(rem)),] # 839  37

# selecting common samples between proteome and annotation file
metaclin2 = droplevels(metaclin2)
commonsamps = intersect(colnames(data.imp.slopes), rownames(metaclin2))
metaclin2 = metaclin2[commonsamps,]
data.imp.slopes = data.imp.slopes[,commonsamps]

# batch effect correction
data.comb.slopes <- ComBat(data.imp.slopes,metaclin2$column)
data.comb.slopes.pca <- prcomp(scale(t(data.comb.slopes)))

datacombPCAplotCPcbt = plot(data.comb.slopes.pca, metaclin2, size = 3,
                            colour = "column", shape = "factor(ProjID)") + 
  theme_bw() +
  scale_color_discrete(name = "Batch", 
                       type = c("#03045e", "#00b4d8","#90e0ef")) +
  scale_shape_discrete(name = "Cohorts", 
                       labels = c("itc" = "Internal/Temporal (ITC)", 
                                  "ec" = "External (EC)"))

library(ggpubr)
ggarrange(dataimpPCAplotCP, datacombPCAplotCPcbt, ncol = 2,  nrow = 1, 
          labels = c("A","B"), common.legend = T,  legend = "right")
file = "~/Data/slope-models-adpkd/data/results/"
ggsave(paste0(file,"FigureS11.png"),width = 15, height = 10 , units = "in")

data.slopesF <- t(scale(t(data.comb.slopes))) #338 839

metaclin2 <- metaclin2[colnames(data.slopesF),] #839  37
metaclin2 = droplevels(metaclin2)

#n patient
paste0("npatientitc=", length(unique(grep("identifier", metaclin2$identifiers, value = T)))) # 462
paste0("npatientec=", length(unique(grep("xxx", metaclin2$identifiers, value = T)))) # 219
#n observation
paste0("nobservationitc=", length((grep("identifier", metaclin2$identifiers, value = T)))) # 620
paste0("nobservationec=", length((grep("xxx", metaclin2$identifiers, value = T)))) # 219
#n followup range
paste0("nfollowupRangeitc=[", max(table(grep("identifier", metaclin2$identifiers, value = T))),"-", # 3
       min(table(grep("identifier", metaclin2$identifiers, value = T))), "]") # 1

# getting uniprotkb to match protein IDs with gene names
library(readxl)
uniprotkb <- read_excel("~/Data/slope-models-adpkd/data/uniprotkb_Human_AND_reviewed_true_AND_m_2023_11_16.xlsx")
library(limma)
uniprotkb$gene = as.character(sapply(uniprotkb$`Gene Names`, 
                                     function(x) strsplit2(x, " ")[1]))
