options(warn = 0)
library(diann) 
library(vsn) 
library(openxlsx) 
library(sva)  
library(readxl)
source("~/Data/slope-models-adpkd/code/00functions.R", echo = TRUE)

# Clinical Data -----
library(openxlsx)
cannot <- read.xlsx("~/Data/slope-models-adpkd/data/ADPKDComplete_v2.xlsx")

# removing some variables
cannot[,c("Monat", "CKD.Stage.G", "Control.vs..Affected", "eGFR", 
          "Alter..in.Jahren.", "Jahr")] = NULL

# patients who don't have identifier ID are removed
cannot = cannot[-which(cannot$identifier == cannot$Labornummer.Serum),]

# preparing the annotation data
cannot$Zuordnung <- gsub("Kontrolle","AB",cannot$Zuordnung) 
cannot$Labornummer.Serum = as.character(cannot$Labornummer.Serum)
cannot$ID.sci <- gsub("xxx","xxx",cannot$ID.sci) # IDs are masked
cannot$Geschlecht = gsub('\\s+', '',cannot$Geschlecht)
cannot$eGFR_v2 = as.numeric(cannot$eGFR_v2)

# converting dates
cannot$Einschluss[grep("10", nchar(cannot$Einschluss))] = as.character(as.Date(cannot$Einschluss[grep("10", nchar(cannot$Einschluss))], format = '%d/%m/%Y'))
cannot$Einschluss[grep("5", nchar(cannot$Einschluss))] = as.character(as.Date(as.numeric(cannot$Einschluss[grep("5", nchar(cannot$Einschluss))]), origin = "1899-12-30"))

cannot$Geburtsdatum = as.Date(cannot$Geburtsdatum, format = '%d/%m/%Y')
cannot$Datum.Untersuchung = as.Date(cannot$Datum.Untersuchung, format = '%d/%m/%Y')
cannot$Einschluss = as.Date(cannot$Einschluss)

cannot$CKD = eGFRtoStageAB(cannot$eGFR_v2)
cannot$CKD_v2 = NULL

cannot$Mayo.class = gsub("^-$", NA, cannot$Mayo.class)

wbfactor = names(which(unlist(lapply(cannot, function(x) length(unique(na.omit(x))))) < 15))
cannot[,wbfactor] = lapply(cannot[,wbfactor], factor)

# removing the ones that don't have batch number
cannot = cannot[which(!is.na(cannot$ID.sci)),]
rownames(cannot) = cannot$ID.sci
cannot$batch = gsub("_.*","",cannot$ID.sci)

# Proteome Data -----
setwd("~/Data/slope-models-adpkd/data/")
df.sci <- diann_load("sc_from_sc_24FracLib.tsv")

# selecting the significant ones and melting the df
protein.groups.sci <- diann_maxlfq(df.sci[df.sci$Q.Value <= 0.01 & df.sci$PG.Q.Value <= 0.01,],
                                    #significant ones are selected here
                                    group.header = "Protein.Group", #melt
                                    id.header = "Precursor.Id",
                                    quantity.header = "Precursor.Normalised")

colnames(protein.groups.sci) <- gsub(".*\\\\","",colnames(protein.groups.sci))
colnames(protein.groups.sci) <- gsub(".raw","",colnames(protein.groups.sci))

# adding pool samples
pools.comin = as.data.frame(matrix(nrow = length(grep("pool", colnames(protein.groups.sci), 
                                                      ignore.case = T)), 
                                   ncol = ncol(cannot)))
colnames(pools.comin) = colnames(cannot)
rownames(pools.comin) = grep("pool", colnames(protein.groups.sci), 
                             ignore.case = T, value = T)
cannot = rbind(cannot, pools.comin)

cannot$Zuordnung[grep("pool", rownames(cannot), ignore.case = T)] = "pool"
cannot$batch[grep("pool", rownames(cannot), ignore.case = T)] = gsub("_.*", "", rownames(cannot)[grep("pool", rownames(cannot), ignore.case = T)])

# selecting only IgAN samples
onlyIgan = grep("xxxxx",cannot$batch) #36 - 6 pools + 30 IgAN samples # IDs are masked
cannot = cannot[onlyIgan,]
commonsmpy = intersect(rownames(cannot), colnames(protein.groups.sci))

# selecting common samples between proteome and annotation file
protein.groups.sci2 = protein.groups.sci[,commonsmpy] #  689  36
cannot = cannot[commonsmpy,] # 36 20

data.na.sci <- apply(protein.groups.sci2,1,function(x) length(which(is.na(x)))/length(x))
data.s.sci <- protein.groups.sci2[which(data.na.sci <= 0.8),] # if it has NA more than 80%, remove it

# normalization & imputation
data.norm = justvsn(as.matrix(data.s.sci))
set.seed(12345)
data.imp <- imputeProteomics(data.norm)

# PCA
data.imp.ca = prcomp(t(data.imp), scale. = T)
plot(data.imp.ca, cannot, colour = "batch", shape = "Zuordnung")

# outlier removals
data.imp = data.imp[,-which(data.imp.ca$x[,"PC1"] > 40)]
cannot = cannot[-which(data.imp.ca$x[,"PC1"] > 40),]

# PCA
data.imp.ca = prcomp(t(data.imp), scale. = T)
plot(data.imp.ca, cannot, colour = "batch", shape = "Zuordnung")

# removing the pools
cannot$Zuordnung = droplevels(cannot$Zuordnung)
data.final <- data.imp[,rownames(cannot)[which(cannot$Zuordnung != "pool")]]
cannot.final <- cannot[colnames(data.final),]

proteinIDtoGenes = cbind(df.sci$Protein.Ids, df.sci$Genes)
proteinIDtoGenes_v2 = as.data.frame(unique(proteinIDtoGenes))
colnames(proteinIDtoGenes_v2) <- c( "ProteinID","Genes")
cannot.final = droplevels(cannot.final)

# merging proteome and annotation data
cannot.final[,c("CKD.Stage.A", "Protein")] = NULL
protannot = cbind.data.frame(t(data.final), cannot.final)

characters = which(unlist(lapply(protannot, class)) == "character")
dates = which(unlist(lapply(protannot, class)) == "Date")
protannot = protannot[,-c(characters, dates)]
dim(protannot) # 29 389

library(readr)
all29protsel <- read_csv("~/Data/slope-models-adpkd/data/all29protsel.csv")

# Correlation of 29 proteins ----
library(makeunique)
igan_data = protannot[which(cannot.final$Zuordnung == "IgA"),all29protsel$proteinID]
colnames(igan_data) = all29protsel$geneName
igan_data_all = cbind.data.frame(igan_data, "eGFR" = cannot.final[which(cannot.final$Zuordnung == "IgA"),c("eGFR_v2")])

library(Hmisc)
allcompcorIGAN <- rcorr(as.matrix(igan_data_all))
cor_all_igan <- allcompcorIGAN$r
p_all_igan <- allcompcorIGAN$P

allcompcorIGAN = data.frame("corGFR" = cor_all_igan[,"eGFR"], 
                            "pGFR" = p_all_igan[,"eGFR"],
                            "sigstarGFR" = sigstar(p_all_igan[,"eGFR"]))
allcompcorIGAN
write.csv(allcompcorIGAN, file = "~/Data/slope-models-adpkd/data/results/TableS5.csv")

# Characteristics of IgAN Cohort ----
library(dplyr)
# clinical characteristics
igan_clin = cannot.final %>% filter(Zuordnung == "IgA") %>% 
  reframe(age, eGFR_v2, Geschlecht)
# shapiro.test(igan_clin$age)
# shapiro.test(igan_clin$eGFR_v2)
igan_clin_ch = igan_clin %>% 
  reframe(n = n(),
          sex = round(length(which(Geschlecht == "w"))*100/length(na.omit(Geschlecht)), 1),
          age = paste0(round(median(age),1), " (", round(IQR(age), 1), ")"),
          eGFR = paste0(round(median(eGFR_v2), 1), 
                        " (", round(IQR(eGFR_v2), 1), ")")) %>% t()

library(ggpubr)
stable.p <- ggtexttable(igan_clin_ch, 
                        theme = ttheme("classic", 
                                       rownames.style = rownames_style(face = "plain",linecolor = "black")),
                        rows = c("Patients, n","Sex (Female %)", "Age (years), median (IQR)", 
                                 bquote("eGFR (ml/min/" ~ m^2 ~ "), median (IQR)")),
                        cols = "IgAN Cohort")
stable.p


# PCA
library(circlize)
set.seed(209902)
indigan = cannot.final$Zuordnung == "IgA"

library(ggrepel)
midd = mean(cannot.final[indigan,"eGFR_v2"], na.rm = T)
pcaFig1 = plot.spca(prcomp(t(data.final[,indigan])), cannot.final[indigan,],
                    colour = "eGFR_v2", shape = "Geschlecht",
                    size = 2, pca.x = "PC1", pca.y = "PC2", slot = "x",
                    slot.loadings = "rotation", slot.names = "center") +
  theme_bw() + theme(legend.title = element_text(face = "bold")) +
  scale_color_gradient2(midpoint = midd, low = "#4b2991", mid = "#f7667c",
                        high = "#edd9a3", space = "Lab", name = "eGFR") +
  scale_shape_manual(name = "Sex", labels = c("F", "M"), values = c(16, 17))
pcaFig1

# Heatmap
library(ComplexHeatmap)

scaledDat = t(data.final[,indigan])
annot = cannot.final[indigan,]
eGFR = "eGFR_v2"
kmR = 3
kmC = 4
meth = "average"
show_column_names = F
valsOrdR = "eGFR_v2"

scaledDat = scale(scaledDat)
intersect = intersect(rownames(annot), rownames(scaledDat))

annot = annot[intersect,]
scaledDat = scaledDat[intersect,]

dist = function(x) as.dist(1 - cor(t(x)))
ha = rowAnnotation(
  eGFR = annot[,"eGFR_v2"],
  Age = annot[,"age"],
  Sex = annot[,"Geschlecht"],
  
  col = list(Age = colorRamp2(c(18.7, 40, 79.1), 
                              c("#e4c7f1", "#9999EA", "#4545D9")),
             Sex = c("w" = "pink", "m" = "cornflowerblue"),
             eGFR = colorRamp2(c(12, 72.5, 133), 
                               c("#4b2991", "#f7667c", "#edd9a3"))
  ))

igan_hm = Heatmap(scaledDat, row_names_gp = gpar(cex = 0.3), 
                  column_names_gp = gpar(cex = 0.5),
                  row_order = match(rownames(reorderHeatUPP(scaledDat, annot, 
                                                            valsOrdR)), 
                                    rownames(scaledDat)),
                  cluster_rows = F,
                  row_dend_reorder = F,
                  column_dend_reorder = T,
                  clustering_distance_columns = dist,
                  clustering_method_columns = "average",
                  #col = colorRamp2(c(-3, 0, 3), c("green", "black", "red")),
                  heatmap_legend_param = list( title = "Intensity"),
                  show_row_names = F, show_column_names = F,
                  right_annotation = ha)

ggarrange(ggarrange(stable.p, pcaFig1, nrow = 2,labels = c("A", "B"), 
                    heights = c(0.4, 1)), 
          ggarrange(NULL, grid.grabExpr(draw(igan_hm)), ncol = 2, 
                    widths = c(0.025, 1)), labels = c("", "C"), ncol = 2)

ggsave(file = "~/Data/slope-models-adpkd/data/results/FigureS4.png",
       width = 12, height = 6 , units = "in")


