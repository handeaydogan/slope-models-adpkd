source("~/Data/slope-models-adpkd/code/00functions.R", echo = TRUE)
source("~/Data/slope-models-adpkd/code/SC/SC_01_annot.half.R", echo = TRUE)
source("~/Data/slope-models-adpkd/code/SC/SC_02_interventions.R", echo = TRUE)
source("~/Data/slope-models-adpkd/code/SC/SC_03_mergedCreasUpd.R", echo = TRUE)
source("~/Data/slope-models-adpkd/code/SC/SC_04_annots_merging.R", echo = TRUE)
source("~/Data/slope-models-adpkd/code/SC/SC_05_prots.R", echo = TRUE)

# Get Size Distribution of ALL proteins ----

# # generating the mass data of all proteins (protMassAll) and 
# # detected proteins (protMasssc)
# protMasssc = getMass(rownames(data.slopesFold), 
#                        proteinIDtoGenes_v2$ProteinID, 
#                        proteinIDtoGenes_v2$Genes)
# protMasssc = na.omit(protMasssc[,1:3])
# protMasssc$`Mass (kDa)` = as.numeric(protMasssc$`Mass (kDa)`) # 534   3
# 
# protMassAll = readFASTA("~/Data/slope-models-adpkd/data/uniprotkb_taxonomy_id_9606_2024_02_01.fasta")
# protMassAll = mw(protMassAll)
# protMassAll = format(round(protMassAll/1000, 2), nsmall = 2)
# protMassAll = as.data.frame(as.numeric(protMassAll)) # 204141 1
# colnames(protMassAll) = "allprots"
# # save(protMasssc, protMassAll, 
# #      file = paste0("~/Data/slope-models-adpkd/data/massDistributionProts",
# #                    gsub("-", "", Sys.Date()), ".RData"))

load("~/Data/slope-models-adpkd/data/massDistributionProts20240726.RData")

# density plots 
denistyMass = ggplot() +
  geom_density(data = protMassAll, aes(x = allprots, colour = 'gray40'), 
               fill = 'gray60', alpha = 0.3) +
  geom_density(data = protMasssc, aes(x = `Mass (kDa)`, colour = '#11A64A'), 
               fill = 'green', alpha = 0.3) +
  scale_colour_manual('Proteins', limits = c("ADPKD Serum \n(Proteome)", 
                                             "Entire Human  \nProteins (UniProtKB)"), 
                      values = c("green", "gray60")) +
  guides(colour = guide_legend(override.aes = list(fill = c("green", "gray60"), alpha = 0.3))) +
  scale_colour_manual("Proteins", values = c("#11A64A", "gray40"),
                      labels = c("ADPKD Serum \n(Proteome)",
                                 "Entire Human  \nProteins (UniProtKB)")) +
  scale_x_log10() + theme_bw() + 
  xlab(expression("log"["10"]~"Mass (kDa)")) + ylab("Density") +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold"))

# Get Enriched Pathways and GO Terms ---- 
# | #  Preprocessed proteins (from data.slopesFold) ----
library(gprofiler2)
# generating an ID list similar to enrichment
keyenrichALL = proteinIDtoGenes_v2[match(rownames(data.slopesFold),
                                         proteinIDtoGenes_v2$ProteinID),] # 398 2

# splitting the protein IDs containing ;
split_pro_all_enr <- strsplit(keyenrichALL$ProteinID, ";")
split_pro_all_enr_unls <- unlist(split_pro_all_enr)
keyenrichALL <- data.frame(protID = split_pro_all_enr_unls, 
                           geneName = rep(keyenrichALL$Genes, 
                                          lengths(split_pro_all_enr))) # 566 2

enrich_all_ensgID = gconvert(keyenrichALL$protID, 
                             organism = "hsapiens", 
                             target = "ENSG")  # 603   7
enrich_all_ensgID$matched = sapply(enrich_all_ensgID$input, function(x) keyenrichALL$geneName[grep(paste0("\\b",x,"\\b"), keyenrichALL$protID)])
enrich_all_ensgID$tf = NA
enrich_all_ensgID$tf = mapply(function(name, matched){# creating a true false variable - tf
  # find a match between the name came up from converter and 
  # the name from the proteome data containing ";"
  grpfnc = grep(paste0("\\b",name,"\\b"), matched) #this check for matching whole words
  ifelse(length(grpfnc) > 0, T, F) #if there is a match, assign tf to true or false
},name = enrich_all_ensgID$name, enrich_all_ensgID$matched) # 603   9

#only novels contain first occurrence and they are only occurrence
none_novels = which(grepl("^None$", enrich_all_ensgID$name, ignore.case = T) & grepl("novel transcript|novel gene|novel protein",enrich_all_ensgID$description))  #28
linc_orf_as = which(grepl("^LINC", enrich_all_ensgID$name) | grepl("open reading frame|antisense",enrich_all_ensgID$description)) #24

enrich_all_ensgIDrefined = enrich_all_ensgID[-c(unique(none_novels, linc_orf_as)),] 
# 575   9

changereq_df = enrich_all_ensgIDrefined[enrich_all_ensgIDrefined$tf == T,] # 496   9
changereq = changereq_df[which(gsub(".*\\.", "", changereq_df$target_number) != "1"),"input"] # 3

wbchange_enr = grep(paste(changereq, collapse = "|"), enrich_all_ensgIDrefined$input) # 6
wbmerged_enr = enrich_all_ensgIDrefined[wbchange_enr,] # getting those 6 #6 9
enrich_all_ensgIDrefined = enrich_all_ensgIDrefined[-wbchange_enr,] # 569 9 
# and removing them from data

# changing those 3 proteins to get additional matches 
wbmerged_enr$tf = c(F,T,T,T,F,T)
wbmerged_enr = wbmerged_enr[wbmerged_enr$tf == T,] #4 9 
# and remove the false matches #4 matches stayed 

enrich_all_ensgIDrefined = enrich_all_ensgIDrefined[grep("\\.1", enrich_all_ensgIDrefined$target_number),] #523   9 #getting the first hits
enrich_all_ensgIDrefined = rbind.data.frame(enrich_all_ensgIDrefined, wbmerged_enr) 
#527 9 # merging it with the 3 protein hits

# | #  Enrichment ----
# setting archived url for gprofiler
set_base_url("https://biit.cs.ut.ee/gprofiler_archive3/e111_eg58_p18")

# performing enrichment
set.seed(123)
all_enrich555_new = gost(query = unique(enrich_all_ensgIDrefined$target),
                         organism = "hsapiens", ordered_query = F,
                         user_threshold = 0.001, correction_method = "fdr",
                         evcodes = T)

all_enrich555_new_res = as.data.frame(all_enrich555_new$result)
all_enrich555_new_res = apply(all_enrich555_new_res, 2, as.character)
file = "~/Data/slope-models-adpkd/data/results/"
write.csv(x = all_enrich555_new_res,
          file = paste0(file, "DataS1.csv"))

# | #  Parent Enriched GO:BP terms ----
all_enrich555_new.res.sim = getrrvgo(all_enrich555_new)
sort(table(all_enrich555_new.res.sim$reducedTerms$parentTerm),decreasing = T) 

gotermsparent = unique(cbind.data.frame(all_enrich555_new.res.sim$reducedTerms$parent,
                                        all_enrich555_new.res.sim$reducedTerms$parentTerm))
write.csv(x = gotermsparent, file = paste0(file, "DataS2.csv"))

# plotting parent enriched terms
set.seed(029510)
enrichPAfig1 = scp(all_enrich555_new.res.sim$simMatrix,
                   all_enrich555_new.res.sim$reducedTerms,
                   max.overlaps = 42, algorithm = "pca", 
                   addLabel = T, onlyParents = F, flip = F,
                   max.time = 5, force = 15, #force_pull = 3, 
                   max.iter = 1000000000, segment.curvature = -1e-20,
                   arrow = arrow(length = unit(0.005, "npc"), type = "closed"),
                   segment.size = 0.5)
enrichPAfig1 = enrichPAfig1 + coord_flip(clip = "off", ylim = c(-0.4, 0.6))

# changing colors
colnmbr = length(unique(all_enrich555_new.res.sim$reducedTerms$parentTerm))
colnmbr2 = colnmbr + 8

plasma_n = hcl.colors(colnmbr2, palette = "plasma")[1:colnmbr]
enrichPAfig1 = enrichPAfig1 + 
  scale_color_manual(values = plasma_n, guide = "none")

enrichPAfig1
ggsave(paste0(file,"FigureS3.png"),width = 12, height = 15 , units = "in")
# PCA ----
library(ComplexHeatmap)
library(circlize)
set.seed(209902)
indadpkd = cannot.slopesF$Zuordnung == "ADPKD"

library(ggrepel)
midd = mean(cannot.slopesF[indadpkd,"merged_eGFR"], na.rm = T)
pcaFig1 = plot.spca(prcomp(t(data.slopesFold[,indadpkd])),
                    cannot.slopesF[indadpkd,],
                    colour = "merged_eGFR", shape = "merged_gender",
                    size = 2, pca.x = "PC1", pca.y = "PC2", slot = "x",
                    slot.loadings = "rotation", slot.names = "center") +
  theme_bw() + theme(legend.title = element_text(face = "bold")) +
  scale_color_gradient2(midpoint = midd, low = "#4b2991", mid = "#f7667c",
                        high = "#edd9a3", space = "Lab", name = "eGFR") +
  scale_shape_manual(name = "Sex", labels = c("F", "M"), values = c(16,17))
pcaFig1

library(ggpubr)
ggarrange(ggarrange(denistyMass, pcaFig1, nrow = 2,labels = c("A", "B")), 
          ggarrange(NULL,grid.grabExpr(draw(wholeptoreomeHM)), 
                    ncol = 2, widths = c(0.025, 1)), labels = c("","C"), 
          ncol = 2, widths = c(1.3,2))

ggsave(paste0(file,"FigureS2.png"),width = 12, height = 8 , units = "in")
