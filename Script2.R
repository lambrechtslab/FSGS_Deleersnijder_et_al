#### PODOCYTE SUBCLUSTERING AND ANALYSIS ####

##### Script info #####

# This script starts from loading in the podocyte (POD) subcluster, extracted from the main Seurat object in part 2.

# 1. PART 1: We re-run Seurat pipeline for the first time on the podocyte subcluster with lowQ cells still included. We perform first QC checks with first removal of lowQ cells, based on visualization of QC-parameters, doublet visualization, marker gene expression, KPMP annotation. 
# Result of part 1 is "seurat_POD_HQ"-object with annotation "POD_annot".

# 2. PART 2: We re-run Seurat pipeline for the second time after first removal of lowQ cells in part 1. We perform additional QC checks and remove additional lowQ cells.
# Result of part 2 is "seurat_POD_HQ2"-object with annotation "POD_HQ_annot".

# 3. PART 3: We re-run Seurat pipeline for a final 3rd time after removal of lowQ cells in part 2. We re-annotate cell subclusters, but do not remove additional cells
# Result of part 3 is the final podocyte-object: "seurat_POD_HQ2"-object with annotation "POD_HQ2_annot". 

# 4. PART 4: We perform DGE analysis (using EdgeR package) of one diagnostic group vs. all other three groups.
# Output is visualized using heatmaps, dotplots, volcanoplots
# Subsequently, the output from EdgeR DGE analysis is used for gene set enrichment analysis (GSEA) using the fgsea package.

# 5. PART 5: We perform DGE analysis (using EdgeR package) of one diagnostic group vs. only the two control groups.

# 6. PART 6: We perform DGE analysis (using EdgeR package) of primary FSGS vs. maladaptive FSGS.



##### 1. PART 1 - Running pipeline for first time - lowQ cells still included #####

## Setting wd and projectfolder:
setwd("path")
getwd()
projectFolder <- "path/"

# Loading packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)
library(biomaRt)
library(pheatmap)
library(RColorBrewer)
library(gridExtra)
library(DoubletFinder)
library(tidyverse)
library(cowplot)
library(future)
library(patchwork)
library(harmony)
library(celldex)
library(SeuratDisk)
library(EnhancedVolcano)
library(hypeR)

# Read in the dataset
seurat_POD <- readRDS(file = "seurat_POD.Rds")

seurat_POD
# An object of class Seurat 
# 36780 features across 3441 samples within 1 assay 
# Active assay: RNA (36780 features, 2000 variable features)
# 3 dimensional reductions calculated: pca, umap, harmony

###### 1.1: Re-Running Seurat pipeline ###### 

## Normalizing the data
seurat_POD <- NormalizeData(seurat_POD, normalization.method = "LogNormalize", scale.factor = 10000)

## Highly variable features
seurat_POD <- FindVariableFeatures(seurat_POD, selection.method = "vst", nfeatures = 2000)

## Scaling the data + regressing for mt.RNA, nCount, nFeature
seurat_POD <- ScaleData(seurat_POD, vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt")) # we only use variable features

## PCA - linear dimensional reduction
seurat_POD <- RunPCA(seurat_POD, features = VariableFeatures(object = seurat_POD), npcs = 50)

pdf(paste0(projectFolder, "seurat_POD_noHarmony_elbowPlot.pdf"), width = 10, height = 10)
ElbowPlot(seurat_POD, n = 50)
dev.off()

## Clustering
numPCs <- 18
seurat_POD <- FindNeighbors(seurat_POD, dims = 1:numPCs)
res = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1, 1.2, 1.5, 1.7, 2)
seurat_POD <- FindClusters(seurat_POD, resolution = res)
seurat_POD <- RunUMAP(seurat_POD, dims = 1:numPCs)

pdf(paste0(projectFolder, "seurat_POD_noHarmony_dimPlot_umap_multiRes_nPC",numPCs,".pdf"), width = 10, height = 8)
for (i in 1:length(res)){
  resolution <- paste0("RNA_snn_res." , res[i])
  print(DimPlot(seurat_POD, reduction = "umap", group.by = resolution, label = T, pt.size = 0.7, raster.dpi = c(1024,1024)) + ggtitle(paste0("Resolution: ", res[i])))
}
dev.off()

## Highlighting samples on the plot (i.e. plotting cells from individual samples) ##
seurat_POD <- SetIdent(seurat_POD, value = "orig.ident")

Lorigin <- list(
  snrFSGS001="snrFSGS001",
  snrFSGS002="snrFSGS002",
  snrFSGS003="snrFSGS003",
  snrFSGS004="snrFSGS004",
  snrFSGS005="snrFSGS005",
  snrFSGS006="snrFSGS006",
  snrFSGS007="snrFSGS007",
  snrFSGS008="snrFSGS008",
  snrFSGS009="snrFSGS009",
  snrFSGS010="snrFSGS010",
  snrFSGS011="snrFSGS011",
  snrFSGS012="snrFSGS012",
  snrFSGS013="snrFSGS013",
  snrFSGS014="snrFSGS014",
  snrFSGS015="snrFSGS015",
  snrFSGS016="snrFSGS016",
  snrFSGS017="snrFSGS017",
  snrFSGS018="snrFSGS018",
  snrNEPH027="snrNEPH027",
  snrNEPH028="snrNEPH028",
  snrNEPH029="snrNEPH029",
  snrNEPH030="snrNEPH030",
  snrNEPH031="snrNEPH031",
  snrNEPH032="snrNEPH032",
  snrNEPH033="snrNEPH033")

# Highlight the samples on the dimplot 
pdf(paste0(projectFolder, "seurat_POD_noHarmony_samplesHighlighted_Dimplot_nPC",numPCs,".pdf"), width = 10, height = 10)
lapply(seq_along(Lorigin),function(i) 
  DimPlot(seurat_POD, cells.highlight = WhichCells(seurat_POD, idents = Lorigin[[i]])) & ggtitle(names(Lorigin)[[i]]))
dev.off()

pdf(paste0(projectFolder, "seurat_POD_dimplot_stress_hypoxia_cycling_features_counts_percentMT_scores_nPC",numPCs,".pdf"), width = 15, height = 15)
FeaturePlot(object = seurat_POD,
            features = c("Stress.Score1", "Hypoxia.Score1", "S.Score", "G2M.Score",
                         "nCount_RNA", "nFeature_RNA", "percent.mt", "percent.MALAT1"),
            cols = c("grey", "blue"), pt.size = 1,
            reduction = "umap")
dev.off()

###### 1.2: Harmony integration ###### 

## Using Harmony
seurat_POD <- RunHarmony(seurat_POD, "orig.ident")

## Clustering
numPCs
seurat_POD <- FindNeighbors(seurat_POD, reduction = "harmony", dims = 1:numPCs)
res = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1, 1.2, 1.5, 1.7, 2)
seurat_POD <- FindClusters(seurat_POD, resolution = res)
seurat_POD <- RunUMAP(seurat_POD, reduction = "harmony", dims = 1:numPCs)

pdf(paste0(projectFolder, "seurat_POD_harmony_dimPlot_umap_multiRes_nPC",numPCs,".pdf"), width = 10, height = 8)
for (i in 1:length(res)){
  resolution <- paste0("RNA_snn_res." , res[i])
  print(DimPlot(seurat_POD, reduction = "umap", group.by = resolution, label = T, pt.size = 0.7, raster.dpi = c(1024,1024)) + ggtitle(paste0("Resolution: ", res[i])))
}
dev.off()

## Dimplot - (to check after Harmony)
pdf(paste0(projectFolder, "seurat_POD_harmony_dimplot_sample_diagnosis_processing_nPC",numPCs,".pdf"), width = 10, height = 8)
DimPlot(seurat_POD, group.by = "orig.ident", label = T, raster=FALSE)
DimPlot(seurat_POD, group.by = "Diagnosis", label = T, raster=FALSE)
DimPlot(seurat_POD, group.by = "Date_processing", label = T, raster=FALSE)
dev.off()

## Feature plot of features, counts and %mtRNA (after Harmony)
pdf(paste0(projectFolder, "seurat_POD_harmony_featureplot_stress_hypoxia_cycling_features_counts_percentMT_nPC",numPCs,".pdf"), width = 15, height = 15)
FeaturePlot(object = seurat_POD,
            features = c("Stress.Score1", "Hypoxia.Score1", "S.Score", "G2M.Score",
                         "nCount_RNA", "nFeature_RNA", "percent.mt", "percent.MALAT1"),
            cols = c("grey", "blue"), pt.size = 1,
            reduction = "umap")
dev.off()


###### 1.3: Checking cell annotation which was performed at main annotation level ("main_annot") ###### 

# See part 1 of main annotation script for info on "main_annot"-annotation

seurat_POD <- SetIdent(seurat_POD, value = "main_annot")

pdf(paste0(projectFolder, "seurat_POD_DimPlot_main_annot.pdf"), width = 11, height = 10)
DimPlot(seurat_POD, group.by = "main_annot", label = T, pt.size = 0.7, raster.dpi = c(1024,1024))
dev.off()

## Plotting the "main annot" clusters
Lcluster <- list(
  cluster_keepfornow="keepfornow",
  cluster_CNT="CNT",
  cluster_PEC="PEC",
  cluster_IC_B="IC_B",
  cluster_T_NKT="T_NKT",
  cluster_POD="POD",
  cluster_FIB="FIB",
  cluster_IC_A="IC_A",
  cluster_TAL="TAL")

# Highlight the clusters on the dimplot 
pdf(paste0(projectFolder, "seurat_POD_harmony_main_annot_clustersHighlighted_Dimplot.pdf"), width = 10, height = 10)
lapply(seq_along(Lcluster),function(i) 
  DimPlot(seurat_POD, cells.highlight = WhichCells(seurat_POD, idents = Lcluster[[i]])) & ggtitle(names(Lcluster)[[i]]))
dev.off()

# Note that most cells - as expected - are PODs according to this annotation


###### 1.4: Checking kidney cell atlas annotation (SingleR) ###### 

# See part 2 of main annotation script for info on the used publicly available dataset and SingleR pipeline

# SOURCE: https://www.kpmp.org/doi-collection/10-48698-yyvc-ak78
# FILENAME: WashU-UCSD_HuBMAP_KPMP-Biopsy_10X-R_12032021.h5Seurat
# NAME: ATLAS EXPLORER V1.3 SINGLE-NUCLEUS RNA-SEQ DATA
# DESCRIPTION: Aggregated, clustered single-nucleus RNA-seq data used in the KPMP Atlas Explorer v1.3 - Dataset published 2021 via Kidney Precision Medicine Project
# CREATOR(S): Kidney Precision Medicine Project
# DATE: December 8, 2021
# LICENSE: Creative Commons Attribution 4.0 International (CC BY 4.0)
# DOI: https://doi.org/10.48698/yyvc-ak78
# HOW TO CITE THIS DATA:Kidney Precision Medicine Project. (2021). Aggregated, clustered single-nucleus RNA-seq data used in the KPMP Atlas Explorer v1.3. Kidney Precision Medicine Project. https://doi.org/10.48698/yyvc-ak78

# Checking kidney atlas projection l1
levels(seurat_POD$kidneyAtlas_prediction_l1)

pdf(paste0(projectFolder, "seurat_POD_harmony_dimplot_KidneyAtlas_prediction_l1.pdf"), width = 10, height = 10)
DimPlot(seurat_POD, group.by = "kidneyAtlas_prediction_l1", label = T, raster=FALSE)
dev.off()

# Checking kidney atlas projection l2
levels(seurat_POD$kidneyAtlas_prediction_l2)

pdf(paste0(projectFolder, "seurat_POD_harmony_dimplot_KidneyAtlas_prediction_l2.pdf"), width = 14, height = 10)
DimPlot(seurat_POD, group.by = "kidneyAtlas_prediction_l2", label = T, raster=FALSE)
dev.off()


###### 1.5: Choosing resolution + QC per cluster + checking doublets ###### 
ChosenRes <- "res04" #proceed with res04
levels(seurat_POD$RNA_snn_res.0.4)
seurat_POD <- SetIdent(seurat_POD, value = "RNA_snn_res.0.4")
DimPlot(seurat_POD, label=T)

# QC per cluster - VlnPlots, tables and number of cells per cluster - for code: see main annotation script
# Checking for doublets - for code: see main annotation script

###### 1.6: FindAllMarkers + plotting of DEGs ###### 

## Finding all markers ##
ChosenRes <- "res04" #proceed with res04
seurat_POD <- SetIdent(seurat_POD, value = "RNA_snn_res.0.4")

allMarkers <- FindAllMarkers(seurat_POD, only.pos = F, min.pct = 0.1, logfc.threshold = 0.25,
                             max.cells.per.ident = 10000)
allMarkers <- group_by(allMarkers, cluster)
write.table(allMarkers, file = paste0(projectFolder, "seurat_POD_harmony_", ChosenRes, "_FindAllMarkers_minPCT01_logFC025.csv"),
            row.names = F, col.names = T, sep = "\t", quote = F)

###### 1.7: Plotting markers genes w/ featureplots and dotplots (KidneyMarkersV2) ###### 

# for code: see main annotation script

###### 1.8: Plotting markers genes w/ featureplots and dotplots (KidneyMarkersKPMP) ###### 

# for code: see main annotation script

###### 1.9: Re-annotate res0.4 to proceed with subcluster clean-up ###### 
seurat_POD <- SetIdent(seurat_POD, value = "RNA_snn_res.0.4")

seurat_POD <- RenameIdents(seurat_POD,'0'='POD') # keep.
seurat_POD <- RenameIdents(seurat_POD,'1'='POD') # Possibly some DTL and TAL contamination, but not as much as other clusters. Also less healthy controls in this group. Also, not so much lowQ in this cluster (from prior annot). So keep cluster for now.
seurat_POD <- RenameIdents(seurat_POD,'2'='POD') # keep.
seurat_POD <- RenameIdents(seurat_POD,'3'='aPOD') # keep, these PODs seem more activated (although not higher nFeat).
seurat_POD <- RenameIdents(seurat_POD,'4'='rm') # high doublet score, lowQ, high PT signature so likely PT-doublet. Remove.
seurat_POD <- RenameIdents(seurat_POD,'5'='aPOD') # keep, these PODs seem more activated (although not higher nFeat).
seurat_POD <- RenameIdents(seurat_POD,'6'='rm') # high doublet score, lowQ, CNT-contamination? Remove.
seurat_POD <- RenameIdents(seurat_POD,'7'='rm') # lower nFeat, very high mtRNA. Remove.
seurat_POD <- RenameIdents(seurat_POD,'8'='POD_IMM') # keep for now, but possibly partial with immune. To check at second round of QC.
seurat_POD <- RenameIdents(seurat_POD,'9'='POD_PEC') # Possible PEC contamination, possible POD/PEC overlap, keep for now
seurat_POD <- RenameIdents(seurat_POD,'10'='rm') #  high doublet score, lowQ, TAL-doublets. Remove.
seurat_POD <- RenameIdents(seurat_POD,'11'='rm') # lowQ, EC-doublet. Remove.
seurat_POD <- RenameIdents(seurat_POD,'12'='POD') # lower nFeat.
seurat_POD <- RenameIdents(seurat_POD,'13'='rm') # lowQ, IC doublet. Remove.

DimPlot(seurat_POD, label=T)

## Cluster IDs in 'POD_annot'-column
seurat_POD$POD_annot <- Idents(seurat_POD)
seurat_POD$POD_annot <- factor(seurat_POD$POD_annot)
seurat_POD <- SetIdent(seurat_POD, value = "POD_annot")

pdf(paste0(projectFolder, "seurat_POD_harmony_", ChosenRes, "_POD_annot_dimPlot_umap.pdf"), width = 10, height = 8)
DimPlot(seurat_POD, group.by = "POD_annot", reduction = "umap", label = TRUE, pt.size = 0.7, raster.dpi = c(1024,1024)) + ggtitle(paste0("seurat_POD_harmony_", ChosenRes, "_POD_annot"))
dev.off()


###### 1.10: Selecting HQ-PODs ###### 
seurat_POD <- SetIdent(seurat_POD, value = "POD_annot")
DimPlot(seurat_POD, label = T)
seurat_POD_HQ <- subset(seurat_POD, ident = c("POD", "aPOD", "POD_IMM","POD_PEC"))

pdf(paste0(projectFolder, "seurat_POD_harmony_", ChosenRes, "_POD_annot_afterRemovalLowQ_DimPlot.pdf"), width = 10, height = 8)
DimPlot(seurat_POD_HQ, group.by = "POD_annot", reduction = "umap", label = TRUE, pt.size = 0.7, raster.dpi = c(1024,1024))
dev.off()

rm(seurat_POD)

# This concludes part 1.


##### 2. PART 2 - Re-Running pipeline on seurat_POD_HQ - first lowQ already exluded #####

## Setting projectfolder:
projectFolder <- "path/seurat_POD_HQ/"

seurat_POD_HQ
# An object of class Seurat 
# 36780 features across 2715 samples within 1 assay 
# Active assay: RNA (36780 features, 2000 variable features)
# 3 dimensional reductions calculated: pca, umap, harmony

###### 2.1: Re-Running Seurat pipeline ###### 

## Normalizing the data
seurat_POD_HQ <- NormalizeData(seurat_POD_HQ, normalization.method = "LogNormalize", scale.factor = 10000)

## Highly variable features
seurat_POD_HQ <- FindVariableFeatures(seurat_POD_HQ, selection.method = "vst", nfeatures = 2000)

## Scaling the data + regressing for mt.RNA, nCount, nFeature
seurat_POD_HQ <- ScaleData(seurat_POD_HQ, vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt")) # we only use variable features

## PCA - linear dimensional reduction
seurat_POD_HQ <- RunPCA(seurat_POD_HQ, features = VariableFeatures(object = seurat_POD_HQ), npcs = 50)

pdf(paste0(projectFolder, "seurat_POD_HQ_noHarmony_elbowPlot.pdf"), width = 10, height = 10)
ElbowPlot(seurat_POD_HQ, n = 50)
dev.off()

## Clustering
numPCs <- 11
seurat_POD_HQ <- FindNeighbors(seurat_POD_HQ, dims = 1:numPCs)
res = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1, 1.2, 1.5, 1.7, 2)
seurat_POD_HQ <- FindClusters(seurat_POD_HQ, resolution = res)
seurat_POD_HQ <- RunUMAP(seurat_POD_HQ, dims = 1:numPCs)

pdf(paste0(projectFolder, "seurat_POD_HQ_noHarmony_dimPlot_umap_multiRes_nPC",numPCs,".pdf"), width = 10, height = 8)
for (i in 1:length(res)){
  resolution <- paste0("RNA_snn_res." , res[i])
  print(DimPlot(seurat_POD_HQ, reduction = "umap", group.by = resolution, label = T, pt.size = 0.7, raster.dpi = c(1024,1024)) + ggtitle(paste0("Resolution: ", res[i])))
}
dev.off()

###### 2.2: Harmony integration ###### 

## Using Harmony
seurat_POD_HQ <- RunHarmony(seurat_POD_HQ, "orig.ident")

## Clustering
numPCs
seurat_POD_HQ <- FindNeighbors(seurat_POD_HQ, reduction = "harmony", dims = 1:numPCs)
res = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1, 1.2, 1.5, 1.7, 2)
seurat_POD_HQ <- FindClusters(seurat_POD_HQ, resolution = res)
seurat_POD_HQ <- RunUMAP(seurat_POD_HQ, reduction = "harmony", dims = 1:numPCs)

pdf(paste0(projectFolder, "seurat_POD_HQ_harmony_dimPlot_umap_multiRes_nPC",numPCs,".pdf"), width = 10, height = 8)
for (i in 1:length(res)){
  resolution <- paste0("RNA_snn_res." , res[i])
  print(DimPlot(seurat_POD_HQ, reduction = "umap", group.by = resolution, label = T, pt.size = 0.7, raster.dpi = c(1024,1024)) + ggtitle(paste0("Resolution: ", res[i])))
}
dev.off()

###### 2.3: Choosing resolution + QC per cluster ###### 
ChosenRes <- "CustomRes_02_20" #proceed with combination of res0.2 and res2.0
seurat_POD_HQ <- SetIdent(seurat_POD_HQ, value = "RNA_snn_res.2")
DimPlot(seurat_POD_HQ, label=T)

seurat_POD_HQ <- RenameIdents(seurat_POD_HQ,'0'='POD') # 
seurat_POD_HQ <- RenameIdents(seurat_POD_HQ,'1'='1') # 
seurat_POD_HQ <- RenameIdents(seurat_POD_HQ,'2'='2') # 
seurat_POD_HQ <- RenameIdents(seurat_POD_HQ,'3'='POD') # 
seurat_POD_HQ <- RenameIdents(seurat_POD_HQ,'4'='POD') # 
seurat_POD_HQ <- RenameIdents(seurat_POD_HQ,'5'='5') # 
seurat_POD_HQ <- RenameIdents(seurat_POD_HQ,'6'='6') # 
seurat_POD_HQ <- RenameIdents(seurat_POD_HQ,'7'='POD') # 
seurat_POD_HQ <- RenameIdents(seurat_POD_HQ,'8'='POD') # 
seurat_POD_HQ <- RenameIdents(seurat_POD_HQ,'9'='POD') # 
seurat_POD_HQ <- RenameIdents(seurat_POD_HQ,'10'='10') # 
seurat_POD_HQ <- RenameIdents(seurat_POD_HQ,'11'='11') # 
seurat_POD_HQ <- RenameIdents(seurat_POD_HQ,'12'='POD') # 
seurat_POD_HQ <- RenameIdents(seurat_POD_HQ,'13'='13') # 
seurat_POD_HQ <- RenameIdents(seurat_POD_HQ,'14'='14') # 
seurat_POD_HQ <- RenameIdents(seurat_POD_HQ,'15'='15') # 
seurat_POD_HQ <- RenameIdents(seurat_POD_HQ,'16'='16') # 
seurat_POD_HQ <- RenameIdents(seurat_POD_HQ,'17'='17') # 
seurat_POD_HQ <- RenameIdents(seurat_POD_HQ,'18'='18') # 
seurat_POD_HQ <- RenameIdents(seurat_POD_HQ,'19'='POD') #

DimPlot(seurat_POD_HQ, label=T)

## Cluster IDs in 'POD_annot'-column
seurat_POD_HQ$CustomRes_02_20 <- Idents(seurat_POD_HQ)
seurat_POD_HQ$CustomRes_02_20 <- factor(seurat_POD_HQ$CustomRes_02_20, levels = c("POD", "1", "2", "5", "6", "10", "11", "13", "14", "15", "16", "17", "18"))
levels(seurat_POD_HQ$CustomRes_02_20)
seurat_POD_HQ <- SetIdent(seurat_POD_HQ, value = "CustomRes_02_20")

pdf(paste0(projectFolder, "seurat_POD_HQ_harmony_", ChosenRes, "_dimPlot_umap.pdf"), width = 10, height = 8)
DimPlot(seurat_POD_HQ, group.by = "CustomRes_02_20", reduction = "umap", label = TRUE, pt.size = 0.7, raster.dpi = c(1024,1024)) + ggtitle(paste0("seurat_POD_HQ_harmony_", ChosenRes))
dev.off()

# QC per cluster - VlnPlots, tables and number of cells per cluster - for code: see main annotation script
# Checking for doublets - for code: see main annotation script

###### 2.4: FindAllMarkers + plotting of DEGs ###### 

## Finding all markers ##
ChosenRes <- "CustomRes_02_20" 

seurat_POD_HQ <- SetIdent(seurat_POD_HQ, value = "CustomRes_02_20")
allMarkers <- FindAllMarkers(seurat_POD_HQ, only.pos = F, min.pct = 0.1, logfc.threshold = 0.25,
                             max.cells.per.ident = 10000)
allMarkers <- group_by(allMarkers, cluster)
write.table(allMarkers, file = paste0(projectFolder, "seurat_POD_HQ_harmony_", ChosenRes, "_FindAllMarkers_minPCT01_logFC025.csv"),
            row.names = F, col.names = T, sep = "\t", quote = F)

###### 2.5: Plotting markers genes w/ featureplots and dotplots (KidneyMarkersV2) ###### 

# for code: see main annotation script

###### 2.6: Plotting markers genes w/ featureplots and dotplots (KidneyMarkersKPMP) ###### 

# for code: see main annotation script

###### 2.7: Re-annotate CustomRes_02_20 to proceed with subcluster clean-up ###### 
seurat_POD_HQ <- SetIdent(seurat_POD_HQ, value = "CustomRes_02_20")

seurat_POD_HQ <- RenameIdents(seurat_POD_HQ,'POD'='POD') # 
seurat_POD_HQ <- RenameIdents(seurat_POD_HQ,'1'='POD') # very few DEGs! But otherwise seems normal POD, so we do not remove. Currently not yet clear why they cluster separately.
seurat_POD_HQ <- RenameIdents(seurat_POD_HQ,'2'='POD') # low nFeat, but otherwise normal POD, so we do not remove (a small fraction with high ribo, but most seem normal PODs).
seurat_POD_HQ <- RenameIdents(seurat_POD_HQ,'5'='POD') # keep.
seurat_POD_HQ <- RenameIdents(seurat_POD_HQ,'6'='aPOD') # keep, seems activated POD.
seurat_POD_HQ <- RenameIdents(seurat_POD_HQ,'10'='aPOD') # keep, seems activated POD.
seurat_POD_HQ <- RenameIdents(seurat_POD_HQ,'11'='POD_PEC') # PEC markers.
seurat_POD_HQ <- RenameIdents(seurat_POD_HQ,'13'='rm') # best to remove (doublet/ambient), marker genes of immune cells.
seurat_POD_HQ <- RenameIdents(seurat_POD_HQ,'14'='rm') # best to remove (doublet/ambient), marker genes of PT.
seurat_POD_HQ <- RenameIdents(seurat_POD_HQ,'15'='rm') # very few DEGs, remove.
seurat_POD_HQ <- RenameIdents(seurat_POD_HQ,'16'='POD') # low nFeat, high MALAT1, keep for now, but no DEGs.
seurat_POD_HQ <- RenameIdents(seurat_POD_HQ,'17'='rm') # high MALAT1, remove.
seurat_POD_HQ <- RenameIdents(seurat_POD_HQ,'18'='rm') # high doublet, remove.

DimPlot(seurat_POD_HQ, label=T)

## Cluster IDs in 'POD_HQ_annot'-column
seurat_POD_HQ$POD_HQ_annot <- Idents(seurat_POD_HQ)
seurat_POD_HQ$POD_HQ_annot <- factor(seurat_POD_HQ$POD_HQ_annot)
seurat_POD_HQ <- SetIdent(seurat_POD_HQ, value = "POD_HQ_annot")

pdf(paste0(projectFolder, "seurat_POD_HQ_harmony_", ChosenRes, "_POD_HQ_annot_dimPlot_umap.pdf"), width = 10, height = 8)
DimPlot(seurat_POD_HQ, group.by = "POD_HQ_annot", reduction = "umap", label = TRUE, pt.size = 0.7, raster.dpi = c(1024,1024)) + ggtitle(paste0("seurat_POD_HQ_harmony_", ChosenRes, "_POD_HQ_annot"))
dev.off()

###### 2.8: Selecting HQ2-PODs ###### 
seurat_POD_HQ <- SetIdent(seurat_POD_HQ, value = "POD_HQ_annot")
seurat_POD_HQ2 <- subset(seurat_POD_HQ, ident = c("POD", "aPOD","POD_PEC"))

pdf(paste0(projectFolder, "seurat_POD_HQ_harmony_", ChosenRes, "_POD_HQ_annot_afterRemovalLowQ_DimPlot.pdf"), width = 10, height = 8)
DimPlot(seurat_POD_HQ2, group.by = "POD_HQ_annot", reduction = "umap", label = TRUE, pt.size = 0.7, raster.dpi = c(1024,1024))
dev.off()

rm(seurat_POD_HQ)

# This concludes part 2.


##### 3. PART 3 - Re-Running pipeline on seurat_POD_HQ2 - lowQ already exluded (2x) #####

## Setting projectfolder:
projectFolder <- "path/seurat_POD_HQ2/"

seurat_POD_HQ2
# An object of class Seurat 
# 36780 features across 2471 samples within 1 assay 
# Active assay: RNA (36780 features, 2000 variable features)
# 3 dimensional reductions calculated: pca, umap, harmony

###### 3.1: Re-Running Seurat pipeline ###### 

## Normalizing the data
seurat_POD_HQ2 <- NormalizeData(seurat_POD_HQ2, normalization.method = "LogNormalize", scale.factor = 10000)

## Highly variable features
seurat_POD_HQ2 <- FindVariableFeatures(seurat_POD_HQ2, selection.method = "vst", nfeatures = 2000)

## Scaling the data + regressing for mt.RNA, nCount, nFeature
seurat_POD_HQ2 <- ScaleData(seurat_POD_HQ2, vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt")) # we only use variable features

## PCA - linear dimensional reduction
seurat_POD_HQ2 <- RunPCA(seurat_POD_HQ2, features = VariableFeatures(object = seurat_POD_HQ2), npcs = 50)

pdf(paste0(projectFolder, "seurat_POD_HQ2_noHarmony_elbowPlot.pdf"), width = 10, height = 10)
ElbowPlot(seurat_POD_HQ2, n = 50)
dev.off()

## Clustering
numPCs <- 13
seurat_POD_HQ2 <- FindNeighbors(seurat_POD_HQ2, dims = 1:numPCs)
res = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1, 1.2, 1.5, 1.7, 2)
seurat_POD_HQ2 <- FindClusters(seurat_POD_HQ2, resolution = res)
seurat_POD_HQ2 <- RunUMAP(seurat_POD_HQ2, dims = 1:numPCs)

pdf(paste0(projectFolder, "seurat_POD_HQ2_noHarmony_dimPlot_umap_multiRes_nPC",numPCs,".pdf"), width = 10, height = 8)
for (i in 1:length(res)){
  resolution <- paste0("RNA_snn_res." , res[i])
  print(DimPlot(seurat_POD_HQ2, reduction = "umap", group.by = resolution, label = T, pt.size = 0.7, raster.dpi = c(1024,1024)) + ggtitle(paste0("Resolution: ", res[i])))
}
dev.off()

## First change some meta.data columns (numeric vs. categorical data)
class(seurat_POD_HQ2$Age_bx)
seurat_POD_HQ2$Age_bx <- as.numeric(seurat_POD_HQ2$Age_bx)
class(seurat_POD_HQ2$Age_bx)

class(seurat_POD_HQ2$Prot_UPCR_bx)
seurat_POD_HQ2$Prot_UPCR_bx <- as.numeric(seurat_POD_HQ2$Prot_UPCR_bx)
class(seurat_POD_HQ2$Prot_UPCR_bx)

class(seurat_POD_HQ2$sAlb_bx)
seurat_POD_HQ2$sAlb_bx <- as.numeric(seurat_POD_HQ2$sAlb_bx)
class(seurat_POD_HQ2$sAlb_bx)

class(seurat_POD_HQ2$eGFR_bx)
seurat_POD_HQ2$eGFR_bx <- as.numeric(seurat_POD_HQ2$eGFR_bx)
class(seurat_POD_HQ2$eGFR_bx)

class(seurat_POD_HQ2$NS_bx)
seurat_POD_HQ2$NS_bx <- factor(seurat_POD_HQ2$NS_bx, levels = c("Yes","No"))
levels(seurat_POD_HQ2$NS_bx)

## Dimplot of clinical metadata
pdf(paste0(projectFolder, "seurat_POD_HQ2_noHarmony_clinical_metadata_DimFeaturePlot.pdf"), width = 7, height = 6)
DimPlot(seurat_POD_HQ2, group.by = "orig.ident", label = T, raster=FALSE)
DimPlot(seurat_POD_HQ2, group.by = "Diagnosis", label = T, raster=FALSE)
DimPlot(seurat_POD_HQ2, group.by = "Date_processing", label = T, raster=FALSE)
DimPlot(seurat_POD_HQ2, group.by = "Sex", label = T, raster=FALSE)
FeaturePlot(seurat_POD_HQ2, features = "Age_bx")
DimPlot(seurat_POD_HQ2, group.by = "Weight_cat_bx", label = T, raster=FALSE)
DimPlot(seurat_POD_HQ2, group.by = "ACE_ARB_bx", label = T, raster=FALSE)
FeaturePlot(seurat_POD_HQ2, features = "Prot_UPCR_bx")
DimPlot(seurat_POD_HQ2, group.by = "Prot_UPCR_bx_cat", label = T, raster=FALSE)
DimPlot(seurat_POD_HQ2, group.by = "Hematuria_bx", label = T, raster=FALSE)
FeaturePlot(seurat_POD_HQ2, features = "sAlb_bx")
DimPlot(seurat_POD_HQ2, group.by = "sAlb_bx_cat", label = T, raster=FALSE)
DimPlot(seurat_POD_HQ2, group.by = "NS_bx", label = T, raster=FALSE)
FeaturePlot(seurat_POD_HQ2, features = "eGFR_bx")
DimPlot(seurat_POD_HQ2, group.by = "eGFR_bx_cat", label = T, raster=FALSE)
DimPlot(seurat_POD_HQ2, group.by = "GS_score", label = T, raster=FALSE)
DimPlot(seurat_POD_HQ2, group.by = "IF_score", label = T, raster=FALSE)
DimPlot(seurat_POD_HQ2, group.by = "TA_score", label = T, raster=FALSE)
DimPlot(seurat_POD_HQ2, group.by = "CV_score", label = T, raster=FALSE)
DimPlot(seurat_POD_HQ2, group.by = "MCCS", label = T, raster=FALSE)
DimPlot(seurat_POD_HQ2, group.by = "Relapse", label = T, raster=FALSE)
DimPlot(seurat_POD_HQ2, group.by = "Death", label = T, raster=FALSE)
DimPlot(seurat_POD_HQ2, group.by = "RNA_snn_res.0.2", label = T, raster=FALSE)
dev.off()

###### 3.2: Harmony integration ###### 

## Using Harmony
seurat_POD_HQ2 <- RunHarmony(seurat_POD_HQ2, "orig.ident")

## Clustering
numPCs
seurat_POD_HQ2 <- FindNeighbors(seurat_POD_HQ2, reduction = "harmony", dims = 1:numPCs)
res = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1, 1.2, 1.5, 1.7, 2)
seurat_POD_HQ2 <- FindClusters(seurat_POD_HQ2, resolution = res)
seurat_POD_HQ2 <- RunUMAP(seurat_POD_HQ2, reduction = "harmony", dims = 1:numPCs)

pdf(paste0(projectFolder, "seurat_POD_HQ2_harmony_dimPlot_umap_multiRes_nPC",numPCs,".pdf"), width = 10, height = 8)
for (i in 1:length(res)){
  resolution <- paste0("RNA_snn_res." , res[i])
  print(DimPlot(seurat_POD_HQ2, reduction = "umap", group.by = resolution, label = T, pt.size = 0.7, raster.dpi = c(1024,1024)) + ggtitle(paste0("Resolution: ", res[i])))
}
dev.off()

###### 3.3: Choosing resolution + QC per cluster ###### 
ChosenRes <- "res02" #proceed with res0.2
seurat_POD_HQ2$RNA_snn_res.0.2 <- factor(seurat_POD_HQ2$RNA_snn_res.0.2)
levels(seurat_POD_HQ2$RNA_snn_res.0.2)
seurat_POD_HQ2 <- SetIdent(seurat_POD_HQ2, value = "RNA_snn_res.0.2")
DimPlot(seurat_POD_HQ2, label=T)

pdf(paste0(projectFolder, "seurat_POD_HQ2_harmony_", ChosenRes, "_dimPlot_umap.pdf"), width = 10, height = 8)
DimPlot(seurat_POD_HQ2, group.by = "RNA_snn_res.0.2", reduction = "umap", label = TRUE, pt.size = 0.7, raster.dpi = c(1024,1024)) + ggtitle(paste0("seurat_POD_HQ2_harmony_", ChosenRes))
dev.off()


###### 3.4: Re-annotate RNA_snn_res.0.2 ###### 
seurat_POD_HQ2 <- SetIdent(seurat_POD_HQ2, value = "RNA_snn_res.0.2")

seurat_POD_HQ2 <- RenameIdents(seurat_POD_HQ2,'0'='POD1') # 
seurat_POD_HQ2 <- RenameIdents(seurat_POD_HQ2,'1'='POD2') # higher nFeat
seurat_POD_HQ2 <- RenameIdents(seurat_POD_HQ2,'2'='POD1') # aPOD (high ribo, HLA-upreg, ...), note that cluster does not have higher nFeat
seurat_POD_HQ2 <- RenameIdents(seurat_POD_HQ2,'3'='POD1') # higher nFeat, aPOD
seurat_POD_HQ2 <- RenameIdents(seurat_POD_HQ2,'4'='POD3') # PEC phenotype, difficult to know whether POD-PEC overlap or doublet

DimPlot(seurat_POD_HQ2, label=T)

## Cluster IDs in 'POD_HQ2_annot'-column
seurat_POD_HQ2$POD_HQ2_annot <- Idents(seurat_POD_HQ2)
seurat_POD_HQ2$POD_HQ2_annot <- factor(seurat_POD_HQ2$POD_HQ2_annot, levels = c("POD1", "POD2", "POD3"))
seurat_POD_HQ2 <- SetIdent(seurat_POD_HQ2, value = "POD_HQ2_annot")

pdf(paste0(projectFolder, "seurat_POD_HQ2_harmony_", ChosenRes, "_POD_HQ2_annot_dimPlot_umap.pdf"), width = 10, height = 8)
DimPlot(seurat_POD_HQ2, group.by = "POD_HQ2_annot", reduction = "umap", label = TRUE, pt.size = 0.7, raster.dpi = c(1024,1024)) + ggtitle(paste0("seurat_POD_HQ2_harmony_", ChosenRes, "_POD_HQ2_annot"))
dev.off()

###### 3.5: QC on SAMPLES (output included in manuscript) ###### 

# QC Per sample - VlnPlot
seurat_POD_HQ2 <- SetIdent(seurat_POD_HQ2, value = "orig.ident")

pdf(paste0(projectFolder, "seurat_POD_HQ2_QC_perSample_VlnPlot.pdf"), width = 14, height = 5)
VlnPlot(seurat_POD_HQ2, features = c("nCount_RNA"),
        pt.size = 0, 
        cols = rep("#1F78B4", times = 25)) &
  geom_boxplot(width=0.1, fill="white", outlier.size = 0.5) & 
  NoLegend() &
  stat_summary(fun.y = median, geom='point', size = 2, colour = "black")

VlnPlot(seurat_POD_HQ2, features = c("nFeature_RNA"), 
        pt.size = 0, 
        cols = rep("#1F78B4", times = 25)) &
  geom_boxplot(width=0.1, fill="white", outlier.size = 0.5) & 
  NoLegend() &
  stat_summary(fun.y = median, geom='point', size = 2, colour = "black")

VlnPlot(seurat_POD_HQ2, features = c("percent.mt"), 
        pt.size = 0, 
        cols = rep("#1F78B4", times = 25)) &
  geom_boxplot(width=0.1, fill="white", outlier.size = 0.5) & 
  NoLegend() &
  stat_summary(fun.y = median, geom='point', size = 2, colour = "black")

VlnPlot(seurat_POD_HQ2, features = c("percent.MALAT1"), 
        pt.size = 0, 
        cols = rep("#1F78B4", times = 25)) &
  geom_boxplot(width=0.1, fill="white", outlier.size = 0.5) & 
  NoLegend() &
  stat_summary(fun.y = median, geom='point', size = 2, colour = "black")

dev.off()

## QC per sample, table dplyr style
df <- data.frame(matrix(ncol = 5, nrow = length(seurat_POD_HQ2$orig.ident)))
x <- c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "Diagnosis")
colnames(df) <- x

df$orig.ident <- seurat_POD_HQ2$orig.ident
df$nCount_RNA <- seurat_POD_HQ2$nCount_RNA
df$nFeature_RNA <- seurat_POD_HQ2$nFeature_RNA
df$percent.mt <- seurat_POD_HQ2$percent.mt
df$percent.MALAT1 <- seurat_POD_HQ2$percent.MALAT1
df$Diagnosis <- seurat_POD_HQ2$Diagnosis

seuratDF <- df %>%
  group_by(orig.ident) %>%
  summarize(  Diagnosis = unique(Diagnosis),
              Median_nCounts = round(median(nCount_RNA),0), 
              Median_nFeatures = round(median(nFeature_RNA),0),
              Median_percent.mt = round(median(percent.mt),2),
              nCells = n()
  )

colnames(seuratDF) <- c("orig.ident","Diagnosis","Median\nnCounts", "Median\nnFeatures", "Median\npercent.mt", "Nuclei")
seuratDF

pdf(paste0(projectFolder, "seurat_POD_HQ2_QC_perSample_table_numberOfCellsCountsFeaturesMito_withoutMALAT1.pdf"), height=8, width=13)
p <- tableGrob(seuratDF)
grid.arrange(p)
dev.off()

colnames(seuratDF) <- c("orig.ident","Diagnosis","Median nCounts", "Median nFeatures", "Median percent.mt", "Nuclei")
seuratDF

write.table(seuratDF, file = paste0(projectFolder, "seurat_POD_HQ2_QC_perSample_table_numberOfCellsCountsFeaturesMito_withoutMALAT1.csv"),
            quote = F, sep = ";", col.names = T, row.names = T)

###### 3.6: QC on DIAGNOSES (output included in manuscript) ###### 

# QC per Diagnosis - VlnPlot
seurat_POD_HQ2 <- SetIdent(seurat_POD_HQ2, value = "Diagnosis")

pdf(paste0(projectFolder, "seurat_POD_HQ2_QC_perDiagnosis_VlnPlot.pdf"), width = 14, height = 5)
VlnPlot(seurat_POD_HQ2, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.MALAT1"),
        pt.size = 0, ncol=5,
        cols = rep("#1F78B4", times = 4)) &
  geom_boxplot(width=0.1, fill="white", outlier.size = 0.5) & 
  NoLegend() &
  stat_summary(fun.y = median, geom='point', size = 2, colour = "black")

dev.off()

## QC per sample - table
seuratDF <- data.frame(Median_nCounts = round(tapply(seurat_POD_HQ2$nCount_RNA, seurat_POD_HQ2$Diagnosis, median)),
                       Median_nFeatures = round(tapply(seurat_POD_HQ2$nFeature_RNA, seurat_POD_HQ2$Diagnosis, median)),
                       Median_percent.mt = round(tapply(seurat_POD_HQ2$percent.mt, seurat_POD_HQ2$Diagnosis, median), 2),
                       nCells = table(seurat_POD_HQ2$Diagnosis))
rownames(seuratDF) <- seuratDF$nCells.Var1
seuratDF$nCells.Var1 <- NULL
colnames(seuratDF) <- c("Median\nnCounts", "Median\nnFeatures", "Median\npercent.mt", "Nuclei")

pdf(paste0(projectFolder, "seurat_POD_HQ2_QC_perDiagnosis_table_numberOfCellsCountsFeaturesMito_withoutMALAT1.pdf"), height=4, width=12)
p <- tableGrob(seuratDF)
grid.arrange(p)
dev.off()

colnames(seuratDF) <- c("Median nCounts", "Median nFeatures", "Median percent.mt", "Nuclei")
seuratDF

write.table(seuratDF, file = paste0(projectFolder, "seurat_POD_HQ2_QC_perDiagnosis_table_numberOfCellsCountsFeaturesMito_withoutMALAT1.csv"),
            quote = F, sep = ";", col.names = T, row.names = T)


###### 3.7: QC on CLUSTERS (POD_HQ2_annot, i.e. annotated RNA_snn_res.0.2) ###### 

# QC per cluster - VlnPlot
seurat_POD_HQ2 <- SetIdent(seurat_POD_HQ2, value = "POD_HQ2_annot")
ChosenRes <- "POD_HQ2_annot"

pdf(paste0(projectFolder, "seurat_POD_HQ2_QC_perCluster_VlnPlot_", ChosenRes, ".pdf"), width = 4, height = 7)
VlnPlot(seurat_POD_HQ2, features = c("nCount_RNA"),
        pt.size = 0, 
        cols = rep("#1F78B4", times = 3)) &
  geom_boxplot(width=0.1, fill="white", outlier.size = 0.5) & NoLegend() &
  stat_summary(fun.y = median, geom='point', size = 2, colour = "black")

VlnPlot(seurat_POD_HQ2, features = c("nFeature_RNA"), 
        pt.size = 0, 
        cols = rep("#1F78B4", times = 3)) &
  geom_boxplot(width=0.1, fill="white", outlier.size = 0.5) & NoLegend() &
  stat_summary(fun.y = median, geom='point', size = 2, colour = "black")

VlnPlot(seurat_POD_HQ2, features = c("percent.mt"), 
        pt.size = 0, 
        cols = rep("#1F78B4", times = 3)) &
  geom_boxplot(width=0.1, fill="white", outlier.size = 0.5) & NoLegend() &
  stat_summary(fun.y = median, geom='point', size = 2, colour = "black")

VlnPlot(seurat_POD_HQ2, features = c("percent.MALAT1"), 
        pt.size = 0, 
        cols = rep("#1F78B4", times = 3)) &
  geom_boxplot(width=0.1, fill="white", outlier.size = 0.5) & 
  NoLegend() &
  stat_summary(fun.y = median, geom='point', size = 2, colour = "black")

dev.off()


# QC per cluster - table
seuratDF <- data.frame(Median_nCounts = round(tapply(seurat_POD_HQ2$nCount_RNA, seurat_POD_HQ2$POD_HQ2_annot, median)),
                       Mean_nCounts = round(tapply(seurat_POD_HQ2$nCount_RNA, seurat_POD_HQ2$POD_HQ2_annot, mean)),
                       Median_nFeatures = round(tapply(seurat_POD_HQ2$nFeature_RNA, seurat_POD_HQ2$POD_HQ2_annot, median)),
                       Mean_nFeatures = round(tapply(seurat_POD_HQ2$nFeature_RNA, seurat_POD_HQ2$POD_HQ2_annot, mean)),
                       Median_percent.mt = round(tapply(seurat_POD_HQ2$percent.mt, seurat_POD_HQ2$POD_HQ2_annot, median), 2),
                       Mean_percent.mt = round(tapply(seurat_POD_HQ2$percent.mt, seurat_POD_HQ2$POD_HQ2_annot, mean), 2),
                       Median_percent.MALAT1 = round(tapply(seurat_POD_HQ2$percent.MALAT1, seurat_POD_HQ2$POD_HQ2_annot, median), 2),
                       Mean_percent.MALAT1 = round(tapply(seurat_POD_HQ2$percent.MALAT1, seurat_POD_HQ2$POD_HQ2_annot, mean), 2),
                       nCells = table(seurat_POD_HQ2$POD_HQ2_annot))
rownames(seuratDF) <- seuratDF$nCells.Var1
seuratDF$nCells.Var1 <- NULL
colnames(seuratDF) <- c("Median\nnCounts", "Mean\nnCounts", "Median\nnFeatures", "Mean\nnFeatures", "Median\npercent.mt", "Mean\npercent.mt", "Median\npercent.MALAT1", "Mean\npercent.MALAT1", "Nuclei")

pdf(paste0(projectFolder, "seurat_POD_HQ2_QC_perCluster_table_", ChosenRes, "_numberOfCellsCountsFeaturesMito.pdf"), height=4, width=10)
p <- tableGrob(seuratDF)
grid.arrange(p)
dev.off()

colnames(seuratDF) <- c("Median nCounts", "Mean nCounts", "Median nFeatures", "Mean nFeatures", "Median percent.mt", "Mean percent.mt", "Median percent.MALAT1", "Mean percent.MALAT1", "Nuclei")
seuratDF

write.table(seuratDF, file = paste0(projectFolder, "seurat_POD_HQ2_QC_perCluster_table_", ChosenRes, "_numberOfCellsCountsFeaturesMito.csv"),
            quote = F, sep = ";", col.names = T, row.names = T)


## Number of nuclei per cluster per sample
table1 <- table(seurat_POD_HQ2$POD_HQ2_annot, seurat_POD_HQ2$orig.ident)
total <- colSums(table1)
table1_with_total <- rbind(table1, Total = total)
write.table(table1_with_total, 
            file = paste0(projectFolder, "seurat_POD_HQ2_cells_perCluster_table.csv"),
            quote = FALSE, sep = ";", col.names = TRUE, row.names = TRUE)


###### 3.8: FindAllMarkers + plotting of DEGs on POD_HQ2_annot (i.e. annotated POD_HQ2_annot) ###### 

## Finding all markers ##
ChosenRes <- "POD_HQ2_annot" 

seurat_POD_HQ2 <- SetIdent(seurat_POD_HQ2, value = "POD_HQ2_annot")
allMarkers <- FindAllMarkers(seurat_POD_HQ2, only.pos = F, min.pct = 0.1, logfc.threshold = 0.25,
                             max.cells.per.ident = 10000)
allMarkers <- group_by(allMarkers, cluster)
write.table(allMarkers, file = paste0(projectFolder, "seurat_POD_HQ2_harmony_", ChosenRes, "_FindAllMarkers_minPCT01_logFC025.csv"),
            row.names = F, col.names = T, sep = "\t", quote = F)

## Heatmap of TOP markers per cluster ##
## Reading in marker genes - TOP10
allSignMarkers <- allMarkers[allMarkers$p_val_adj < 0.05,]
top_10 <- top_n(allSignMarkers, n = 10, wt = avg_log2FC)
top_10 <- top_10[order(top_10$cluster),]
top10Expr <- as.data.frame(AverageExpression(seurat_POD_HQ2, assays = "RNA", features = unique(top_10$gene), slot = "data")) #calculates average expression in clusters

top10Zscore <- t(scale(t(top10Expr)))
top10Zscore <- top10Zscore[!is.na(top10Zscore[,1]),]

# Heatmap of z-scores for top 10 marker genes
breaks <- c(seq(min(top10Zscore), 0, -min(top10Zscore)/20),
            seq(max(top10Zscore)/21, max(top10Zscore), max(top10Zscore)/20))
myColor <- rev(colorRampPalette(brewer.pal(n = 9, name = "RdBu"))(40))

pdf(paste0(projectFolder, "seurat_POD_HQ2_harmony_", ChosenRes, "_top10_heatmap_zscore.pdf"), width = 4, height = 10)
pheatmap(top10Zscore, cluster_rows = F, cluster_cols = F, color = myColor, fontsize_row = 5,
         breaks = breaks) #, scale = "column"
dev.off()

# dotplot of top10 genes per cluster
pdf(paste0(projectFolder, "seurat_POD_HQ2_harmony_", ChosenRes, "_top10_dotplot.pdf"), width = 8, height = 15)
DotPlot(seurat_POD_HQ2, features = rev(unique(top_10$gene))) + coord_flip()
dev.off()

# Save the annotated final object
saveRDS(seurat_POD_HQ2, file = paste0(projectFolder, "seurat_POD_HQ2_harmony_annot.Rds"))



##### 4. PART 4 - DGE ANALYSIS BETWEEN DIAGNOSES - EdgeR (pseudobulk)  #####

# Switch to R-version with installation of EdgeR

# Step 1: run R 4.3
library(Seurat)
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(edgeR)
library(EnhancedVolcano)

###### 4.1: Run general EdgeR DEG analysis ###### 

## Setting projectfolder:
projectFolder <- "path/seurat_POD_HQ2/Diagnosis_analysis/EdgeR/"

# Read in the dataset
seurat_POD_HQ2 <- readRDS(file = "path/seurat_POD_HQ2/seurat_POD_HQ2_harmony_annot.Rds")
seurat_obj <- seurat_POD_HQ2 #insert object
ChosenCellGroup <- "seurat_POD_HQ2"

# Check conditions
table(seurat_obj$Diagnosis)
table(seurat_obj$orig.ident,seurat_obj$Diagnosis)

# Create pseudo-bulk samples
y <- Seurat2PB(seurat_obj, sample="orig.ident", cluster="Diagnosis")

dim(y) #36780    25
head(y$samples, n=10L)
View(y$samples)

write.table(y$samples, file = paste0(projectFolder, ChosenCellGroup, "_EdgeR_pseudobulk_librarySizes.csv"),
            row.names = T, col.names = T, sep = "\t", quote = F)

# We first examine the library sizes of all the pseudo-bulk samples and filter out those below the threshold of 50,000.
summary(y$samples$lib.size)

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 19127  123440  303602  612920  937913 4107317

# Using aggregate to summarize lib.size by cluster
summary_by_cluster <- aggregate(lib.size ~ cluster, data = y$samples, FUN = summary)
print(summary_by_cluster)

keep.samples <- y$samples$lib.size > 5e4 #5e4
table(keep.samples)
# FALSE  TRUE 
# 2    23

y <- y[, keep.samples]

# We then filter out lowly expressed genes.
# Default values: min.count=10, min.total.count=20
keep.genes <- filterByExpr(y, group=y$samples$cluster, min.count=10, min.total.count=20)

table(keep.genes)
# FALSE  TRUE 
# 29515  7265

y <- y[keep.genes, , keep=FALSE]

# TMM normalization is performed to estimate effective library sizes.
y <- normLibSizes(y)

summary(y$samples$norm.factors)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.8184  0.9479  0.9721  1.0074  1.0226  1.4377

cluster <- as.factor(y$samples$cluster)
levels(cluster)
cluster <- factor(cluster, levels = c('Primary FSGS', 'Maladaptive FSGS', 'PLA2R+ MN', 'Healthy control'))

pdf(paste0(projectFolder, ChosenCellGroup, "_EdgeR_MDS.pdf"))
plotMDS(y, pch=16, col=c(2:8)[cluster], main="MDS")
legend("bottomright", legend=paste0("cluster",levels(cluster)),pch=16, col=2:8, cex=0.8)
dev.off()

# To perform differential expression analysis between cell clusters, we create a design matrix.
design <- model.matrix(object = ~cluster)

colnames(design)[1] <- "Int"
head(design)
dim(design) #23 4

# The NB dispersion can be estimated using the estimateDisp function and visualized with plotBCV.
y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion

pdf(paste0(projectFolder, ChosenCellGroup, "_EdgeR_plotBCV_dispersion.pdf"))
plotBCV(y)
dev.off()

fit <- glmQLFit(y, design, robust=TRUE)

pdf(paste0(projectFolder, ChosenCellGroup, "_EdgeR_plotQLDispFit.pdf"))
plotQLDisp(fit)
dev.off()

# To confirm the identities of cell clusters, we perform differential expression analysis to identify marker genes of each cluster. 
# In particular, we compare each cluster with all the other clusters. 
# We construct a contrast matrix as follows so that each column of the contrast matrix represents a testing contrast for one cell cluster.
ncls <- nlevels(cluster)
contr <- rbind( matrix(1/(1-ncls), ncls, ncls), matrix(0, ncol(design)-ncls, ncls) )
diag(contr) <- 1 # diagonal gets 1
contr[1,] <- 0 # first row gets 0
rownames(contr) <- colnames(design)
colnames(contr) <- paste0("cluster", levels(cluster))
contr

#                                 clusterPrimary FSGS clusterMaladaptive FSGS clusterPLA2R+ MN clusterHealthy control
# Int                               0.0000000               0.0000000        0.0000000              0.0000000
# clusterMaladaptive FSGS          -0.3333333               1.0000000       -0.3333333             -0.3333333
# clusterPLA2R+ MN                 -0.3333333              -0.3333333        1.0000000             -0.3333333
# clusterHealthy control           -0.3333333              -0.3333333       -0.3333333              1.0000000

# We then perform quasi-likelihood F-test for each testing contrast. (glm=general linear models)
# The results are stored as a list of DGELRT objects, one for each comparison.
qlf <- list()
for(i in 1:ncls){
  qlf[[i]] <- glmQLFTest(fit, contrast=contr[,i])
  qlf[[i]]$comparison <- paste0("cluster", levels(cluster)[i], "_vs_others") 
}

#The top most significant DE genes of cluster 0 vs other clusters can be examined with topTags.
topTags(qlf[[1]], n=10L) #clusterPrimary FSGS_vs_others
topTags(qlf[[2]], n=10L) #clusterMaladaptive FSGS_vs_others 
topTags(qlf[[3]], n=10L) #clusterPLA2R+ MN_vs_others 
topTags(qlf[[4]], n=10L) #clusterHealthy control_vs_others

#The numbers of DE genes under each comparison are shown below.
dt <- lapply(lapply(qlf, decideTestsDGE), summary)
dt.all <- do.call("cbind", dt)
dt.all

#               clusterPrimary FSGS_vs_others         clusterMaladaptive FSGS_vs_others   clusterPLA2R+ MN_vs_others        clusterHealthy control_vs_others
# Down                               2                                 0                          0                                1
# NotSig                          7257                              7264                       7264                             7231
# Up                                 6                                 1                          1                               33

qlf_table1 <- qlf[[1]]$table

# A heat map is produced to visualize the top marker genes across all the pseudo-bulk samples.
lcpm <- edgeR::cpm(y, log=TRUE)
annot <- data.frame(cluster=paste0("cluster ", cluster))
rownames(annot) <- colnames(y)
ann_colors <- list(cluster=2:5)
names(ann_colors$cluster) <- paste0("cluster ", levels(cluster))

# Save the objects in a single file
save(lcpm, annot, ann_colors, file = paste0(projectFolder, "EdgeR_DEG_analysis_conditions_lcpm.RData"))


# Using TopTags to identify top DEGs and save datasheets

# Extract clusterPrimary FSGS_vs_others
topPrimary <- topTags(qlf[[1]], n=nrow(qlf[[1]]))$table
topPrimary <- topPrimary[order(topPrimary$FDR),]
write.csv2(topPrimary, file = paste0(projectFolder, ChosenCellGroup,"_EdgeR_topTags_primary_vs_others.csv"), row.names = F)

# Extract clusterMaladaptive control_vs_others
topMaladaptive <- topTags(qlf[[2]], n=nrow(qlf[[2]]))$table 
topMaladaptive <- topMaladaptive[order(topMaladaptive$FDR),]
write.csv2(topMaladaptive, file = paste0(projectFolder, ChosenCellGroup,"_EdgeR_topTags_maladaptive_vs_others.csv"), row.names = F)

# Extract clusterPLA2R+ MN_vs_others 
topPLA2R <- topTags(qlf[[3]], n=nrow(qlf[[3]]))$table
topPLA2R <- topPLA2R[order(topPLA2R$FDR),]
write.csv2(topPLA2R, file = paste0(projectFolder, ChosenCellGroup,"_EdgeR_topTags_PLA2R_vs_others.csv"), row.names = F)

# Extract clusterHealthy control_vs_others
topHealthy <- topTags(qlf[[4]], n=nrow(qlf[[4]]))$table
topHealthy <- topHealthy[order(topHealthy$FDR),] #ordering according to increasing FDR
write.csv2(topHealthy, file = paste0(projectFolder, ChosenCellGroup,"_EdgeR_topTags_healthy_vs_others.csv"), row.names = F)



###### 4.2: Heatmap: top 20 significant genes per condition - log2FC up > 0.5 and down < -0.5 - FDR < 0.10 - sorted on (unadj) PValue (output included in manuscript) ###### 

# Aim: to select the top 20 most significantly differentially expressed genes per condition. 
# Here, we applied a logFC cut-off AND a FDR cut-off AND sort on unadjusted P-value  (final sorting based upon unadjusted P-value and not FDR because there were some draws between genes using FDR), 
# We want to visualize the clinically relevant genes, either up or down (i.e. > 0.5 or < -0.5).
# Here, we are stringent (i.e., FDR of < 0.10).

# Using TopTags to identify top DEGs and visualize the heatmap

# Set the settings

chosenLog2FC_up <- 0.5
chosenLog2FC_up_filename <- "050"

chosenLog2FC_down <- -0.5
chosenLog2FC_down_filename <- "-050"

chosenFDR <- 0.10
chosenFDR_filename <- "010"

NumberTopGenes <- 20

# Extract clusterPrimary FSGS_vs_others and filter
topPrimary <- read.csv2(file = paste0(projectFolder, ChosenCellGroup,"_EdgeR_topTags_primary_vs_others.csv"))
topPrimary <- topPrimary %>% filter(logFC < -0.5 | logFC > 0.5)
topPrimary <- topPrimary[topPrimary$FDR < chosenFDR,]
topPrimary <- topPrimary[order(topPrimary$PValue),]
topPrimary <- topPrimary$gene # less than 20 genes remain

# Extract clusterMaladaptive control_vs_others and filter
topMaladaptive <- read.csv2(file = paste0(projectFolder, ChosenCellGroup,"_EdgeR_topTags_maladaptive_vs_others.csv"))
topMaladaptive <- topMaladaptive %>% filter(logFC < -0.5 | logFC > 0.5)
topMaladaptive <- topMaladaptive[topMaladaptive$FDR < chosenFDR,]
topMaladaptive <- topMaladaptive[order(topMaladaptive$PValue),]
topMaladaptive <- topMaladaptive$gene

# Extract clusterPLA2R+ MN_vs_others and filter
topPLA2R <- read.csv2(file = paste0(projectFolder, ChosenCellGroup,"_EdgeR_topTags_PLA2R_vs_others.csv"))
topPLA2R <- topPLA2R %>% filter(logFC < -0.5 | logFC > 0.5)
topPLA2R <- topPLA2R[topPLA2R$FDR < chosenFDR,]
topPLA2R <- topPLA2R[order(topPLA2R$PValue),]
topPLA2R <- topPLA2R$gene

# Extract clusterHealthy control_vs_others and filter
topHealthy <- read.csv2(file = paste0(projectFolder, ChosenCellGroup,"_EdgeR_topTags_healthy_vs_others.csv"))
topHealthy <- topHealthy %>% filter(logFC < -0.5 | logFC > 0.5)
topHealthy <- topHealthy[topHealthy$FDR < chosenFDR,]
topHealthy <- topHealthy[order(topHealthy$PValue),]
topHealthy <- topHealthy[1:NumberTopGenes,]$gene

topAll <- unique(c(topPrimary, topMaladaptive, topPLA2R, topHealthy))

topAll
# [1] "ECM1"       "DMD"        "MPP3"       "EHD4"       "AC004879.1" "SEMA4G"     "AC073359.2" "ZNF44"      "ZNF250"     "EPN2"       "SPSB1"      "AC011239.1" "MAGI3"      "PHLDB1"    
# [15] "LRP1B"      "TNFRSF19"   "WASHC1"     "SFMBT2"     "DLG5"       "LPIN1"      "PLEKHA2"    "ZBTB16"     "CEMIP2"     "FOXO3"      "PTPN13"     "FKBP5"      "IP6K3"      "KLF9"      
# [29] "PPM1L"      "LNX2"       "BTBD9"      "FBXO21"     "XYLT1"      "INSR"       "MAMLD1" 


# First re-organize the columns (i.e. samples) of heatmap
colnames(lcpm)

# Specify the new order of columns
new_column_order <- c(1,4,5,7,9,10,13,14,16,2,3,6,8,11,12,15,19,21,23,17,18,20,22) # samples sorted according to diagnoses

# Rearrange the columns
lcpm <- lcpm[, new_column_order]
colnames(lcpm)

# Make heatmap of top genes per condition
pdf(paste0(projectFolder, ChosenCellGroup,"_forPublication_EdgeR_Pheatmap_topTags_allContrasts_withSampleNames_Log2FCup_", chosenLog2FC_up_filename, "_Log2FCdown_", chosenLog2FC_down_filename, "_FDR_", chosenFDR_filename, "_NumberTopGenes_", NumberTopGenes, "_sortedPvalue.pdf"), width = 7, height = 7)
pheatmap::pheatmap(lcpm[topAll, ], breaks=seq(-2,2,length.out=101),
                   color=colorRampPalette(c("blue","white","red"))(100), scale="row", 
                   cluster_cols= F, 
                   cluster_rows=T,
                   border_color="NA", 
                   fontsize_row=6,
                   treeheight_row=30, 
                   treeheight_col=30, 
                   cutree_cols=7,
                   clustering_method="ward.D2", show_colnames=T, fontsize_col = 6,
                   annotation_col=annot, annotation_colors=ann_colors,
                   angle_col = 45)
dev.off()


###### 4.3: Heatmap of selected genes/gene sets ###### 

######  4.3.1: Plotting genes of immune pathways (output included in manuscript)  ###### 

# REACTOME_INNATE_IMMUNE_SYSTEM
genes_REACTOME_INNATE_IMMUNE_SYSTEM <- c("PSMA3", "HLA-E", "HLA-B", "MYO1C", "TIMP2", "CD63", "B2M", 
                                         "NLRP1", "HLA-C", "KLRD1", "C1R", "ARPC5", "CD81", "CHI3L1", 
                                         "CAPZA2", "CMTM6", "LRRFIP1", "CST3", "HSP90B1", "PSMB1", "VCP", 
                                         "ATP6V0B", "C1S", "MAGT1", "ATG5", "ATP6V1D", "MOSPD2", "NDUFC2", 
                                         "NPC2", "PRKCSH", "PGAM1", "PRDX6", "SHC1", "PLD3", "LAMP2", 
                                         "CYSTM1", "UBB", "ATOX1", "CYFIP2", "PPP2CB", "DYNLT1", "BTRC", 
                                         "MIF", "JUP", "S100A11", "DNAJC3", "TMBIM1", "LAMP1", "PGRMC1", 
                                         "SURF4", "ACTG1", "SRP14", "SDCBP", "PSMA1", "ATP6AP2", "TANK", 
                                         "HSPA8", "TUBB4B", "CTSA", "ATP6V0D1", "CAP1", "RAB5C", "UBE2N", 
                                         "HRAS", "TRIM56", "DUSP3", "ATP6V1G1", "TXNDC5", "GNS", "CD46", 
                                         "PPP2R1A", "POLR2L", "NCSTN", "FUCA2", "CD55", "RAB31", "GSTP1", "PPIA", 
                                         "RAB7A", "PSMC2", "RAP1B", "PSMB5", "BST2", "VAV2", "MEF2A", "BCL2", 
                                         "PSAP", "AP1M1", "WASF2", "ATP6V0A1", "GSN", "MEF2C", "CKAP4", "UBC", 
                                         "ATF1", "RHOA", "TMEM179B", "ARSA", "CRK", "CD47", "DTX4", "ACTR2", 
                                         "DHX9", "FOS", "SLC15A4", "CD59", "TXNIP", "ITGAV", "NFKB1", "ACTB", 
                                         "CREBBP", "CAPN1", "KRAS", "RELB", "PIK3R4", "YPEL5", "SCAMP1", "ASAH1", 
                                         "FTH1", "LYN", "CYB5R3", "DNM1", "BRK1", "DYNLL1", "PRKDC", "RAB14", 
                                         "CYBA", "GRN", "SKP1", "PSME1", "OTUD5", "NFATC1", "GOLGA7", "SUGT1", 
                                         "UBA3", "NOD1", "UBR4", "PAFAH1B2", "ARPC1A", "VAMP8", "CSTB", "IGF2R", 
                                         "PDAP1", "CTSL", "FTL", "MAP2K1", "RPS6KA3")

# REACTOME_NEUTROPHIL_DEGRANULATION
genes_REACTOME_NEUTROPHIL_DEGRANULATION <-  c("HLA-B", "TIMP2", "CD63", "B2M", "HLA-C", "ARPC5", "CHI3L1", "CMTM6", "CST3", "PSMB1", 
                                              "VCP", "MAGT1", "ATP6V1D", "MOSPD2", "NDUFC2", "NPC2", "PGAM1", "PRDX6", "LAMP2", 
                                              "CYSTM1", "DYNLT1", "MIF", "JUP", "S100A11", "DNAJC3", "TMBIM1", "LAMP1", 
                                              "PGRMC1", "SURF4", "SRP14", "SDCBP", "ATP6AP2", "HSPA8", "TUBB4B", "CTSA",
                                              "CAP1", "RAB5C", "TXNDC5", "GNS", "NCSTN", "FUCA2", "CD55", "RAB31", "GSTP1", 
                                              "PPIA", "RAB7A", "PSMC2", "RAP1B", "BST2", "PSAP", "AP1M1", "ATP6V0A1", "GSN", 
                                              "CKAP4", "RHOA", "TMEM179B", "ARSA", "CD47", "ACTR2", "SLC15A4", "CD59", "ITGAV", 
                                              "NFKB1", "CAPN1", "YPEL5", "SCAMP1", "ASAH1", "FTH1", "CYB5R3", "DYNLL1", "RAB14", 
                                              "CYBA", "GRN", "GOLGA7", "UBR4", "PAFAH1B2", "VAMP8", "CSTB", "IGF2R", "PDAP1",
                                              "FTL", "PDXK", "TUBB", "FAF2", "TMEM30A", "HSP90AA1", "PSEN1", "HSPA1A", "MLEC")

# REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION
genes_REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION <- c("PSMA3", "HLA-E", "HLA-B", "HLA-A", "CALR", "B2M", "HLA-C", "PDIA3", "PSMB1", 
                                                          "UBB", "PSMA1", "PSMC2", "SEC61B", "PSMB5", "UBC", "TAPBP", "ITGAV", 
                                                          "CYBA", "PSME1", "VAMP8", "CTSL")

# REACTOME_ANTIGEN_PRESENTATION_FOLDING_ASSEMBLY_AND_PEPTIDE_LOADING_OF_CLASS_I_MHC
genes_REACTOME_ANTIGEN_PRESENTATION_FOLDING_ASSEMBLY_AND_PEPTIDE_LOADING_OF_CLASS_I_MHC <- c("HLA-E", "HLA-B", "HLA-A", "CALR", "B2M", "HLA-C", "HSPA5", "PDIA3", "CANX", "BCAP31", "ERAP2")

# REACTOME_INTERFERON_ALPHA_BETA_SIGNALING
genes_REACTOME_INTERFERON_ALPHA_BETA_SIGNALING <- c("HLA-E", "HLA-B", "HLA-A", "ADAR", "HLA-C", "IFI6", "UBB", "OAS1", "IRF9", "IFITM3", "STAT1", "KPNA1", "BST2", "UBC")

# REACTOME_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL
genes_REACTOME_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL <- c("HLA-E", "HLA-B", "HLA-A", "B2M", "HLA-C", "KLRD1", "CD81", "ITGB1", "CXADR", "CDH1")

# HALLMARK_COMPLEMENT
genes_HALLMARK_COMPLEMENT <- c("TIMP1", "TIMP2", "CA2", "C1R", "HSPA5", "C1S", "LAMP2", "ATOX1", "PPP2CB", "ANXA5", "GNAI2", "CTSO", "ERAP2", "S100A13", "CD46", "SPOCK2", "CD55", "FDX1", "USP15", "CSRP1", "GNB2", "CD59", "LIPA", "DPP4", "IRF1", "GPD2", "LYN")

# HALLMARK_INTERFERON_GAMMA_RESPONSE
genes_HALLMARK_INTERFERON_GAMMA_RESPONSE <- c("PSMA3", "HLA-B", "HLA-A", "ADAR", "B2M", "C1R", "C1S", "TNFAIP2", "SPPL2A", "IRF9", "LGALS3BP", "IL7", "IFITM3", "STAT1", "APOL6", "BST2", "EPSTI1", "TAPBP", "AUTS2", "TXNIP", "DDX60", "NFKB1", "CMTR1", "IRF1", "TRIM26", "PDE4B", "PSME1", "VAMP5", "NOD1", "PNPT1", "VAMP8", "RBCK1")

# HALLMARK_INTERFERON_ALPHA_RESPONSE
genes_HALLMARK_INTERFERON_ALPHA_RESPONSE <- c("PSMA3", "ADAR", "B2M", "HLA-C", "C1S", "OAS1", "IRF9", "LGALS3BP", "IL7", "IFITM3", "BST2", "EPSTI1", "TRIM5", "CD47", "TXNIP", "DDX60", "CMTR1", "IRF1", "TRIM26", "PSME1")

immune_genes <- unique(c(genes_REACTOME_INNATE_IMMUNE_SYSTEM,
                         genes_REACTOME_NEUTROPHIL_DEGRANULATION,
                         genes_REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION,
                         genes_REACTOME_ANTIGEN_PRESENTATION_FOLDING_ASSEMBLY_AND_PEPTIDE_LOADING_OF_CLASS_I_MHC,
                         genes_REACTOME_INTERFERON_ALPHA_BETA_SIGNALING,
                         genes_REACTOME_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL,
                         genes_HALLMARK_COMPLEMENT,
                         genes_HALLMARK_INTERFERON_GAMMA_RESPONSE,
                         genes_HALLMARK_INTERFERON_ALPHA_RESPONSE))

# Now let's check whether these genes are significantly upreg in primary FSGS PODs
chosenLog2FC_up <- 0.25
chosenPValue <- 0.05
topPrimary <- read.csv2(file = paste0("path/seurat_POD_HQ2/Diagnosis_analysis/EdgeR/seurat_POD_HQ2_EdgeR_topTags_primary_vs_others.csv"))
topPrimary <- topPrimary[topPrimary$logFC > chosenLog2FC_up,]
topPrimary <- topPrimary[topPrimary$PValue < chosenPValue,]
allGenes <- topPrimary$gene

ChosenGenes <- immune_genes[immune_genes %in% allGenes] # 45 genes are present

# Make heatmap of these selected genes
pdf(paste0(projectFolder, ChosenCellGroup,"_EdgeR_Pheatmap_immune-genes.pdf"), width = 8, height = 7)
pheatmap::pheatmap(lcpm[ChosenGenes, ], breaks=seq(-2,2,length.out=101),
                   color=colorRampPalette(c("blue","white","red"))(100), scale="row", 
                   cluster_cols=F, border_color="NA", fontsize_row=6,
                   treeheight_row=70, treeheight_col=70, cutree_cols=7,
                   clustering_method="ward.D2", show_colnames=T, fontsize_col = 6,
                   annotation_col=annot, annotation_colors=ann_colors,
                   angle_col = 45)
dev.off()


######  4.3.2: Plotting genes of mTORC1-signaling (output included in manuscript)  ###### 

# HALLMARK_MTORC1_SIGNALING
selected_genes <- c("PSMA3", "CALR", "SQSTM1", "ELOVL5", "HSPA5", "USO1", "RPN1", "HSP90B1", "ATP2A2", "BTG2", "PRDX1", "CANX",
                    "CYP51A1", "ATP6V1D", "ACSL3", "CD9", "SQLE", "HMGCS1", "GBE1", "RDH11", "UCHL5", "PPIA", "PSMC2", "PGK1",
                    "SKAP2", "PSMB5", "UFM1", "ACTR2", "LDHA", "SSR1", "GAPDH", "SERP1", "NMT1", "M6PR", "CCT6A", "TXNRD1",
                    "SEC11A", "PDAP1", "SLC6A6", 
                    "SREBF1", # added manually based on literature
                    "VEGFA" # added manually based on literature 
)

# Now let's check whether these genes are significantly upreg in primary FSGS PODs
chosenLog2FC_up <- 0.25
chosenPValue <- 0.05
topPrimary <- read.csv2(file = paste0("path/seurat_POD_HQ2/Diagnosis_analysis/EdgeR/seurat_POD_HQ2_EdgeR_topTags_primary_vs_others.csv"))
topPrimary <- topPrimary[topPrimary$logFC > chosenLog2FC_up,]
topPrimary <- topPrimary[topPrimary$PValue < chosenPValue,]
allGenes <- topPrimary$gene

ChosenGenes <- selected_genes[selected_genes %in% allGenes] # 15 genes are present

# Make heatmapof these selected genes
pdf(paste0(projectFolder, ChosenCellGroup,"_EdgeR_Pheatmap_mTORC1-genes.pdf"), width = 8, height = 4)
pheatmap::pheatmap(lcpm[ChosenGenes, ], breaks=seq(-2,2,length.out=101),
                   color=colorRampPalette(c("blue","white","red"))(100), scale="row", 
                   cluster_cols=F, border_color="NA", fontsize_row=6,
                   treeheight_row=30, treeheight_col=70, cutree_cols=7,
                   clustering_method="ward.D2", show_colnames=T, fontsize_col = 6,
                   annotation_col=annot, annotation_colors=ann_colors,
                   angle_col = 45)
dev.off()



######  4.3.3: Plotting ECM genes (output included in manuscript)  ###### 

# We see REACTOME_ECM organization in fgsea for primary FSGS AND maladaptive FSGS podocytes.
# Let's plot the significant genes from the leading edge of both.

leadingEdgePrimary <- c("TIMP1",
                        "TIMP2",
                        "DAG1",
                        "CD151",
                        "BMP2",
                        "PPIB",
                        "ITGB1",
                        # "HSPG2", # also present in maladaptive FSGS
                        "CAPN12",
                        "P4HB",
                        "BGN",
                        "HTRA1",
                        "P3H2",
                        "PLEC",
                        "NCAM1",
                        "SPARC",
                        "CAST",
                        "LAMB1",
                        "LOX",
                        "ADAM17",
                        "CDH1",
                        "DCN",
                        "ITGB3",
                        "COL4A4",
                        "AGRN",
                        "NCSTN",
                        "COL7A1",
                        "LAMB2",
                        "LTBP4",
                        "CD47",
                        "SDC4",
                        "ITGAV",
                        "CAPN1")

# Check whether these genes are significantly upreg in primary FSGS PODs
chosenLog2FC_up <- 0.25
chosenPValue <- 0.05
topPrimary <- read.csv2(file = "path/seurat_POD_HQ2/Diagnosis_analysis/EdgeR/seurat_POD_HQ2_EdgeR_topTags_primary_vs_others.csv")
topPrimary <- topPrimary[topPrimary$logFC > chosenLog2FC_up,]
topPrimary <- topPrimary[topPrimary$PValue < chosenPValue,]
allGenes <- topPrimary$gene

ChosenGenes_primary <- leadingEdgePrimary[leadingEdgePrimary %in% allGenes] # 7 genes are present

# Leading edge from maladaptive vs. others: following significant pathways:  
# REACTOME_NON_INTEGRIN_MEMBRANE_ECM_INTERACTIONS
# REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION
# REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS

leadingEdgeMaladaptive <- c("HSPG2",
                            "DMD", 
                            "COL1A2", 
                            "COL3A1", 
                            "THBS1", 
                            "COL4A1",
                            "COL4A2",
                            "NTN4",
                            "PRKCA",
                            "COL5A2",
                            "LAMB1",
                            "AGRN",
                            "COLGALT2", 
                            "CDH1", 
                            "CD44", 
                            "ITGA1",
                            "P4HA1",
                            "LTBP3",
                            "BMP2",
                            "KLK7",
                            "EFEMP1",
                            "FBN1",
                            "PTPRS",
                            "ADAM17",
                            "CTSD",
                            "LAMA3",
                            "MMP16",
                            "DCN",
                            "LOX",
                            "P3H2",
                            "CTSL",
                            "ITGA9",
                            "BGN",
                            "ITGA3",
                            "BSG",
                            "LTBP1",
                            "PHYKPL",
                            "APP",
                            "COL27A1",
                            "CD151",
                            "HTRA1",
                            "CAPN1",
                            "DDR1",
                            "ACTN1",
                            "LAMB2",
                            "SDC4",
                            "PPIB",
                            "SDC2",
                            "CAPN7") 

# Check whether these genes are significantly upreg in maladaptive FSGS PECs
chosenLog2FC_up <- 0.25
chosenPValue <- 0.05
topMaladaptive <- read.csv2(file = "path/seurat_POD_HQ2/Diagnosis_analysis/EdgeR/seurat_POD_HQ2_EdgeR_topTags_maladaptive_vs_others.csv")
topMaladaptive <- topMaladaptive[topMaladaptive$logFC > chosenLog2FC_up,]
topMaladaptive <- topMaladaptive[topMaladaptive$PValue < chosenPValue,]
allGenes <- topMaladaptive$gene

ChosenGenes_maladaptive <- leadingEdgeMaladaptive[leadingEdgeMaladaptive %in% allGenes] # 9 genes are present

# Combination
ChosenGenes <- unique(c(ChosenGenes_primary, ChosenGenes_maladaptive))

# Plot heatmap (output included in manuscript)
pdf(paste0(projectFolder, ChosenCellGroup,"_EdgeR_Pheatmap_leadingEdge_REACTOME_ECM-genes_noReclustering.pdf"), width = 6, height = 3.5)
pheatmap::pheatmap(lcpm[ChosenGenes, ], breaks=seq(-2,2,length.out=101),
                   color=colorRampPalette(c("blue","white","red"))(100), scale="row", 
                   cluster_cols=F, border_color="NA", fontsize_row=6,
                   cluster_rows = F,
                   treeheight_row=30, treeheight_col=30, cutree_cols=7,
                   clustering_method="ward.D2", show_colnames=T, fontsize_col = 6,
                   annotation_col=annot, annotation_colors=ann_colors,
                   angle_col = 45)
dev.off()


###### 4.4: Volcano plots of DGE per condition ###### 

# clusterPrimary FSGS_vs_others
topPrimary <- read.csv2(file = paste0(projectFolder, ChosenCellGroup,"_EdgeR_topTags_primary_vs_others.csv"))

# Adjusted P-val (FDR)

chosenLog2FC_up <- 0.5
chosenLog2FC_up_filename <- "050"

chosenFDR <- 0.10
chosenFDR_filename <- "010"

pdf(paste0(projectFolder, ChosenCellGroup,"_EdgeR_VolcanoPlot_topTags_primary_vs_others_chosenLog2FCup_", chosenLog2FC_up_filename, "_chosenFDRvalue_", chosenFDR_filename, ".pdf"), width = 5, height = 5.5)
EnhancedVolcano(topPrimary,
                lab = topPrimary$gene,
                x = 'logFC',
                y = 'FDR',
                title = 'Primary vs. all others',
                subtitle = 'EdgeR DEG analysis, adjusted P-value (FDR)',
                xlim = c(min(topPrimary$logFC) * 1.2, max(topPrimary$logFC) * 1.2),
                xlab = bquote('Log2(fold change)'),
                pCutoff = chosenFDR,
                FCcutoff = chosenLog2FC_up,
                pointSize = 2,
                labSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                colConnectors = 'grey30',
                max.overlaps = Inf,
                col = c("grey30", "royalblue", "royalblue", "red2")) +
  scale_y_continuous(breaks = seq(0, 10, by = 1)) +
  scale_x_continuous(breaks = seq(-10, 10, by = 1))
dev.off()

# Analogous analyses for "maladaptive vs. others", "PLA2R+ MN vs. others", "healthy vs. others"


###### 4.5: GSEA with fsgea - Primary vs. others ###### 

# Load packages
library(tidyverse)
library(RColorBrewer)
library(fgsea)
library(ggplot2)

# Load previous DEG-files with log2FC and P-values (here from EdgeR)
topPrimary <- read.csv2(file = paste0("path/seurat_POD_HQ2/Diagnosis_analysis/EdgeR/seurat_POD_HQ2_EdgeR_topTags_primary_vs_others.csv"))
topMaladaptive <- read.csv2(file = paste0("path/seurat_POD_HQ2/Diagnosis_analysis/EdgeR/seurat_POD_HQ2_EdgeR_topTags_maladaptive_vs_others.csv"))
topPLA2R <- read.csv2(file = paste0("path/seurat_POD_HQ2/Diagnosis_analysis/EdgeR/seurat_POD_HQ2_EdgeR_topTags_PLA2R_vs_others.csv"))
topHealthy <- read.csv2(file = paste0("path/seurat_POD_HQ2/Diagnosis_analysis/EdgeR/seurat_POD_HQ2_EdgeR_topTags_healthy_vs_others.csv"))

# Set seed
set.seed(123456)

# Set path for the GSEA gene sets
bg_path <- "path/GSEA_GeneSets/"

# Load two custom functions
## Function: Adjacency matrix to list 
matrix_to_list <- function(pws){
  pws.l <- list()
  for (pw in colnames(pws)) {
    pws.l[[pw]] <- rownames(pws)[as.logical(pws[, pw])]
  }
  return(pws.l)
}

## Function: prepare_gmt 
prepare_gmt <- function(gmt_file, genes_in_data, savefile = FALSE){
  # for debug
  #file <- gmt_files[1]
  #genes_in_data <- df$gene_symbol
  
  # Read in gmt file
  gmt <- gmtPathways(gmt_file)
  hidden <- unique(unlist(gmt))
  
  # Convert gmt file to a matrix with the genes as rows and for each go annotation (columns) the values are 0 or 1
  mat <- matrix(NA, dimnames = list(hidden, names(gmt)),
                nrow = length(hidden), ncol = length(gmt))
  for (i in 1:dim(mat)[2]){
    mat[,i] <- as.numeric(hidden %in% gmt[[i]])
  }
  
  #Subset to the genes that are present in our data to avoid bias
  hidden1 <- intersect(genes_in_data, hidden)
  mat <- mat[hidden1, colnames(mat)[which(colSums(mat[hidden1,])>5)]] # filter for gene sets with more than 5 genes annotated
  # And get the list again
  final_list <- matrix_to_list(mat) # for this we use the function we previously defined
  
  if(savefile){
    saveRDS(final_list, file = paste0(gsub('.gmt', '', gmt_file), '_subset_', format(Sys.time(), '%d%m'), '.RData'))
  }
  
  print('.gmt conversion successfull')
  return(final_list)
}

# Analysis

# Set the comparison and read in the data:
projectFolder <- "path/seurat_POD_HQ2/Diagnosis_analysis/EdgeR/fgsea/primary_vs_others/"
name_of_comparison <- 'primary_vs_others'
df <- topPrimary
head(df)

### 1. Rank the DEG-list
rankings <- sign(df$logFC)*(-log10(df$PValue))
names(rankings) <- df$gene # genes as names#

head(rankings)

rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking
plot(rankings)

max(rankings)
min(rankings)

## 2. Prepare background genes
# For GSEA
# Filter out the gmt files for KEGG, Reactome and GOBP
my_genes <- df$gene

list.files(bg_path)
gmt_files <- list.files(path = bg_path, pattern = '.gmt', full.names = TRUE)
gmt_files

######### REACTOME #

# Choose the Gene set
background_genes <- 'REACTOME'
bg_genes <- prepare_gmt(gmt_files[2], my_genes, savefile = FALSE) # for C2 REACTOME

## 3. Run GSEA 
GSEAres <- fgsea(pathways = bg_genes, # List of gene sets to check
                 stats = rankings,
                 scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                 minSize = 10,
                 maxSize = 500,
                 nproc = 1) # for parallelisation

GSEAres <- GSEAres[order(padj), ]

## 4. Save the results 
filename <- paste0(projectFolder, 'seurat_POD_HQ2_GSEA_', name_of_comparison, '_', background_genes, '_uncollapsed') 
data.table::fwrite(GSEAres, file = paste0(filename, '_gsea_results.csv'), sep = "\t", sep2 = c("", " ", ""))

## 5. Check results
# Top enriched pathways (ordered by p-val)
head(GSEAres[order(padj), ])

sum(GSEAres[, padj < 0.05])

signPathwaysUP <- GSEAres[GSEAres$ES > 0 & GSEAres$padj < 0.05,]
signPathwaysUP <- signPathwaysUP[order(signPathwaysUP$padj),]     # order on increasing Padj            
signPathwaysUP <- signPathwaysUP$pathway

signPathwaysDOWN <- GSEAres[GSEAres$ES < 0 & GSEAres$padj < 0.05,]
signPathwaysDOWN <- signPathwaysDOWN[order(signPathwaysDOWN$padj),]     # order on increasing Padj            
signPathwaysDOWN <- signPathwaysDOWN$pathway

signPathways <- c(signPathwaysUP, rev(signPathwaysDOWN))

pdf(file = paste0(projectFolder, 'seurat_POD_HQ2_GSEA_', name_of_comparison, '_', background_genes, '_signPathways.pdf'), width = 20, height = 15)
plotGseaTable(bg_genes[signPathways], 
              stats = rankings, 
              fgseaRes = GSEAres, 
              gseaParam = 0.5)
dev.off()

# Select only independent pathways, removing redundancies/similar pathways
collapsedPathways <- collapsePathways(GSEAres[order(padj)][padj < 0.05],
                                      bg_genes,
                                      rankings)

# Write table of significant collapsed pathways
TopCollapsedPathways <- GSEAres[pathway %in% collapsedPathways$mainPathways][order(padj),]

filename <- paste0(projectFolder, 'seurat_POD_HQ2_GSEA_', name_of_comparison, '_', background_genes, '_collapsed') 
data.table::fwrite(TopCollapsedPathways, file = paste0(filename, '_gsea_results.csv'), sep = "\t", sep2 = c("", " ", ""))

# Make fgsea plot of significant collapsed pathways
mainPathways <- GSEAres[pathway %in% collapsedPathways$mainPathways][order(-NES), pathway]

pdf(file = paste0(projectFolder, 'seurat_POD_HQ2_GSEA_', name_of_comparison, '_', background_genes, '_CollapsedPathways.pdf'), width = 20, height = 15)
plotGseaTable(bg_genes[mainPathways], 
              rankings, 
              GSEAres, 
              gseaParam = 0.5)
dev.off()

# Plot the NES of the significant collapsed pathways
TopCollapsedPathways <- GSEAres[pathway %in% collapsedPathways$mainPathways][order(-NES),]
TopCollapsedPathwaysLevels <- TopCollapsedPathways[order(NES), pathway]
TopCollapsedPathways$pathway <- factor(TopCollapsedPathways$pathway, levels = TopCollapsedPathwaysLevels)
levels(TopCollapsedPathways$pathway)

pdf(file = paste0(projectFolder, 'seurat_POD_HQ2_GSEA_', name_of_comparison, '_', background_genes, '_CollapsedPathways_NES_barchart.pdf'), width = 12, height = 10)
gsea_NES_barplot <- ggplot(TopCollapsedPathways, aes(y=NES, x=pathway, fill = padj)) +
  geom_bar(stat= "identity", position="dodge") +
  scale_fill_gradient(low = "#76C6FB", high = "black") +
  ylab("Normalized Enrichment Score") +
  xlab("Pathway") +
  labs(title= paste0("Podocytes: ", name_of_comparison), subtitle= paste0(background_genes, " pathway")) +
  theme_bw() +
  theme(axis.title.x = element_text(size=10), axis.text.x = element_text(size=10, angle = 0, hjust = 0.5, vjust=0.5, colour="black"),
        axis.title.y = element_text(size=10), axis.text.y = element_text(size=6, angle = 0, hjust = 1, vjust=0.5, colour="black"),
        plot.title = element_text(hjust=0.5),
        plot.subtitle = element_text(hjust=0.5, size=10, face="italic"),
        axis.ticks = element_blank(),
        legend.key = element_rect(fill = "white"),
        legend.title=element_text(size=10),
        legend.text = element_text(size=8),
        legend.position="right") + coord_flip()
gsea_NES_barplot
dev.off()

# plot the most significantly enriched pathway
pdf(file = paste0(projectFolder, 'seurat_POD_HQ2_GSEA_', name_of_comparison, '_', background_genes, '_EnrichmentPlot_topPathway.pdf'), width = 10, height = 10)
plotEnrichment(bg_genes[[head(GSEAres[order(padj), ], 1)$pathway]],
               rankings) +
  labs(title = head(GSEAres[order(padj), ], 1)$pathway)
dev.off()

# plot for publication
TopCollapsedPathways <- GSEAres[pathway %in% collapsedPathways$mainPathways][order(padj),] # order on most significant (collapsed) pathways
TopCollapsedPathways <- TopCollapsedPathways[1:15,] # select top 15
TopCollapsedPathwaysLevels <- TopCollapsedPathways[order(NES), pathway]
TopCollapsedPathways$pathway <- factor(TopCollapsedPathways$pathway, levels = TopCollapsedPathwaysLevels)
levels(TopCollapsedPathways$pathway)

pdf(file = paste0(projectFolder, 'forPublication_seurat_POD_HQ2_GSEA_', name_of_comparison, '_', background_genes, '_top15_significant_CollapsedPathways_NES_barchart.pdf'), width = 10, height = 5)
ggplot(TopCollapsedPathways, aes(y=NES, x=pathway, fill = padj)) +
  geom_bar(stat= "identity", position="dodge") +
  scale_fill_gradient(low = "#76C6FB", high = "black") +
  ylab("Normalized Enrichment Score") +
  xlab("Pathway") +
  labs(title= paste0("Podocytes: ", name_of_comparison), subtitle= paste0(background_genes, " pathway")) +
  theme_bw() +
  theme(axis.title.x = element_text(size=10), axis.text.x = element_text(size=10, angle = 0, hjust = 0.5, vjust=0.5, colour="black"),
        axis.title.y = element_text(size=10), axis.text.y = element_text(size=6, angle = 0, hjust = 1, vjust=0.5, colour="black"),
        plot.title = element_text(hjust=0.5),
        plot.subtitle = element_text(hjust=0.5, size=10, face="italic"),
        axis.ticks = element_blank(),
        legend.key = element_rect(fill = "white"),
        legend.title=element_text(size=10),
        legend.text = element_text(size=8),
        legend.position="right") + coord_flip()
dev.off()

######### Hallmark #

# Choose the Gene set
background_genes <- 'Hallmark'
bg_genes <- prepare_gmt(gmt_files[4], my_genes, savefile = FALSE) # for H Hallmark

## 3. Run GSEA 
GSEAres <- fgsea(pathways = bg_genes, # List of gene sets to check
                 stats = rankings,
                 scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                 minSize = 10,
                 maxSize = 500,
                 nproc = 1) # for parallelisation

GSEAres <- GSEAres[order(padj), ]

## 4. Save the results 
filename <- paste0(projectFolder, 'seurat_POD_HQ2_GSEA_', name_of_comparison, '_', background_genes, '_uncollapsed') 
data.table::fwrite(GSEAres, file = paste0(filename, '_gsea_results.csv'), sep = "\t", sep2 = c("", " ", ""))

## 5. Check results
# Top enriched pathways (ordered by p-val)
head(GSEAres[order(padj), ])

sum(GSEAres[, padj < 0.05])

signPathwaysUP <- GSEAres[GSEAres$ES > 0 & GSEAres$padj < 0.05,]
signPathwaysUP <- signPathwaysUP[order(signPathwaysUP$padj),]     # order on increasing Padj            
signPathwaysUP <- signPathwaysUP$pathway

signPathwaysDOWN <- GSEAres[GSEAres$ES < 0 & GSEAres$padj < 0.05,]
signPathwaysDOWN <- signPathwaysDOWN[order(signPathwaysDOWN$padj),]     # order on increasing Padj            
signPathwaysDOWN <- signPathwaysDOWN$pathway

signPathways <- c(signPathwaysUP, rev(signPathwaysDOWN))

pdf(file = paste0(projectFolder, 'seurat_POD_HQ2_GSEA_', name_of_comparison, '_', background_genes, '_signPathways.pdf'), width = 20, height = 15)
plotGseaTable(bg_genes[signPathways], 
              stats = rankings, 
              fgseaRes = GSEAres, 
              gseaParam = 0.5)
dev.off()

# Select only independent pathways, removing redundancies/similar pathways
collapsedPathways <- collapsePathways(GSEAres[order(padj)][padj < 0.05],
                                      bg_genes,
                                      rankings)

# Write table of significant collapsed pathways
TopCollapsedPathways <- GSEAres[pathway %in% collapsedPathways$mainPathways][order(padj),]

filename <- paste0(projectFolder, 'seurat_POD_HQ2_GSEA_', name_of_comparison, '_', background_genes, '_collapsed') 
data.table::fwrite(TopCollapsedPathways, file = paste0(filename, '_gsea_results.csv'), sep = "\t", sep2 = c("", " ", ""))

# Make fgsea plot of significant collapsed pathways
mainPathways <- GSEAres[pathway %in% collapsedPathways$mainPathways][order(-NES), pathway]

pdf(file = paste0(projectFolder, 'seurat_POD_HQ2_GSEA_', name_of_comparison, '_', background_genes, '_CollapsedPathways.pdf'), width = 20, height = 15)
plotGseaTable(bg_genes[mainPathways], 
              rankings, 
              GSEAres, 
              gseaParam = 0.5)
dev.off()

# Plot the NES of the significant collapsed pathways
TopCollapsedPathways <- GSEAres[pathway %in% collapsedPathways$mainPathways][order(-NES),]
TopCollapsedPathwaysLevels <- TopCollapsedPathways[order(NES), pathway]
TopCollapsedPathways$pathway <- factor(TopCollapsedPathways$pathway, levels = TopCollapsedPathwaysLevels)
levels(TopCollapsedPathways$pathway)

pdf(file = paste0(projectFolder, 'seurat_POD_HQ2_GSEA_', name_of_comparison, '_', background_genes, '_CollapsedPathways_NES_barchart.pdf'), width = 6, height = 8)
gsea_NES_barplot <- ggplot(TopCollapsedPathways, aes(y=NES, x=pathway, fill = padj)) +
  geom_bar(stat= "identity", position="dodge") +
  scale_fill_gradient(low = "#76C6FB", high = "black") +
  ylab("Normalized Enrichment Score") +
  xlab("Pathway") +
  labs(title= paste0("Podocytes: ", name_of_comparison), subtitle= paste0(background_genes, " pathway")) +
  theme_bw() +
  theme(axis.title.x = element_text(size=10), axis.text.x = element_text(size=10, angle = 0, hjust = 0.5, vjust=0.5, colour="black"),
        axis.title.y = element_text(size=10), axis.text.y = element_text(size=6, angle = 0, hjust = 1, vjust=0.5, colour="black"),
        plot.title = element_text(hjust=0.5),
        plot.subtitle = element_text(hjust=0.5, size=10, face="italic"),
        axis.ticks = element_blank(),
        legend.key = element_rect(fill = "white"),
        legend.title=element_text(size=10),
        legend.text = element_text(size=8),
        legend.position="right") + 
  coord_flip()
gsea_NES_barplot
dev.off()

# plot the most significantly enriched pathway
pdf(file = paste0(projectFolder, 'seurat_POD_HQ2_GSEA_', name_of_comparison, '_', background_genes, '_EnrichmentPlot_topPathway.pdf'), width = 10, height = 10)
plotEnrichment(bg_genes[[head(GSEAres[order(padj), ], 1)$pathway]],
               rankings) +
  labs(title = head(GSEAres[order(padj), ], 1)$pathway)
dev.off()

# plot for publication
TopCollapsedPathways <- GSEAres[pathway %in% collapsedPathways$mainPathways][order(padj),] # order on most significant (collapsed) pathways
TopCollapsedPathways <- TopCollapsedPathways[1:15,] # select top 15
TopCollapsedPathwaysLevels <- TopCollapsedPathways[order(NES), pathway]
TopCollapsedPathways$pathway <- factor(TopCollapsedPathways$pathway, levels = TopCollapsedPathwaysLevels)
levels(TopCollapsedPathways$pathway)

pdf(file = paste0(projectFolder, 'forPublication_seurat_POD_HQ2_GSEA_', name_of_comparison, '_', background_genes, '_top15_significant_CollapsedPathways_NES_barchart.pdf'), width = 5.5, height = 5)
ggplot(TopCollapsedPathways, aes(y=NES, x=pathway, fill = padj)) +
  geom_bar(stat= "identity", position="dodge") +
  scale_fill_gradient(low = "#76C6FB", high = "black") +
  ylab("Normalized Enrichment Score") +
  xlab("Pathway") +
  labs(title= paste0("Podocytes: ", name_of_comparison), subtitle= paste0(background_genes, " pathway")) +
  theme_bw() +
  theme(axis.title.x = element_text(size=10), axis.text.x = element_text(size=10, angle = 0, hjust = 0.5, vjust=0.5, colour="black"),
        axis.title.y = element_text(size=10), axis.text.y = element_text(size=6, angle = 0, hjust = 1, vjust=0.5, colour="black"),
        plot.title = element_text(hjust=0.5),
        plot.subtitle = element_text(hjust=0.5, size=10, face="italic"),
        axis.ticks = element_blank(),
        legend.key = element_rect(fill = "white"),
        legend.title=element_text(size=10),
        legend.text = element_text(size=8),
        legend.position="right") + coord_flip()
dev.off()

# Analogous analyses for "maladaptive vs. others", "PLA2R+ MN vs. others", "healthy vs. others"



##### 5. PART 5 - DGE ANALYSIS OF ONE GROUP VS. CONTROLS - EdgeR (pseudobulk)  #####

# Switch to R-version with installation of EdgeR

# Step 1: run R 4.3
library(Seurat)
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(edgeR)
library(EnhancedVolcano)

###### 5.1: Run general EdgeR DEG analysis ###### 

projectFolder <- "path/seurat_POD_HQ2/Diagnosis_analysis_primary_vs_control/EdgeR/"

# Read in the dataset
seurat_POD_HQ2 <- readRDS(file = "path/seurat_POD_HQ2/seurat_POD_HQ2_harmony_annot.Rds")
seurat_obj <- seurat_POD_HQ2 #insert your object here
ChosenCellGroup <- "seurat_POD_HQ2"

#check conditions
table(seurat_obj$Diagnosis)
table(seurat_obj$orig.ident,seurat_obj$Diagnosis)

# Exclude maladaptive, i.e.: use only primary FSGS and the controls.
seurat_obj <- SetIdent(seurat_obj, value = "Diagnosis")
DimPlot(seurat_obj, label = T)

seurat_obj <- subset(seurat_obj, ident = c("Primary FSGS",
                                           "PLA2R+ MN",
                                           "Healthy control"))

seurat_obj$Diagnosis <- factor(seurat_obj$Diagnosis, levels = c("Primary FSGS",
                                                                "PLA2R+ MN",
                                                                "Healthy control"))

seurat_obj$orig.ident <- factor(seurat_obj$orig.ident, levels=c("snrFSGS001", 
                                                                "snrFSGS006",
                                                                "snrFSGS007",
                                                                "snrFSGS009",
                                                                "snrFSGS011",
                                                                "snrFSGS012",
                                                                "snrFSGS015",
                                                                "snrFSGS016",
                                                                "snrFSGS018",
                                                                "snrNEPH029",
                                                                "snrNEPH031",
                                                                "snrNEPH033",
                                                                "snrNEPH027",
                                                                "snrNEPH028",
                                                                "snrNEPH030",
                                                                "snrNEPH032"))
levels(seurat_obj$orig.ident)

table(seurat_obj$Diagnosis) #ok
table(seurat_obj$orig.ident,seurat_obj$Diagnosis) #ok

# create new metadata-column: FSGS vs control
seurat_obj <- SetIdent(seurat_obj, value = "Diagnosis")
seurat_obj <- RenameIdents(seurat_obj, "PLA2R+ MN" = "Control")
seurat_obj <- RenameIdents(seurat_obj, "Healthy control" = "Control")

DimPlot(seurat_obj, label = T)

seurat_obj$FSGS_vs_control <- Idents(seurat_obj)
seurat_obj$FSGS_vs_control <- factor(seurat_obj$FSGS_vs_control, levels = c("Primary FSGS", "Control"))
seurat_obj <- SetIdent(seurat_obj, value = "FSGS_vs_control")

#Create pseudo-bulk samples
y <- Seurat2PB(seurat_obj, sample="orig.ident", cluster="FSGS_vs_control")

dim(y) #36780    16
head(y$samples, n=10L)
View(y$samples)

write.table(y$samples, file = paste0(projectFolder, ChosenCellGroup, "_EdgeR_pseudobulk_librarySizes.csv"),
            row.names = T, col.names = T, sep = "\t", quote = F)

# We first examine the library sizes of all the pseudo-bulk samples and filter out those below the threshold of 50,000.
summary(y$samples$lib.size)

# Using aggregate to summarize lib.size by cluster
summary_by_cluster <- aggregate(lib.size ~ cluster, data = y$samples, FUN = summary)
print(summary_by_cluster)

# cluster lib.size.Min. lib.size.1st Qu. lib.size.Median lib.size.Mean lib.size.3rd Qu. lib.size.Max.
# 1      Control       72952.0         258983.5        699377.0      722433.3         974269.5     1818198.0
# 2 Primary FSGS       74173.0         291676.0        361414.0      806851.3         490445.0     4107317.0

keep.samples <- y$samples$lib.size > 5e4 #5e4
table(keep.samples)
# TRUE 
# 16

y <- y[, keep.samples]

# We then filter out lowly expressed genes.
# Default values: min.count=10, min.total.count=20
keep.genes <- filterByExpr(y, group=y$samples$cluster, min.count=10, min.total.count=20)

table(keep.genes)
# FALSE  TRUE 
# 31190  5590 

y <- y[keep.genes, , keep=FALSE]

#TMM normalization is performed to estimate effective library sizes.
y <- normLibSizes(y)

head(y$samples, n=10L)

summary(y$samples$norm.factors)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.8554  0.9506  0.9722  1.0097  1.0183  1.5802

cluster <- as.factor(y$samples$cluster)
levels(cluster)
cluster <- factor(cluster, levels = c('Primary FSGS', 'Control'))

pdf(paste0(projectFolder, ChosenCellGroup, "_EdgeR_MDS.pdf"))
plotMDS(y, pch=16, col=c(2:8)[cluster], main="MDS")
legend("bottomright", legend=paste0("cluster",levels(cluster)),pch=16, col=2:8, cex=0.8)
dev.off()

#To perform differential expression analysis between cell clusters, we create a design matrix using both cluster and donor information.
design <- model.matrix(object = ~cluster)

colnames(design)[1] <- "Int"
design
dim(design) #16 2

#The NB dispersion can be estimated using the estimateDisp function and visualized with plotBCV.
y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion

pdf(paste0(projectFolder, ChosenCellGroup, "_EdgeR_plotBCV_dispersion.pdf"))
plotBCV(y)
dev.off()

#  No residual df: setting dispersion to NA ?
fit <- glmQLFit(y, design, robust=TRUE)

pdf(paste0(projectFolder, ChosenCellGroup, "_EdgeR_plotQLDispFit.pdf"))
plotQLDisp(fit)
dev.off()

# To confirm the identities of cell clusters, we perform differential expression analysis to identify marker genes of each cluster. 
# In particular, we compare each cluster with all the other clusters. 
# We construct a contrast matrix as follows so that each column of the contrast matrix represents a testing contrast for one cell cluster.
ncls <- nlevels(cluster)
contr <- rbind( matrix(1/(1-ncls), ncls, ncls), matrix(0, ncol(design)-ncls, ncls) )
diag(contr) <- 1 # diagonal gets 1
contr[1,] <- 0 # first row gets 0
rownames(contr) <- colnames(design)
colnames(contr) <- paste0("cluster", levels(cluster))
contr

# clusterPrimary FSGS clusterControl
# Int                              0              0
# clusterControl                  -1              1

# We then perform quasi-likelihood F-test for each testing contrast. (glm=general linear models)
# The results are stored as a list of DGELRT objects, one for each comparison.
qlf <- list()
for(i in 1:ncls){
  qlf[[i]] <- glmQLFTest(fit, contrast=contr[,i])
  qlf[[i]]$comparison <- paste0("cluster", levels(cluster)[i], "_vs_others") 
}

#The top most significant DE genes of cluster 0 vs other clusters can be examined with topTags.
topTags(qlf[[1]], n=10L) #clusterPrimary FSGS_vs_others
topTags(qlf[[2]], n=10L) #clusterControl_vs_others 

#The numbers of DE genes under each comparison are shown below.
dt <- lapply(lapply(qlf, decideTestsDGE), summary)
dt.all <- do.call("cbind", dt)
dt.all

# clusterPrimary FSGS_vs_others clusterControl_vs_others
# Down                               0                        0
# NotSig                          5590                     5590
# Up                                 0                        0

#A heat map is produced to visualize the top marker genes across all the pseudo-bulk samples.
lcpm <- edgeR::cpm(y, log=TRUE)
annot <- data.frame(cluster=paste0("cluster ", cluster))
rownames(annot) <- colnames(y)
ann_colors <- list(cluster=2:5)
names(ann_colors$cluster) <- paste0("cluster ", levels(cluster))

# Save the objects in a single file
save(lcpm, annot, ann_colors, file = paste0(projectFolder, ChosenCellGroup, "_EdgeR_primary_vs_controls_DEG_analysis_conditions_lcpm.RData"))

# Follow-up visualizations with heatmaps & volcanoplots are analogous to upstream R-code.
# Follow-up fgsea analysis is analogous to upstream R-code.


##### 6. PART 6 - DGE ANALYSIS OF PRIMARY FSGS VS. MALADAPTIVE FSGS - EdgeR (pseudobulk)  #####

# Switch to R-version with installation of EdgeR

# Step 1: run R 4.3
library(Seurat)
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(edgeR)
library(EnhancedVolcano)

projectFolder <- "path/seurat_POD_HQ2/Diagnosis_analysis_primary_vs_maladaptive/EdgeR/"

###### 6.1: Run general EdgeR DEG analysis ###### 

# Read in the dataset
seurat_POD_HQ2 <- readRDS(file = "path/seurat_POD_HQ2/seurat_POD_HQ2_harmony_annot.Rds")
seurat_obj <- seurat_POD_HQ2 #insert object
ChosenCellGroup <- "seurat_POD_HQ2"

#check conditions
table(seurat_obj$Diagnosis)
table(seurat_obj$orig.ident,seurat_obj$Diagnosis)

# Exclude controls, i.e.: use only primary FSGS and maladaptive FSGS
seurat_obj <- SetIdent(seurat_obj, value = "Diagnosis")
DimPlot(seurat_obj, label = T)

seurat_obj <- subset(seurat_obj, ident = c("Primary FSGS",
                                           "Maladaptive FSGS"))

seurat_obj$Diagnosis <- factor(seurat_obj$Diagnosis, levels = c("Primary FSGS",
                                                                "Maladaptive FSGS"))

seurat_obj$orig.ident <- factor(seurat_obj$orig.ident, levels=c("snrFSGS001", 
                                                                "snrFSGS006",
                                                                "snrFSGS007",
                                                                "snrFSGS009",
                                                                "snrFSGS011",
                                                                "snrFSGS012",
                                                                "snrFSGS015",
                                                                "snrFSGS016",
                                                                "snrFSGS018",
                                                                "snrFSGS002",
                                                                "snrFSGS003",
                                                                "snrFSGS004",
                                                                "snrFSGS005",
                                                                "snrFSGS008",
                                                                "snrFSGS010",
                                                                "snrFSGS013",
                                                                "snrFSGS014",
                                                                "snrFSGS017"))
levels(seurat_obj$orig.ident)

table(seurat_obj$Diagnosis) #ok
table(seurat_obj$orig.ident,seurat_obj$Diagnosis) #ok

seurat_obj <- SetIdent(seurat_obj, value = "Diagnosis")

#Create pseudo-bulk samples
y <- Seurat2PB(seurat_obj, sample="orig.ident", cluster="Diagnosis")

dim(y) #36780    18
head(y$samples, n=10L)
View(y$samples)

write.table(y$samples, file = paste0(projectFolder, ChosenCellGroup, "_EdgeR_pseudobulk_librarySizes.csv"),
            row.names = T, col.names = T, sep = "\t", quote = F)

# We first examine the library sizes of all the pseudo-bulk samples and filter out those below the threshold of 50,000.
summary(y$samples$lib.size)

# Using aggregate to summarize lib.size by cluster
summary_by_cluster <- aggregate(lib.size ~ cluster, data = y$samples, FUN = summary)
print(summary_by_cluster)

# cluster lib.size.Min. lib.size.1st Qu. lib.size.Median lib.size.Mean lib.size.3rd Qu. lib.size.Max.
# 1 Maladaptive FSGS       19127.0          79504.0        127826.0      333812.4         277133.0     1290836.0
# 2     Primary FSGS       74173.0         291676.0        361414.0      806851.3         490445.0     4107317.0

keep.samples <- y$samples$lib.size > 5e4 #5e4
table(keep.samples)
# FALSE  TRUE 
# 2    16

y <- y[, keep.samples]

# We then filter out lowly expressed genes.
# Default values: min.count=10, min.total.count=20
keep.genes <- filterByExpr(y, group=y$samples$cluster, min.count=10, min.total.count=20)

table(keep.genes)
# FALSE  TRUE 
# 31883  4897

y <- y[keep.genes, , keep=FALSE]

#TMM normalization is performed to estimate effective library sizes.
y <- normLibSizes(y)

head(y$samples, n=10L)

summary(y$samples$norm.factors)

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.8611  0.9366  0.9525  1.0128  0.9854  1.5391

cluster <- as.factor(y$samples$cluster)
levels(cluster)
cluster <- factor(cluster, levels = c('Primary FSGS', 'Maladaptive FSGS'))

pdf(paste0(projectFolder, ChosenCellGroup, "_EdgeR_MDS.pdf"))
plotMDS(y, pch=16, col=c(2:8)[cluster], main="MDS")
legend("bottomright", legend=paste0("cluster",levels(cluster)),pch=16, col=2:8, cex=0.8)
dev.off()

#To perform differential expression analysis between cell clusters, we create a design matrix using both cluster and donor information.
design <- model.matrix(object = ~cluster)

colnames(design)[1] <- "Int"
design
dim(design) #16 2

#The NB dispersion can be estimated using the estimateDisp function and visualized with plotBCV.
y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion

pdf(paste0(projectFolder, ChosenCellGroup, "_EdgeR_plotBCV_dispersion.pdf"))
plotBCV(y)
dev.off()

#  No residual df: setting dispersion to NA ?
fit <- glmQLFit(y, design, robust=TRUE)

pdf(paste0(projectFolder, ChosenCellGroup, "_EdgeR_plotQLDispFit.pdf"))
plotQLDisp(fit)
dev.off()

# To confirm the identities of cell clusters, we perform differential expression analysis to identify marker genes of each cluster. 
# In particular, we compare each cluster with all the other clusters. 
# We construct a contrast matrix as follows so that each column of the contrast matrix represents a testing contrast for one cell cluster.
ncls <- nlevels(cluster)
contr <- rbind( matrix(1/(1-ncls), ncls, ncls), matrix(0, ncol(design)-ncls, ncls) )
diag(contr) <- 1 # diagonal gets 1
contr[1,] <- 0 # first row gets 0
rownames(contr) <- colnames(design)
colnames(contr) <- paste0("cluster", levels(cluster))
contr

# clusterPrimary FSGS clusterMaladaptive FSGS
# Int                                       0                       0
# clusterMaladaptive FSGS                  -1                       1

# We then perform quasi-likelihood F-test for each testing contrast. (glm=general linear models)
# The results are stored as a list of DGELRT objects, one for each comparison.
qlf <- list()
for(i in 1:ncls){
  qlf[[i]] <- glmQLFTest(fit, contrast=contr[,i])
  qlf[[i]]$comparison <- paste0("cluster", levels(cluster)[i], "_vs_others") 
}

#The top most significant DE genes of cluster 0 vs other clusters can be examined with topTags.
topTags(qlf[[1]], n=10L) #clusterPrimary FSGS_vs_others
topTags(qlf[[2]], n=10L) #clusterMaladaptive FSGS_vs_others 

#The numbers of DE genes under each comparison are shown below.
dt <- lapply(lapply(qlf, decideTestsDGE), summary)
dt.all <- do.call("cbind", dt)
dt.all

# clusterPrimary FSGS_vs_others clusterMaladaptive FSGS_vs_others
# Down                               0                                 0
# NotSig                          4897                              4897
# Up                                 0                                 0

#A heat map is produced to visualize the top marker genes across all the pseudo-bulk samples.
lcpm <- edgeR::cpm(y, log=TRUE)
annot <- data.frame(cluster=paste0("cluster ", cluster))
rownames(annot) <- colnames(y)
ann_colors <- list(cluster=2:5)
names(ann_colors$cluster) <- paste0("cluster ", levels(cluster))

# Save the objects in a single file
save(lcpm, annot, ann_colors, file = paste0(projectFolder, ChosenCellGroup, "_EdgeR_primary_vs_maladaptive_DEG_analysis_conditions_lcpm.RData"))

# Follow-up visualizations with heatmaps & volcanoplots are analogous to upstream R-code.
# Follow-up fgsea analysis is analogous to upstream R-code.