#### MAIN CLUSTERING, ANNOTATION and ABUNDANCY ANALYSIS ####

##### Script info #####

# This script starts from loading in Cellbender-corrected indivual samples, creating a merged object, and performing QC.

# 1. PART 1: We run Seurat pipeline for the first time with all lowQ cells still included (total of 151,909 nuclei). 
# We perform first QC checks with first removal of lowQ cells (8,644 lowQ nuclei removed).

# 2. PART 2: We re-run Seurat pipeline for the second time after first removal of lowQ cells in part 1. 
# We perform additional QC checks. We also use SingleR to import the cell annotation of a publicly available snRNA-seq kidney dataset from KPMP, which aids in annotation of our Seurat object. 
# We save all cell subclusters as separate .rds-objects. We perform QC checks on every individual subcluster level and identify additional lowQ cells (22,514 low-quality nuclei) on subcluster level.

# 3. PART 3: We re-map all "clean" (high-quality) subclusters back to the main annotation level. 
# All the lowQ cells that were identified on subcluster level in part 2 are now also visualized on main annotation level and subsequently removed. 
# This concludes QC of the main Seurat object, yielding 120,751 high-quality nuclei.

# 4. PART 4: We re-run Seurat pipeline for a final 3rd time to recreate a clean and reclustered UMAP; no additional cells are removed in this step. 
# The final highQ Seurat-object containing all celltypes is called "seurat_merged_HQ2_after_ReRunUMAP.Rds" (120,751 high-quality nuclei)
# The final main annotation can be found in meta.data "seurat_merged_HQ2$main_annot_HQ2"
# The final subcluster annotation can be found in meta.data "seurat_merged_HQ2$subcluster_annot_V2"

# 5. PART 5: We calculate cell cluster frequencies and abundancy analysis on the final main Seurat object ("seurat_merged_HQ2_after_ReRunUMAP.Rds"), using annotation in meta.data "seurat_merged_HQ2$main_annot_HQ2".


##### 1. PART 1 - Running pipeline for first time - all lowQ still included #####

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
library(SingleR)
library(SingleCellExperiment)
library(SeuratDisk)

###### 1.1: Reading in sample info ######

###### 1.2: Summary files, lower bound filtering, running DoubletFinder ######

## Create summary files of all sequencing metrics
for(i in sampleIDs){
	
	print(i)
	
  ## Finding sequencing saturation and number of reads info in web summary (from raw data)
  summaryHtmlFile <- grep(i, summaryHTMLs, value = T)
  summaryHtml <- readLines(paste0(getwd(), "/path/", summaryHtmlFile))
  if(length(grep("<td>Sequencing Saturation</td>", summaryHtml) ==1) > 0){
    saturation <- strsplit(strsplit(summaryHtml[grep("<td>Sequencing Saturation</td>", summaryHtml) + 1],
                                    "<td>")[[1]][2], "</td>")[[1]][1]
    reads <- strsplit(strsplit(summaryHtml[grep("<td>Number of Reads</td>", summaryHtml) + 1],
                               "<td>")[[1]][2], "</td>")[[1]][1]
  }else{
    saturation <- strsplit(strsplit(summaryHtml[grep("\"Sequencing Saturation\", \"", summaryHtml)],
                                    "\"Sequencing Saturation\", \"")[[1]][2], "\"")[[1]][1]
    reads <- strsplit(strsplit(summaryHtml[grep("\"Number of Reads\", \"", summaryHtml)],
                               "\"Number of Reads\", \"")[[1]][2], "\"")[[1]][1]
  }
  
  ## Creating 10X dataset
  tenX <- Read10X_h5(filename = grep(i, hdf5Files_full, value = T))
  seurat <- CreateSeuratObject(counts = tenX, min.cells = 1, min.features = 1, project = i)
  seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
  
  # Filtering 
  seurat <- subset(seurat, subset = nFeature_RNA > 1000 & nCount_RNA > 1000 & percent.mt < 5) # prior to DoubletFinder, we do not yet filter the upper bounds. 

  # Pre-process Seurat object for doubletFinder
  seurat <- NormalizeData(seurat)
  seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
  seurat <- ScaleData(seurat)
  seurat <- RunPCA(seurat)
  seurat <- RunUMAP(seurat, dims = 1:10)
  
  # Run doubletFinder, pK Identification (no ground-truth)
  sweep.res.list <- paramSweep_v3(seurat, PCs = 1:10, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  seurat_pK <- find.pK(sweep.stats)
  
  png(paste0(projectFolder, "BCmetric_doubletFinder_", i, ".png"), height=550, width=350)
  barplot(seurat_pK$BCmetric, names.arg = seurat_pK$pK, las=2)
  dev.off()
  
  pK_selected <- seurat_pK$pK[seurat_pK$BCmetric == max(seurat_pK$BCmetric)]
  pK_selected <- as.numeric(levels(pK_selected)[pK_selected])
  expDoublets <- ifelse(dim(seurat)[2] < 2000, 0.01, ifelse(dim(seurat)[2] > 10000, 0.08, 0.05))
  nExp <- round(ncol(seurat) * expDoublets)
  seurat <- doubletFinder_v3(seurat, pN = 0.25, pK = pK_selected, nExp = nExp, PCs = 1:10)
  
  seurat_objs[[length(seurat_objs) + 1]] <- seurat
  names(seurat_objs)[length(seurat_objs)] <- i
  
  sampleSummary <- data.frame(sampleID = i,
                              numberOfReads = reads,
                              saturation = saturation,
							  numberOfCells = dim(seurat)[2],
							  expDoublets = expDoublets,
							  nExp = nExp)
  summaryDF <- rbind(summaryDF, sampleSummary)
  
}

write.table(summaryDF, file = paste0(projectFolder, "Metadata_summaryTable1.csv"),
            sep = "\t", quote = F, row.names = F, col.names = T)
			
# We tested different upper bounds for nFeatures cut-off and finally chose 4000 as upper bound filtering treshold.
upperBound <- 4000

if (upperBound == 4000){

	for(i in sampleIDs){
	
		print(i)
	
		seurat_objs[[i]] <- subset(seurat_objs[[i]], subset = nFeature_RNA < upperBound) # now filter upper bound, after running DoubletFinder

	}

}

## DimPlots showing doublets for each sample
pdf(paste0(projectFolder, "dimPlot_doubletFinder_allClassifications.pdf"), height=10, width=10)
for(i in 1:length(seurat_objs)){
  doublets <- colnames(seurat_objs[[i]]@meta.data)[grepl("DF.classification", colnames(seurat_objs[[i]]@meta.data))]
  print(DimPlot(seurat_objs[[i]], group.by = doublets) + ggtitle(doublets))
}
dev.off()

## ViolinPlots showing differences between doublets and singlets for each sample
pdf(paste0(projectFolder, "ViolinPlot_doubletFinder_allClassifications.pdf"), height=10, width=10)
for(i in 1:length(seurat_objs)){
  doublets <- colnames(seurat_objs[[i]]@meta.data)[grepl("DF.classification", colnames(seurat_objs[[i]]@meta.data))]
  print(VlnPlot(seurat_objs[[i]], features = "nFeature_RNA", group.by = doublets, pt.size = 0.1) + ggtitle(doublets))
}
dev.off()

## Merging Seurat objects (when all in seurat_objs object)
seurat_merged <- merge(seurat_objs[[1]], y = c(seurat_objs[2:length(seurat_objs)]),
                       add.cell.ids = names(seurat_objs))


###### 1.3: Adding metadata ###### 
# Factorizing and re-ordering orig.ident according to diagnosis, factorizing and re-ordering 'Diagnosis'

## Re-factor orig.ident according to diagnoses
levels(seurat_merged$orig.ident)
seurat_merged$orig.ident <- factor(seurat_merged$orig.ident, levels=c("snrFSGS001", 
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
                                                                                                "snrFSGS017",
                                                                                                "snrNEPH029",
                                                                                                "snrNEPH031",
                                                                                                "snrNEPH033",
                                                                                                "snrNEPH027",
                                                                                                "snrNEPH028",
                                                                                                "snrNEPH030",
                                                                                                "snrNEPH032"))
levels(seurat_merged$orig.ident)

## Re-factor diagnoses
levels(seurat_merged$Diagnosis)
seurat_merged$Diagnosis <- factor(seurat_merged$Diagnosis, levels= c("Primary FSGS", 
                                                                                               "Maladaptive FSGS", 
                                                                                               "PLA2R+ MN", 
                                                                                               "Healthy control"))
levels(seurat_merged$Diagnosis)

###### 1.4: QC-plots and tables (per sample, per diagnosis) ###### 
seurat_merged <- SetIdent(seurat_merged, value = "orig.ident")

# add MALAT1-expression as an additional QC parameter
seurat_merged[["percent.MALAT1"]] <- PercentageFeatureSet(seurat_merged, pattern = "MALAT1")

# QC Per sample - VlnPlot

pdf(paste0(projectFolder, "upperBoundnFeatures_", upperBound, "_seurat_merged_QC_per_sample_VlnPlot.pdf"), width = 14, height = 7)
VlnPlot(seurat_merged, features = c("nCount_RNA"),
        pt.size = 0, 
        cols = rep("#1F78B4", times = 25)) &
  geom_boxplot(width=0.1, fill="white", outlier.size = 0.5) & 
  NoLegend() &
  stat_summary(fun.y = median, geom='point', size = 2, colour = "black")

VlnPlot(seurat_merged, features = c("nFeature_RNA"), 
        pt.size = 0, 
        cols = rep("#1F78B4", times = 25)) &
  geom_boxplot(width=0.1, fill="white", outlier.size = 0.5) & 
  NoLegend() &
  stat_summary(fun.y = median, geom='point', size = 2, colour = "black")

VlnPlot(seurat_merged, features = c("percent.mt"), 
        pt.size = 0, 
        cols = rep("#1F78B4", times = 25)) &
  geom_boxplot(width=0.1, fill="white", outlier.size = 0.5) & 
  NoLegend() &
  stat_summary(fun.y = median, geom='point', size = 2, colour = "black")

VlnPlot(seurat_merged, features = c("percent.MALAT1"), 
      pt.size = 0, 
      cols = rep("#1F78B4", times = 25)) &
geom_boxplot(width=0.1, fill="white", outlier.size = 0.5) & 
NoLegend() &
stat_summary(fun.y = median, geom='point', size = 2, colour = "black")

dev.off()


## QC per sample, table
df <- data.frame(matrix(ncol = 5, nrow = length(seurat_merged$orig.ident)))
x <- c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "Diagnosis")
colnames(df) <- x

df$orig.ident <- seurat_merged$orig.ident
df$nCount_RNA <- seurat_merged$nCount_RNA
df$nFeature_RNA <- seurat_merged$nFeature_RNA
df$percent.mt <- seurat_merged$percent.mt
df$percent.MALAT1 <- seurat_merged$percent.MALAT1
df$Diagnosis <- seurat_merged$Diagnosis

seuratDF <- df %>%
	group_by(orig.ident) %>%
	summarize(  Diagnosis = unique(Diagnosis),
				Median_nCounts = round(median(nCount_RNA),0), 
				Mean_nCounts = round(mean(nCount_RNA),0),
				Median_nFeatures = round(median(nFeature_RNA),0),
				Mean_nFeatures = round(mean(nFeature_RNA),0),
				Median_percent.mt = round(median(percent.mt),2),
				Mean_percent.mt = round(mean(percent.mt),2),
				Median_percent.MALAT1 = round(median(percent.MALAT1),2),
				Mean_percent.MALAT1 = round(mean(percent.MALAT1),2),
				nCells = n()
			)

colnames(seuratDF) <- c("orig.ident","Diagnosis","Median\nnCounts", "Mean\nnCounts", "Median\nnFeatures", "Mean\nnFeatures", "Median\npercent.mt", "Mean\npercent.mt", "Median\npercent.MALAT1", "Mean\npercent.MALAT1", "Nuclei")
seuratDF

pdf(paste0(projectFolder, "upperBoundnFeatures_", upperBound, "_seurat_merged_QC_perSample_table_numberOfCellsCountsFeaturesMito.pdf"), height=8, width=13)
p <- tableGrob(seuratDF)
grid.arrange(p)
dev.off()

colnames(seuratDF) <- c("orig.ident","Diagnosis","Median nCounts", "Mean nCounts", "Median nFeatures", "Mean nFeatures", "Median percent.mt", "Mean percent.mt", "Median percent.MALAT1", "Mean percent.MALAT1", "Nuclei")
seuratDF

write.table(seuratDF, file = paste0(projectFolder, "upperBoundnFeatures_", upperBound, "_seurat_merged_QC_perSample_table_numberOfCellsCountsFeaturesMito.csv"),
            quote = F, sep = ";", col.names = T, row.names = T)

# QC per Diagnosis - VlnPlot
levels(seurat_merged$Diagnosis)
seurat_merged <- SetIdent(seurat_merged, value = "Diagnosis")

pdf(paste0(projectFolder, "upperBoundnFeatures_", upperBound, "_seurat_merged_QC_perDiagnosis_VlnPlot.pdf"), width = 14, height = 5)
VlnPlot(seurat_merged, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.MALAT1"),
        pt.size = 0, ncol=5,
        cols = rep("#1F78B4", times = 25)) &
  geom_boxplot(width=0.1, fill="white", outlier.size = 0.5) & 
  NoLegend() &
  stat_summary(fun.y = median, geom='point', size = 2, colour = "black")

dev.off()

## QC per sample - table
seuratDF <- data.frame(Median_nCounts = round(tapply(seurat_merged$nCount_RNA, seurat_merged$Diagnosis, median)),
                       Mean_nCounts = round(tapply(seurat_merged$nCount_RNA, seurat_merged$Diagnosis, mean)),
                       Median_nFeatures = round(tapply(seurat_merged$nFeature_RNA, seurat_merged$Diagnosis, median)),
                       Mean_nFeatures = round(tapply(seurat_merged$nFeature_RNA, seurat_merged$Diagnosis, mean)),
                       Median_percent.mt = round(tapply(seurat_merged$percent.mt, seurat_merged$Diagnosis, median), 2),
                       Mean_percent.mt = round(tapply(seurat_merged$percent.mt, seurat_merged$Diagnosis, mean), 2),
                       Median_percent.MALAT1 = round(tapply(seurat_merged$percent.MALAT1, seurat_merged$Diagnosis, median), 2),
                       Mean_percent.MALAT1 = round(tapply(seurat_merged$percent.MALAT1, seurat_merged$Diagnosis, mean), 2),
                       nCells = table(seurat_merged$Diagnosis))
rownames(seuratDF) <- seuratDF$nCells.Var1
seuratDF$nCells.Var1 <- NULL
colnames(seuratDF) <- c("Median\nnCounts", "Mean\nnCounts", "Median\nnFeatures", "Mean\nnFeatures", "Median\npercent.mt", "Mean\npercent.mt", "Median\npercent.MALAT1", "Mean\npercent.MALAT1", "Nuclei")

pdf(paste0(projectFolder, "upperBoundnFeatures_", upperBound, "_seurat_merged_QC_perDiagnosis_table_numberOfCellsCountsFeaturesMito.pdf"), height=4, width=12)
p <- tableGrob(seuratDF)
grid.arrange(p)
dev.off()

colnames(seuratDF) <- c("Median nCounts", "Mean nCounts", "Median nFeatures", "Mean nFeatures", "Median percent.mt", "Mean percent.mt", "Median percent.MALAT1", "Mean percent.MALAT1", "Nuclei")
seuratDF

write.table(seuratDF, file = paste0(projectFolder, "upperBoundnFeatures_", upperBound, "_seurat_merged_QC_perDiagnosis_table_numberOfCellsCountsFeaturesMito.csv"),
            quote = F, sep = ";", col.names = T, row.names = T)


###### 1.5: Running Seurat pipeline ###### 

## Normalizing the data
seurat_merged <- NormalizeData(seurat_merged, normalization.method = "LogNormalize", scale.factor = 10000)

## Highly variable features
seurat_merged <- FindVariableFeatures(seurat_merged, selection.method = "vst", nfeatures = 2000)

## Scaling the data + regressing for mt.RNA, nCount, nFeature
seurat_merged <- ScaleData(seurat_merged, vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt")) # we only use variable features

## PCA - linear dimensional reduction
seurat_merged <- RunPCA(seurat_merged, features = VariableFeatures(object = seurat_merged), npcs = 50)

pdf(paste0(projectFolder, "seurat_merged_noHarmony_elbowPlot.pdf"), width = 10, height = 10)
ElbowPlot(seurat_merged, n = 50)
dev.off()

## Clustering
numPCs <- 25
seurat_merged <- FindNeighbors(seurat_merged, dims = 1:numPCs)
res = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1, 1.2, 1.5, 1.7, 2)
seurat_merged <- FindClusters(seurat_merged, resolution = res)
seurat_merged <- RunUMAP(seurat_merged, dims = 1:numPCs)

pdf(paste0(projectFolder, "seurat_merged_noHarmony_dimPlot_umap_multiRes_nPC",numPCs,".pdf"), width = 10, height = 8)
for (i in 1:length(res)){
  resolution <- paste0("RNA_snn_res." , res[i])
  print(DimPlot(seurat_merged, reduction = "umap", group.by = resolution, label = T) + ggtitle(paste0("Resolution: ", res[i])))
}
dev.off()

## Dimplot - (to check if Harmony needed)
pdf(paste0(projectFolder, "seurat_merged_noHarmony_dimplot_sample_diagnosis_processing_nPC",numPCs,".pdf"), width = 10, height = 8)
DimPlot(seurat_merged, group.by = "orig.ident", label = T, raster=FALSE)
DimPlot(seurat_merged, group.by = "Diagnosis", label = T, raster=FALSE)
DimPlot(seurat_merged, group.by = "Date_processing", label = T, raster=FALSE)
dev.off()

## Highlighting samples on the plot (i.e. plotting cells from individual samples) ##
seurat_merged <- SetIdent(seurat_merged, value = "orig.ident")

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
pdf(paste0(projectFolder, "seurat_merged_noHarmony_samplesHighlighted_Dimplot_nPC",numPCs,".pdf"), width = 10, height = 10)
lapply(seq_along(Lorigin),function(i) 
  DimPlot(seurat_merged, cells.highlight = WhichCells(seurat_merged, idents = Lorigin[[i]])) & ggtitle(names(Lorigin)[[i]]))
dev.off()

## Checking for stress-score, hypoxia-score and cell-cycling score ##
## Stress genes scoring
stress.genes <- read.table("stress.genes.csv", header = T)
stressGenes <- stress.genes$gene[stress.genes$gene %in% rownames(seurat_merged)]

seurat_merged <- AddModuleScore(object = seurat_merged,
                                     features = list(stressGenes),
                                     ctrl = length(stressGenes),
                                     name = 'Stress.Score')

## Hypoxia scoring
genes.hypoxia<- c("VEGFA", "SLC2A1", "PGAM1", "ENO1","LDHA", "TPI1", "P4HA1", "MRPS17",
                  "CDKN3", "ADM", "NDRG1", "TUBB6","ALDOA", "MIF", "ACOT7")

hypoxiaGenes <- genes.hypoxia[genes.hypoxia %in% rownames(seurat_merged)]

seurat_merged <- AddModuleScore(object = seurat_merged,
                                     features = list(hypoxiaGenes),
                                     ctrl = length(hypoxiaGenes),
                                     name = 'Hypoxia.Score')

## Cell cycle scoring
cc.genes <- readLines(con = "regev_lab_cell_cycle_genes.txt")
sGenes <- cc.genes[1:43]
g2mGenes <- cc.genes[44:97]

sGenes <- sGenes[sGenes %in% rownames(seurat_merged)]
g2mGenes <- g2mGenes[g2mGenes %in% rownames(seurat_merged)]

seurat_merged <- CellCycleScoring(object = seurat_merged,
                                       s.features = sGenes,
                                       g2m.features = g2mGenes,
                                       set.ident = TRUE)

pdf(paste0(projectFolder, "seurat_merged_dimplot_stress_hypoxia_cycling_features_counts_percentMT_scores_nPC",numPCs,".pdf"), width = 15, height = 15)
FeaturePlot(object = seurat_merged,
            features = c("Stress.Score1", "Hypoxia.Score1", "S.Score", "G2M.Score",
                         "nCount_RNA", "nFeature_RNA", "percent.mt", "percent.MALAT1"),
            cols = c("grey", "blue"), pt.size = 1,
            reduction = "umap")
dev.off()

# When using 25 PCs, Harmony seems to be necessary as 2 samples seem to cluster separately

###### 1.6: Harmony integration ###### 

## Using Harmony
seurat_merged <- RunHarmony(seurat_merged, "orig.ident")

## Clustering
numPCs <- 25
seurat_merged <- FindNeighbors(seurat_merged, reduction = "harmony", dims = 1:numPCs)
res = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1, 1.2, 1.5, 1.7, 2)
seurat_merged <- FindClusters(seurat_merged, resolution = res)
seurat_merged <- RunUMAP(seurat_merged, reduction = "harmony", dims = 1:numPCs)

pdf(paste0(projectFolder, "seurat_merged_Harmony_dimPlot_umap_multiRes_nPC",numPCs,".pdf"), width = 10, height = 8)
for (i in 1:length(res)){
  resolution <- paste0("RNA_snn_res." , res[i])
  print(DimPlot(seurat_merged, reduction = "umap", group.by = resolution, label = T, pt.size = 0.7, raster.dpi = c(1024,1024)) + ggtitle(paste0("Resolution: ", res[i])))
}
dev.off()

## Dimplot - to check after Harmony
pdf(paste0(projectFolder, "seurat_merged_Harmony_dimplot_sample_diagnosis_processing_nPC",numPCs,".pdf"), width = 10, height = 8)
DimPlot(seurat_merged, group.by = "orig.ident", label = T, raster=FALSE)
DimPlot(seurat_merged, group.by = "Diagnosis", label = T, raster=FALSE)
DimPlot(seurat_merged, group.by = "Date_processing", label = T, raster=FALSE)
dev.off()

## Feature plot of features, counts and %mtRNA (after Harmony)
pdf(paste0(projectFolder, "seurat_merged_harmony_featureplot_stress_hypoxia_cycling_features_counts_percentMT_nPC",numPCs,".pdf"), width = 15, height = 15)
FeaturePlot(object = seurat_merged,
            features = c("Stress.Score1", "Hypoxia.Score1", "S.Score", "G2M.Score",
                         "nCount_RNA", "nFeature_RNA", "percent.mt", "percent.MALAT1"),
            cols = c("grey", "blue"), pt.size = 1,
            reduction = "umap")
dev.off()


###### 1.7: Re-mapping annotation of previous final MainAnnot object ###### 
# We imported the cell annotation from a preliminary analysis of a subset of cells from the current Seurat-object, which can be found in meta.data under "subcluster_annot".

table(seurat_merged$subcluster_annot)

## Plotting
seurat_merged <- SetIdent(seurat_merged, value = "subcluster_annot")

pdf(paste0(projectFolder, "seurat_merged_DimPlot_HQ_mapped.pdf"), width = 11, height = 10)
DimPlot(seurat_merged, group.by = "subcluster_annot", label = T)
dev.off()

###### 1.8: Choosing resolution + QC per cluster ###### 
ChosenRes <- "res12"

seurat_merged$RNA_snn_res.1.2 <- factor(seurat_merged$RNA_snn_res.1.2, levels = c(0:(length(levels(seurat_merged$RNA_snn_res.1.2)) -1)))
seurat_merged <- SetIdent(seurat_merged, value = "RNA_snn_res.1.2")
DimPlot(seurat_merged, label=T)

## DimPlot
pdf(paste0(projectFolder, "seurat_merged_harmony_", ChosenRes, "_dimPlot_umap.pdf"), width = 8, height = 7)
DimPlot(seurat_merged, reduction = "umap", label = TRUE, pt.size = 0.7, raster.dpi = c(1024,1024)) + ggtitle("seurat_merged_harmony_res1.2_dimPlot_umap")
dev.off()

## Plotting the 'unknown' and 'lowQ' from "subcluster_annot"
seurat_merged <- SetIdent(seurat_merged, value = "subcluster_annot")

Lcluster <- list(
  cluster_B_cell="B_cell",
  cluster_CNT="CNT",
  cluster_DCT="DCT",
  cluster_DTL="DTL",
  cluster_EC_AEA_DVR="EC_AEA_DVR",
  cluster_EC_GC="EC_GC",
  cluster_EC_others="EC_others",
  cluster_EC_PTC_AVR="EC_PTC_AVR",
  cluster_FIB="FIB",
  cluster_IC_A="IC_A",
  cluster_IC_B="IC_B",
  cluster_lowQ="lowQ",
  cluster_MES="MES",
  cluster_MYEL="MYEL",
  cluster_PEC="PEC",
  cluster_POD="POD",
  cluster_PT="PT",
  cluster_REN="REN",
  cluster_T_NKT="T_NKT",
  cluster_TAL="TAL",
  cluster_unknown="unknown",
  cluster_VSMCP="VSMC/P")

# Highlight the clusters on the dimplot 
pdf(paste0(projectFolder, "seurat_merged_harmony_subcluster_annot_clustersHighlighted_Dimplot.pdf"), width = 10, height = 10)
lapply(seq_along(Lcluster),function(i) 
  DimPlot(seurat_merged, cells.highlight = WhichCells(seurat_merged, idents = Lcluster[[i]])) & ggtitle(names(Lcluster)[[i]]))
dev.off()

###### 1.9: Checking doublets ###### 

## Visualizing doublets on UMAP
metaDF <- S4ToList(seurat_merged)$meta.data # converts the metadata-slot to a list, called metaDF
seuratObjMetaNames <- names(metaDF) # takes the column names of metaDF and assign to new object
doubletsData <- grep("^DF.classification.", seuratObjMetaNames, value = T) # extracts names that start with "DF.classification"
doubletDF <- metaDF[,colnames(metaDF) %in% doubletsData] # selects the columns in metaDF that have column names in doubletsData
mergedDoubletInfo <- coalesce(!!!doubletDF) # collapses multiple columns in doubletDF into a single column

seurat_merged$mergedDoubletInfo <- mergedDoubletInfo

pdf(paste0(projectFolder, "seurat_merged_harmony_dimPlot_", ChosenRes, "_umap_allMergedDoublets.pdf"), width = 7, height = 7)
DimPlot(seurat_merged, group.by = "mergedDoubletInfo", reduction = "umap", label = F)
dev.off()

## Average doublets score per cluster
metaDF <- S4ToList(seurat_merged)$meta.data
seuratObjMetaNames <- names(metaDF)
doubletsData <- grep("^pANN_", seuratObjMetaNames, value = T)
doubletDF <- metaDF[,colnames(metaDF) %in% doubletsData]
mergedDoubletInfo <- coalesce(!!!doubletDF)

metaDF$mergedDoubletInfo2 <- mergedDoubletInfo

doubletAvg <- data.frame(tapply(metaDF$mergedDoubletInfo2, metaDF$RNA_snn_res.1.2, mean)) # applies the tapply() function to group the values in "mergedDoubletInfo2" by the values in "RNA_snn_res.0.2", and calculates the mean of each group. 

doubletAvg$cluster <- rownames(doubletAvg)
doubletAvg$cluster <- factor(doubletAvg$cluster, levels = 0:nrow(doubletAvg))
colnames(doubletAvg) <- c("Doublet", "cluster")
barColour <- ifelse(doubletAvg$Doublet > 0.3, "#c45c3d", "#48a7c2")

pdf(paste0(projectFolder, "seurat_merged_harmony_", ChosenRes, "_avgOfDoubletsScorePerCell.pdf"), width = 3, height = 6)
ggplot(doubletAvg, aes(x = cluster, y = Doublet)) + geom_bar(stat = "identity", fill = barColour) + theme_bw() +
  geom_hline(yintercept = 0.3) + labs(y="Average doublet score", x = "Cluster") + coord_flip()
dev.off()

## Creating new UMAP with cluster IDs only highlighted for clusters that have high doublet score
doubletAvgTable <- tapply(metaDF$mergedDoubletInfo2, metaDF$RNA_snn_res.1.2, mean)

seurat_merged[["highDoublet_clusters"]] <- as.character(seurat_merged$RNA_snn_res.1.2)
seurat_merged$highDoublet_clusters[seurat_merged$RNA_snn_res.1.2 %in% names(doubletAvgTable)[doubletAvg < 0.3]] <- "Low doublets"

pdf(paste0(projectFolder, "seurat_merged_harmony_", ChosenRes, "_DimPlot_umap_highDoublets.pdf"), width = 7, height = 7)
DimPlot(seurat_merged, group.by = "highDoublet_clusters", reduction = "umap", label = T)
dev.off()

###### 1.10: FindAllMarkers + plotting of DEGs ###### 

## Finding all markers ##
ChosenRes <- "res12"
seurat_merged <- SetIdent(seurat_merged, value = "RNA_snn_res.1.2")
allMarkers <- FindAllMarkers(seurat_merged, only.pos = F, min.pct = 0.1, logfc.threshold = 0.25,
                             max.cells.per.ident = 10000)
allMarkers <- group_by(allMarkers, cluster)
write.table(allMarkers, file = paste0(projectFolder, "seurat_merged_harmony_", ChosenRes, "_FindAllMarkers_minPCT01_logFC025.csv"),
            row.names = F, col.names = T, sep = "\t", quote = F)

###### 1.11: Annotation ###### 

# Re-annotate seurat_merged ## 
seurat_merged <- SetIdent(seurat_merged, value = "RNA_snn_res.1.2")

seurat_merged <- RenameIdents(seurat_merged,'0'='TAL') # highQ TAL 
seurat_merged <- RenameIdents(seurat_merged,'1'='PT') # low nFeat, higher MALAT1, PT
seurat_merged <- RenameIdents(seurat_merged,'2'='CNT') # highQ CNT 
seurat_merged <- RenameIdents(seurat_merged,'3'='DTL') # highQ DTL 
seurat_merged <- RenameIdents(seurat_merged,'4'='DCT') # highQ DCT 
seurat_merged <- RenameIdents(seurat_merged,'5'='TAL') # highQ TAL 
seurat_merged <- RenameIdents(seurat_merged,'6'='PT') # low nFeat, PT
seurat_merged <- RenameIdents(seurat_merged,'7'='CNT') # highQ CNT 
seurat_merged <- RenameIdents(seurat_merged,'8'='PT') # high nFeat, PT
seurat_merged <- RenameIdents(seurat_merged,'9'='EC_PTC_AVR') # highQ EC_PTC_AVR  
seurat_merged <- RenameIdents(seurat_merged,'10'='EC_GC') # highQ EC_GC 
seurat_merged <- RenameIdents(seurat_merged,'11'='IC_A') # highQ IC_A 
seurat_merged <- RenameIdents(seurat_merged,'12'='FIB') # highQ FIB 
seurat_merged <- RenameIdents(seurat_merged,'13'='POD') # high nFeat, highQ POD
seurat_merged <- RenameIdents(seurat_merged,'14'='remove') # high nFeat, lowQ
seurat_merged <- RenameIdents(seurat_merged,'15'='keepfornow') # high mito, likely to remove on subcluster level
seurat_merged <- RenameIdents(seurat_merged,'16'='MYEL_B') # highQ myeloid and B cell 
seurat_merged <- RenameIdents(seurat_merged,'17'='T_NKT') # highQ T_NKT 
seurat_merged <- RenameIdents(seurat_merged,'18'='IC_B') # highQ IC_B 
seurat_merged <- RenameIdents(seurat_merged,'19'='PEC') # highQ PEC 
seurat_merged <- RenameIdents(seurat_merged,'20'='keepfornow') # high nFeat, likely to remove on subcluster level
seurat_merged <- RenameIdents(seurat_merged,'21'='remove') # remove, lowQ EC 
seurat_merged <- RenameIdents(seurat_merged,'22'='VSMC_P') # highQ VSMC_P 
seurat_merged <- RenameIdents(seurat_merged,'23'='EC_AEA_DVR') # highQ EC_AEA_DVR
seurat_merged <- RenameIdents(seurat_merged,'24'='PT') # high mito, keep for now
seurat_merged <- RenameIdents(seurat_merged,'25'='MES') # highQ MES/REN
seurat_merged <- RenameIdents(seurat_merged,'26'='remove') # high nFeat, lowQ
seurat_merged <- RenameIdents(seurat_merged,'27'='PT') # high nFeat, keep for now
seurat_merged <- RenameIdents(seurat_merged,'28'='remove') # high doublet score, high nFeat, doublets
seurat_merged <- RenameIdents(seurat_merged,'29'='CNT') # high mito, keep for now
seurat_merged <- RenameIdents(seurat_merged,'30'='keepfornow') # likely lowQ
seurat_merged <- RenameIdents(seurat_merged,'31'='keepfornow') # high doublet score, high nFeat, keep for now, likely to remove on subcluster level
seurat_merged <- RenameIdents(seurat_merged,'32'='remove') # high mito, high nFeat, lowQ
seurat_merged <- RenameIdents(seurat_merged,'33'='remove') # high doublet score, lowQ
seurat_merged <- RenameIdents(seurat_merged,'34'='PT') # PT
seurat_merged <- RenameIdents(seurat_merged,'35'='keepfornow') # high doublet score, high nFeat, likely to remove on subcluster level
seurat_merged <- RenameIdents(seurat_merged,'36'='EC_others') # highQ EC_others
seurat_merged <- RenameIdents(seurat_merged,'37'='remove') # high doublet score, high nFeat
seurat_merged <- RenameIdents(seurat_merged,'38'='keepfornow') # to check on subcluster level
seurat_merged <- RenameIdents(seurat_merged,'39'='keepfornow') # high doublet score, high nFeat, likely to remove on subcluster level
seurat_merged <- RenameIdents(seurat_merged,'40'='keepfornow') # to check on subcluster level
seurat_merged <- RenameIdents(seurat_merged,'41'='remove') # high nFeat, no marker genes

DimPlot(seurat_merged, label=T)

## Cluster IDs in 'main_annot'-column
seurat_merged$main_annot <- Idents(seurat_merged)
seurat_merged$main_annot <- factor(seurat_merged$main_annot)
seurat_merged <- SetIdent(seurat_merged, value = "main_annot")

## Plotting the result
pdf(paste0(projectFolder, "seurat_merged_harmony_main_annot_dimPlot_umap_withLowQ.pdf"), width = 10, height = 10)
DimPlot(seurat_merged, reduction = "umap", group.by = "main_annot", label = TRUE, pt.size = 0.7, raster.dpi = c(1024,1024))
dev.off()

seurat_merged <- SetIdent(seurat_merged, value = "main_annot")
seurat_merged <- subset(seurat_merged, ident = c("remove"), invert = T) 

pdf(paste0(projectFolder, "seurat_merged_harmony_main_annot_dimPlot_umap_withoutLowQ.pdf"), width = 10, height = 10)
DimPlot(seurat_merged, reduction = "umap", group.by = "main_annot", label = TRUE, pt.size = 0.7, raster.dpi = c(1024,1024))
dev.off()

levels(seurat_merged$main_annot)
table(seurat_merged$main_annot)

## This concludes part 1



##### 2. PART 2 - Re-running pipeline after first lowQ removal #####

###### 2.1: Re-Running Seurat pipeline ###### 
projectFolder <- "path/HQ1/"

## Normalizing the data
seurat_merged <- NormalizeData(seurat_merged, normalization.method = "LogNormalize", scale.factor = 10000)

## Highly variable features
seurat_merged <- FindVariableFeatures(seurat_merged, selection.method = "vst", nfeatures = 2000)

## Scaling the data + regressing for mt.RNA, nCount, nFeature
seurat_merged <- ScaleData(seurat_merged, vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt")) # we only use variable features

## PCA - linear dimensional reduction
seurat_merged <- RunPCA(seurat_merged, features = VariableFeatures(object = seurat_merged), npcs = 50)

pdf(paste0(projectFolder, "seurat_merged_noHarmony_elbowPlot.pdf"), width = 10, height = 10)
ElbowPlot(seurat_merged, n = 50)
dev.off()

## Clustering
numPCs <- 30 
seurat_merged <- FindNeighbors(seurat_merged, dims = 1:numPCs)
res = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1, 1.2, 1.5, 1.7, 2)
seurat_merged <- FindClusters(seurat_merged, resolution = res)

seurat_merged <- RunUMAP(seurat_merged, dims = 1:numPCs)

pdf(paste0(projectFolder, "seurat_merged_noHarmony_dimPlot_umap_multiRes_nPC",numPCs,".pdf"), width = 10, height = 8)
for (i in 1:length(res)){
  resolution <- paste0("RNA_snn_res." , res[i])
  print(DimPlot(seurat_merged, reduction = "umap", group.by = resolution, label = T, pt.size = 0.7, raster.dpi = c(1024,1024)) + ggtitle(paste0("Resolution: ", res[i])))
}
dev.off()

###### 2.2: Harmony integration ###### 

## Using Harmony
seurat_merged <- RunHarmony(seurat_merged, "orig.ident")

## Clustering
numPCs <- 30
seurat_merged <- FindNeighbors(seurat_merged, reduction = "harmony", dims = 1:numPCs)
res = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1, 1.2, 1.5, 1.7, 2)
seurat_merged <- FindClusters(seurat_merged, resolution = res)
seurat_merged <- RunUMAP(seurat_merged, reduction = "harmony", dims = 1:numPCs)

pdf(paste0(projectFolder, "seurat_merged_Harmony_dimPlot_umap_multiRes_nPC",numPCs,".pdf"), width = 10, height = 8)
for (i in 1:length(res)){
  resolution <- paste0("RNA_snn_res." , res[i])
  print(DimPlot(seurat_merged, reduction = "umap", group.by = resolution, label = T, pt.size = 0.7, raster.dpi = c(1024,1024)) + ggtitle(paste0("Resolution: ", res[i])))
}
dev.off()

###### 2.3: Checking previous final MainAnnot object ###### 

## Plotting
seurat_merged <- SetIdent(seurat_merged, value = "subcluster_annot")

pdf(paste0(projectFolder, "seurat_merged_DimPlot_HQ_mapped.pdf"), width = 11, height = 10)
DimPlot(seurat_merged, group.by = "subcluster_annot", label = T, pt.size = 0.7, raster.dpi = c(1024,1024))
dev.off()

###### 2.4: Choosing resolution + QC per cluster ###### 

ChosenRes <- "res08" #proceed with res08
seurat_merged <- SetIdent(seurat_merged, value = "RNA_snn_res.0.8")
DimPlot(seurat_merged, label=T)

###### 2.5: Re-annotate res0.8 to proceed with subcluster clean-up ###### 
seurat_merged <- SetIdent(seurat_merged, value = "RNA_snn_res.0.8")

seurat_merged <- RenameIdents(seurat_merged,'0'='TAL')
seurat_merged <- RenameIdents(seurat_merged,'1'='DTL')
seurat_merged <- RenameIdents(seurat_merged,'2'='PT')
seurat_merged <- RenameIdents(seurat_merged,'3'='DCT')
seurat_merged <- RenameIdents(seurat_merged,'4'='CNT')
seurat_merged <- RenameIdents(seurat_merged,'5'='CNT')
seurat_merged <- RenameIdents(seurat_merged,'6'='PT')
seurat_merged <- RenameIdents(seurat_merged,'7'='PT')
seurat_merged <- RenameIdents(seurat_merged,'8'='EC')
seurat_merged <- RenameIdents(seurat_merged,'9'='TAL')
seurat_merged <- RenameIdents(seurat_merged,'10'='EC')
seurat_merged <- RenameIdents(seurat_merged,'11'='IC_A')
seurat_merged <- RenameIdents(seurat_merged,'12'='TAL')
seurat_merged <- RenameIdents(seurat_merged,'13'='FIB')
seurat_merged <- RenameIdents(seurat_merged,'14'='POD')
seurat_merged <- RenameIdents(seurat_merged,'15'='CNT')
seurat_merged <- RenameIdents(seurat_merged,'16'='B_MYEL')
seurat_merged <- RenameIdents(seurat_merged,'17'='IC_B')
seurat_merged <- RenameIdents(seurat_merged,'18'='T_NKT')
seurat_merged <- RenameIdents(seurat_merged,'19'='PEC')
seurat_merged <- RenameIdents(seurat_merged,'20'='VSMC_P')
seurat_merged <- RenameIdents(seurat_merged,'21'='EC')
seurat_merged <- RenameIdents(seurat_merged,'22'='MES')
seurat_merged <- RenameIdents(seurat_merged,'23'='PT')
seurat_merged <- RenameIdents(seurat_merged,'24'='PT')
seurat_merged <- RenameIdents(seurat_merged,'25'='FIB')
seurat_merged <- RenameIdents(seurat_merged,'26'='IC_A')
seurat_merged <- RenameIdents(seurat_merged,'27'='CNT')
seurat_merged <- RenameIdents(seurat_merged,'28'='unknown')
seurat_merged <- RenameIdents(seurat_merged,'29'='PT')
seurat_merged <- RenameIdents(seurat_merged,'30'='EC')
seurat_merged <- RenameIdents(seurat_merged,'31'='CNT')
seurat_merged <- RenameIdents(seurat_merged,'32'='CNT')

DimPlot(seurat_merged, label=T)

## Cluster IDs in 'main_annot_HQ'-column
seurat_merged$main_annot_HQ <- Idents(seurat_merged)
seurat_merged$main_annot_HQ <- factor(seurat_merged$main_annot_HQ)
seurat_merged <- SetIdent(seurat_merged, value = "main_annot_HQ")

pdf(paste0(projectFolder, "seurat_merged_harmony_main_annot_HQ_dimPlot_umap.pdf"), width = 10, height = 8)
DimPlot(seurat_merged, group.by = "main_annot_HQ", reduction = "umap", label = TRUE, pt.size = 0.7, raster.dpi = c(1024,1024)) + ggtitle("seurat_merged_harmony_res08_main_annot_HQ")
dev.off()


###### 2.6: Plotting markers genes w/ featureplots and dotplots (KidneyMarkersV2) ###### 

# Load in marker genes (KidneyMarkersV2) (Proprietary list, curated from the literature)

# Podocytes (POD):
POD <- c("NPHS1", "NPHS2", "PLA2R1", "PCOLCE2", "WT1", "PODXL", "PTPRO", "PTPRQ", "CLIC5")

# Mesangial cells (MES):
MES <- c("PDGFRB", "AKAP12", "COL1A1", "COL1A2", "LUM", "ACTN1", "EBF1", "ROBO1", "ITGA8", "POSTN", "TAGLN")

# Parietal epithelial cell (PEC):
PEC <- c("PAX2", "CLDN1","CRB2","VCAM1","RBFOX1","ALDH1A2","CFH")

# Proximal tubule cell (PT):
PT <- c("SLC34A1", "LRP2", "CUBN", "SLC13A1", "SLC13A2")
PT_S1 <- c("SLC5A2", "PRODH2", "SLC22A8")
PT_S1S2 <- c("SLC5A12", "SLC13A3", "SLC22A6", "SLC17A3")
PT_S2 <- c("SLC34A1")
PT_S2S3 <- c("SLC22A7")
PT_S3 <- c("MOGAT1", "SLC5A11", "SLC7A13", "SLC5A8", "ABCC3", "SATB2", "AGT", "SLC16A9")

# Loop of Henle, descending thin limb (DTL)
DTL <- c("IRX3", "FOXC1")
DTL_1 <- c("SATB2", "JAG1", "ADGRL3", "ID1")
DTL_2 <- c("VCAM1", "SLC39A8", "AQP1", "LRRC4C", "LRP2", "UNC5D", "SATB2")
DTL_3 <- c("SLC14A2", "SMOC2")

# Loop of Henle, descending thin limb cell type 3 - ascending thin limb (overlap)
DTL3_ATL <- c("CLDN1", "CLDN3", "AKR1B1", "CLDN4", "BCL6", "SH3GL3")

# Loop of Henle, ascending thin limb (ATL)
ATL <- c("BCAS1", "CLCNKA", "PROX1")

# Loop of Henle, ascending thin limb - thick ascending limb (overlap)
ATL_TAL <- c("CLDN10")

# Loop of Henle, thick ascending limb (TAL)
TAL <- c("SLC12A1", "UMOD", "CLDN16", "EGF")

# Macula densa cell (MD)
MD <- c("SLC12A1", "NOS1", "ROBO2", "CALCR", "PPFIA2", "PAPPA2", "OXTR", "PTGS2")

# Distal convoluted tubule (DCT)
DCT <- c("SLC12A3", "CNNM2", "FGF13", "KLHL3", "LHX1", "TRPM6")

# Connecting tubule cell (CNT)
CNT <- c("AQP2", "SCNN1G", "SLC8A1", "CALB1")

# Collecting duct: principal cell - inner medullary collecting duct cell 
PC_IMCD <- c("AQP2", "AQP3", "GATA3", "FXYD4", "SOX5")

# Collecting duct, principal cell (PC)
PC <- c("AQP2", "AQP3", "SCNN1G", "SCNN1B")

# Inner medullary collecting duct (IMCD)
IMCD <- c("SLC14A2", "PHACTR1", "PCDH7","HS3ST5")

# Collecting duct, intercalated cell (IC)
IC <- c("ATP6V1C2", "ATP6V0D2", "TMEM213","CLNK")

# Collecting duct, type A intercalated cell (A_IC)
IC_A <- c("SLC4A1", "SLC26A7", "KIT", "AQP6")

# Collecting duct, type B intercalated cell (B_IC)
IC_B <- c("SLC26A4", "SLC4A9")

# Endothelial cells, NOS (EC); PECAM1 is pan-endothelium
EC <- c("PECAM1", "KDR", "CDH5", "PTPRB", "SEMA3G")

# Glomerular endothelial cells, capillary (EC-GC)
EC_GC <- c("EHD3", "SMAD6", "EMCN", "HECW2", "SOST")

# Afferent / Efferent Arteriole Endothelial Cell - Descending Vasa Recta Endothelial Cell (overlap)
EC_AEA_DVR <- c("FBLN5", "CXCL12", "CLDN5", "SOX17")

# Descending Vasa Recta Endothelial Cell
EC_DVR <- c("AQP1", "SLC14A1")

# Peritubular Capilary Endothelial Cell - Ascending Vasa Recta Endothelial Cell (overlap) i.e. fenestrated endothelium
EC_PTC_AVR <- "PLVAP"

# Ascending Vasa Recta Endothelial Cell
EC_AVR <- "NR2F2"

# Lymphatic Cell
EC_lym <- c("MMRN1", "CD36", "TBX1", "PKHD1L1", "PROX1")

# Vascular smooth muscle cell (VSMC)
VSMC <- c("ACTA2", "TAGLN")

# Fibroblast (FIB)
FIB <- c("COL1A2", "DCN", "LUM", "C7", "TNC", "COL1A1", "EMILIN1", "MMP2")

# Myofibroblast (MYOF)
MYOF <- c("PDGFRB", "ACTA2")

# All immune cells (pan-immune), IMMUNE
IMMUNE <- "PTPRC"

# NK cell, NKT cell, CD8 (NK_T)
NKT <- c("FCGR3A", "GNLY", "GZMA", "PRF1", "GZMH", "NKG7")

# CD8 T cell
CD8_T <- c("CD8A", "GZMK")

# CD4 T cell
CD4_T <- "CD4"

# T cell
T_cell <- c("IL7R", "CD3D", "CD3E", "CD3G")

# B cell
B_cell <- c("MS4A1", "CD19", "IGKC", "IGHA1", "IGHM", "IGHG1", "CD79A", "CD79B")

# Dendritic cell (DC)
DC <- c("CLEC9A", "XCR1", "BATF3", "IDO1", "CLEC10A")

# Plasmacytoid dendritic cell (pDC)
pDC <- c("IL3RA", "CLEC4C" ,"GZMB")

# Macrophage (MAC)
MAC <-c("ITGAM", "CD14", "CD68", "CD163")

# Classical monocyte (MONO)
MONO <- c("S100A8", "S100A9", "S100A12")

# neutrophil/granulocyte 
GRAN <- c("FCGR3B", "AQP9", "MNDA")

# Mast cell
MAST <- c("GATA2","CPA3","KIT")

# Pelvic epithelium and/or transitional epithelium of ureter
URO <- c("SAA2", "UPK1A", "UPK1B", "UPK3A", "S100P", "DHRS2")

# Make a list of the celltypes
KidneyMarkers <- list(
  POD=POD,
  MES=MES,
  PEC=PEC,
  PT=PT,
  PT_S1=PT_S1,
  PT_S1S2=PT_S1S2,
  PT_S2=PT_S2,
  PT_S2S3=PT_S2S3,
  PT_S3=PT_S3,
  DTL=DTL,
  DTL_1=DTL_1,
  DTL_2=DTL_2,
  DTL_3=DTL_3,
  DTL3_ATL=DTL3_ATL,
  ATL=ATL,
  ATL_TAL=ATL_TAL,
  TAL=TAL,
  MD=MD,
  DCT=DCT,
  CNT=CNT,
  PC_IMCD=PC_IMCD,
  PC=PC,
  IMCD=IMCD,
  IC=IC,
  IC_A=IC_A,
  IC_B=IC_B,
  EC=EC,
  EC_GC=EC_GC,
  EC_AEA_DVR=EC_AEA_DVR,
  EC_DVR=EC_DVR,
  EC_PTC_AVR=EC_PTC_AVR,
  EC_AVR=EC_AVR,
  EC_lym=EC_lym,
  VSMC=VSMC,
  FIB=FIB,
  MYOF=MYOF,
  IMMUNE=IMMUNE,
  NKT=NKT,
  CD8_T=CD8_T,
  CD4_T=CD4_T,
  T_cell=T_cell,
  B_cell=B_cell,
  DC=DC,
  pDC=pDC,
  MAC=MAC,
  MONO=MONO,
  GRAN=GRAN,
  MAST=MAST,
  URO=URO)

# Featureplots
pdf(paste0(projectFolder, "seurat_merged_harmony_KidneyMarkers_featureplot.pdf"), width = 10, height = 10)
for (i in 1:length(KidneyMarkers)){
  plot_subtitle <- paste0("Celltype ", names(KidneyMarkers)[[i]])
  
  # Create a new page with the additional title
  plot.new()
  text(x = 0.5, y = 0.95, plot_subtitle, cex = 1.2, font = 2, adj = 0.5)
  
  # Create the feature plot with the original title
  print(FeaturePlot(seurat_merged, features = KidneyMarkers[[i]], 
                    min.cutoff='q5',
                    max.cutoff='q95',
                    order=TRUE,
                    label=TRUE, label.size=3, repel=TRUE,
                    coord.fixed=TRUE ) )
}
dev.off()

# Featureplots with "order=F"
pdf(paste0(projectFolder, "seurat_merged_harmony_KidneyMarkers_featureplot_orderFALSE.pdf"), width = 10, height = 10)
for (i in 1:length(KidneyMarkers)){
  plot_subtitle <- paste0("Celltype ", names(KidneyMarkers)[[i]])
  
  # Create a new page with the additional title
  plot.new()
  text(x = 0.5, y = 0.95, plot_subtitle, cex = 1.2, font = 2, adj = 0.5)
  
  # Create the feature plot with the original title
  print(FeaturePlot(seurat_merged, features = KidneyMarkers[[i]], 
                    min.cutoff='q5',
                    max.cutoff='q95',
                    order=F,
                    label=TRUE, label.size=3, repel=TRUE,
                    coord.fixed=TRUE ) )
}
dev.off()

# DotPlots
pdf(paste0(projectFolder, "seurat_merged_harmony_KidneyMarkers_DotPlot.pdf"), width = 10, height = 10)
for (i in 1:length(KidneyMarkers)){
  plot_subtitle <- paste0("Celltype ", names(KidneyMarkers)[[i]])
  
  # Create a new page with the additional title
  plot.new()
  text(x = 0.5, y = 0.95, plot_subtitle, cex = 1.2, font = 2, adj = 0.5)
  
  # Create the feature plot with the original title
  print(DotPlot(seurat_merged, features = KidneyMarkers[[i]]) + RotatedAxis()) 
}
dev.off()


###### 2.7: Plotting markers genes w/ featureplots and dotplots (KidneyMarkersKPMP) ###### 

# Make a list of the celltypes, based on marker genes from:
# Lake, B.B., Menon, R., Winfree, S. et al. An atlas of healthy and injured cell states and niches in the human kidney. Nature 619, 585–594 (2023). DOI: 10.1038/s41586-023-05769-3

POD <- c("PTPRQ", "WT1", "NTNG1", "NPHS1", "NPHS2", "CLIC5", "PODXL")
dPOD <- c("CDKN1C", "SPOCK2", "B2M", "CD81", "S100A6", "IGFBP2", "PTPRQ") # PTPRQ downregulated in degenerative state

PEC <- c("CLDN1", "VCAM1", "CFH", "RBFOX1", "ALDH1A2")

PT <- c("LRP2", "CUBN", "SLC13A1")
dPT <- c("APOE", "PDZK1IP1", "CCNI", "GATM", "ALDOB")
aPT <- c("ITGB8", "CDH6", "DCDC2", "TPM1", "VCAM1", "DLGAP1", "ACSM3", "KIF26B", "HAVCR1")
cPT <- c("APOLD1", "CENPF", "DIAPH3", "MKI67", "TOP2A")

PT_S1 <- c("SLC5A12", "SLC13A3", "SLC22A6", "PRODH2", "SLC5A2", "SLC22A8")
PT_S2 <- c("SLC5A12", "SLC13A3", "SLC22A6", "SLC34A1", "SLC22A7")
PT_S3 <- c("SLC22A7", "MOGAT1", "SLC5A11", "SLC22A24", "SLC7A13", "SLC5A8", "ABCC3", "SATB2")

TL <- c("CRYAB", "TACSTD2", "SLC44A5", "KLRG2", "COL26A1", "BOC")

DTL1 <- c("SATB2", "JAG1", "ADGRL3", "ID1", "CLDN10", "AQP1") #CLDN10, AQP1 negative markers
DTL2 <- c("VCAM1", "SLC39A8", "AQP1", "LRRC4C", "LRP2", "UNC5D", "SATB2", "CLDN10") #CLDN10 negative marker
DTL3 <- c("CLDN1", "AKR1B1", "CLDN4", "BCL6", "SH3GL3", "SLC14A2", "SMOC2", "CLDN10", "AQP2") # CLDN10, AQP2 negative markers 
dDTL3 <- c("SERPINA1", "SPP1", "NNMT", "USP43")

ATL <- c("CLDN1", "AKR1B1", "CLDN4", "BCL6", "SH3GL3", "BCAS1", "CLCNKA", "CLDN10", "PROX1")
dATL <- c("SERPINA1", "SPP1", "S100A4", "FTL", "FTH1", "DEFB1")

TAL <- c("CASR", "SLC12A1", "UMOD")
dTAL <- c("SPP1", "S100A6", "TPT1", "ACTB", "WFDC2", "UMOD", "EGF") # UMOD, EGF downregulated in degenerative state
aTAL <- c("ITGB6", "PROM1")

M_TAL <- c("NELL1", "ESRRB", "EGF", "CLDN14", "PROX1", "MFSD4A", "KCTD16", "RAP1GAP", "ANK2", "CYFIP2")
dM_TAL <- c("DEFB1", "SFRP1", "KNG1")
C_TAL <- c("NELL1", "ESRRB", "EGF", "PPM1E", "GP2", "ENOX1", "TMEM207", "TMEM52B", "CLDN16", "WNK1")
MD <- c("NOS1", "ROBO2", "CALCR", "PPFIA2", "PAPPA2", "UMOD") # UMOD low in MD

DCT <- c("SLC12A3", "CNNM2", "FGF13", "KLHL3", "LHX1", "TRPM6", "CASR") #CASR low
dDCT <- c("SPP1", "DEFB1", "S100A6", "FTL", "FTH1", "CST3")
cDCT <- c("DIAPH3", "CENPF", "ANLN", "APOLD1", "MKI67", "TOP2A")
DCT1 <- c("TRPM7", "ADAMTS17", "ITPKB", "ZNF385D", "HS6ST2", "SLC8A1") # SLC8A1 negative marker
dDCT1 <- c("SFRP1", "CA2", "FXYD2")
DCT2 <- c("TRPV5", "SLC8A1", "SCN2A", "HSD11B2", "CALB1", "SLC12A3") # SLC12A3 negative marker

CNTg <- c("SLC8A1", "SCN2A", "HSD11B2", "CALB1", "SLC12A3") #g for general, SLC12A3 negative marker
dCNTg <- c("SPP1", "DEFB1", "S100A6", "FTL", "FTH1", "CST3")
cCNTg <- c("DIAPH3", "CENPF", "ANLN", "APOLD1", "MKI67", "TOP2A")
CNT <- c("KITLG", "PCDH7")
dCNT <- c("PIGR", "S100A11", "CLU", "WFDC2")
CNT_PC <- c("RALYL", "TOX", "SGPP1", "SCNN1G", "SCNN1B", "KCNIP1")

PC <- c("GATA3", "AQP2", "AQP3")
dPC <- c("B2M", "FTH1", "FTL", "DEFB1", "CD24")
CCD_PC <- c("SCNN1G", "SCNN1B", "FXYD4", "SOX5", "PDE10A", "SLC25A29", "ST6GAL1", "PAPPA")
OMCD_PC <- c("SCNN1G", "SCNN1B", "FXYD4", "SOX5", "SYK", "FAM81A", "PROM1", "KCNK13")
dOMCD_PC <- c("PTGER1", "PCSK1N", "FGL2", "EPCAM", "CDH16", "PSAP")
IMCD <- c("FXYD4", "SOX5", "PHACTR1", "PCDH7", "SLC14A2", "HS3ST5")
dIMCD <- c("PTGS1", "KCNK3", "AOC1", "CLDN8", "CLU", "MMP7", "SLPI", "CRYAB", "S100A11", "S100A6")
PapE <- c("TACSTD2", "TP63", "GPX2", "FXYD3", "KRT5")

IC <- c("ATP6V0D2", "ATP6V1C2", "TMEM213", "CLNK")
dIC <- c("CLU", "B2M", "FTH1", "FTL", "DEFB1", "CD24")
CCD_IC_A <- c("SLC4A1", "SLC26A7", "HS6ST3", "NXPH2", "LEF1", "ADGRF5", "CALCA") #CALCA negative marker
dCCD_IC_A <- c("APOE", "DSG2", "CKB")
CNT_IC_A <- c("SLC4A1", "SLC26A7", "SLC8A1", "SCN2A", "CALB1")
OMCD_IC_A <- c("SLC4A1", "SLC26A7", "KIT", "AQP6", "STAP1", "FAM184B", "CALCA")
IC_B <- c("SLC4A9", "SLC35F3", "SLC26A4", "INSRR", "TLDC2", "KIT", "SLC4A1", "CALCA") #KIT, SLC4A1, CALCA negative markers

EC <- c("CD34", "PECAM1", "PTPRB", "MEIS2", "FLT1", "EMCN")
dEC <- c("B2M", "TMSB4X", "TMSB10", "HLA-B", "IGFBP5")
cEC <- c("APOLD1", "DIAPH3", "EZH2", "CENPF", "MKI67", "TOP2A")
EC_GC <- c("EMCN", "HECW2", "PLAT", "ITGA8", "EHD3", "KDR", "SOST", "NR2F2") #NR2F2 negative marker
EC_AEA <- c("BTNL9", "ADAMTS6", "PALMD", "AQP1", "TM4SF1", "VEGFC", "CCDC3", "CDH5", "SERPINE2", "FBLN5", "CXCL12", "SOX17", "NR2F2", "SOST") # NR2F2, SOST negative markers
EC_DVR <- c("BTNL9", "ADAMTS6", "PALMD", "AQP1", "TM4SF1", "MCTP1", "SLC14A1", "ENPP2", "LYPD6B", "NR2F2", "IGFBP5", "KDR", "SOST") #NR2F2, IGFBP5, KDR, SOST negative markers
EC_PTC <- c("CEACAM1", "DNASE1L3", "PLVAP", "PITPNC1", "GRB10", "SLCO2A1", "RAPGEF4", "NR2F2", "SOST") #NR2F2, SOST negative markers
dEC_PTC <- c("FTL", "DHFR", "SPP1", "GPX3", "PEBP1")
EC_AVR <- c("CEACAM1", "DNASE1L3", "PLVAP", "GPM6A", "EDIL3", "TLL1", "ZNF385D", "NR2F2", "PALMD", "BTNL9", "SLC14A1", "SOST") # PALMD, BTNL9, SLC14A1, SOST negative markers
EC_LYM <- c("MMRN1", "CD36", "TBX1", "PKHD1L1", "PROX1")

VSM_Pg <- c("NOTCH3", "PDGFRB", "ITGA8")
dVSM_Pg <- c("TAGLN", "FLNA", "MYL9", "ACTA2", "DSTN","PDGFRB") #PDGFRB downregulated in degenerative state
MC <- c("PIP5K1B", "ROBO1", "PIEZO2", "DAAM2", "PHTF2", "GATA3", "POSTN")
REN <- c("PIP5K1B", "ROBO1", "REN", "PDE10A", "ABCC8", "COL13A1", "GRID2")
VSMC <- c("NTRK3", "MYH11", "RGS6", "ADRA1A", "LDB3", "MCAM")
VSMC_P <- c("NTRK3", "CCDC102B", "RGS5", "ABCC9", "ADCY3", "ADGRB3")

FIBg <- c("COL1A1", "COL1A2", "C7", "NEGR1", "FBLN5", "DCN", "CDH11")
FIB <- c("LAMA2", "GGT5", "LUM", "AEBP1", "C1S", "SFRP1", "MEG3", "CXCL12")
dFIB <- c("IGFBP7", "B2M", "TPT1", "VIM", "C1R")
aFIB <- c("FLRT2", "FGF14", "PRRX1", "NAV3", "ABI3BP", "IGF1")

M_FIB <- c("SYT1", "TNC", "PLCXD3", "GABRG3", "GREB1L", "KCNK2") 
dM_FIB <- c("S100A6", "S100A11", "KRT19", "SLPI", "AKR1B1")
MYOF <- c("SYNPO2", "PCDH7", "KCNMA1", "LMOD1", "TTLL7", "DTNA", "COL14A1")
cMYOF <- c("DIAPH3", "EZH2", "BRIP1", "MELK", "TPX2", "MKI67", "TOP2A")

IMM <- "PTPRC"
cIMM <- c("DIAPH3", "APOLD1", "CENPF", "MKI67", "TOP2A")
Bcell <- c("BANK1", "BLK", "MS4A1", "BACH2")
PL <- c("IGKC", "TENT5C", "MZB1", "FCRL5", "CD38", "JCHAIN")
Tcell <- c("CD96", "CD247", "THEMIS", "BCL11B", "CAMK4", "IL7R")
NKT <- c("CD96", "CD247", "RUNX3", "GNLY", "NKG7", "CCL5", "KLRF1", "CCL4", "GZMA")
MAST <- c("MS4A2", "CPA3", "KIT", "IL3RA") # IL3RA negative marker
MAC_M2 <- c("F13A1", "MRC1", "CD163", "STAB1", "SLC1A3", "CD14")
cDC <- c("ITGAX", "HLA-DQA1", "HLA-DRA", "CSF2RA", "CIITA", "WDFY4", "FLT3", "ZNF366", "CADM1", "ZBTB46", "CLEC9A", "CD14") # CD14 negative marker
pDC <- c("IRF8", "CUX2", "P2RY14", "IL3RA", "CLEC4C", "CD14") # CD14 negative marker
ncMON <- c("CTSS", "IRAK3", "TCF7L2", "TNFRSF1B", "FCN1", "HLA-DRA", "FCGR3A")
NC <- c("S100A9", "S100A8", "IFITM2", "FCGR3B", "CD1C")
SC_NEU <- c("CDH19", "NRXN1", "GINS3")

KidneyMarkersKPMP <- list(
  POD=POD,
  dPOD=dPOD,
  PEC=PEC,
  PT=PT,
  dPT=dPT,
  aPT=aPT,
  cPT=cPT,
  PT_S1=PT_S1,
  PT_S2=PT_S2,
  PT_S3=PT_S3,
  TL=TL,
  DTL1=DTL1,
  DTL2=DTL2,
  DTL3=DTL3,
  dDTL3=dDTL3,
  ATL=ATL,
  dATL=dATL,
  TAL=TAL,
  dTAL=dTAL,
  aTAL=aTAL,
  M_TAL=M_TAL,
  dM_TAL=dM_TAL,
  C_TAL=C_TAL,
  MD=MD,
  DCT=DCT,
  dDCT=dDCT,
  cDCT=cDCT,
  DCT1=DCT1,
  dDCT1=dDCT1,
  DCT2=DCT2,
  CNTg=CNTg,
  dCNTg=dCNTg,
  cCNTg=cCNTg,
  CNT=CNT,
  dCNT=dCNT,
  CNT_PC=CNT_PC,
  PC=PC,
  dPC=dPC,
  CCD_PC=CCD_PC,
  OMCD_PC=OMCD_PC,
  dOMCD_PC=dOMCD_PC,
  IMCD=IMCD,
  dIMCD=dIMCD,
  PapE=PapE,
  IC=IC,
  dIC=dIC,
  CCD_IC_A=CCD_IC_A,
  dCCD_IC_A=dCCD_IC_A,
  CNT_IC_A=CNT_IC_A,
  OMCD_IC_A=OMCD_IC_A,
  IC_B=IC_B,
  EC=EC,
  dEC=dEC,
  cEC=cEC,
  EC_GC=EC_GC,
  EC_AEA=EC_AEA,
  EC_DVR=EC_DVR,
  EC_PTC=EC_PTC,
  dEC_PTC=dEC_PTC,
  EC_AVR=EC_AVR,
  EC_LYM=EC_LYM,
  VSM_Pg=VSM_Pg,
  dVSM_Pg=dVSM_Pg,
  MC=MC,
  REN=REN,
  VSMC=VSMC,
  VSMC_P=VSMC_P,
  FIBg=FIBg,
  FIB=FIB,
  dFIB=dFIB,
  aFIB=aFIB,
  M_FIB=M_FIB,
  dM_FIB=dM_FIB,
  MYOF=MYOF,
  cMYOF=cMYOF,
  IMM=IMM,
  cIMM=cIMM,
  Bcell=Bcell,
  PL=PL,
  Tcell=Tcell,
  NKT=NKT,
  MAST=MAST,
  MAC_M2=MAC_M2,
  cDC=cDC,
  pDC=pDC,
  ncMON=ncMON,
  NC=NC,
  SC_NEU=SC_NEU)


# Featureplots
pdf(paste0(projectFolder, "seurat_merged_harmony_KidneyMarkersKPMP_featureplot.pdf"), width = 10, height = 10)
for (i in 1:length(KidneyMarkersKPMP)){
  plot_subtitle <- paste0("Celltype ", names(KidneyMarkersKPMP)[[i]])
  
  # Create a new page with the additional title
  plot.new()
  text(x = 0.5, y = 0.95, plot_subtitle, cex = 1.2, font = 2, adj = 0.5)
  
  # Create the feature plot with the original title
  print(FeaturePlot(seurat_merged, features = KidneyMarkersKPMP[[i]], 
                    min.cutoff='q5',
                    max.cutoff='q95',
                    order=TRUE,
                    label=TRUE, label.size=3, repel=TRUE,
                    coord.fixed=TRUE ) )
}
dev.off()

# Featureplots with "order=F"
pdf(paste0(projectFolder, "seurat_merged_harmony_KidneyMarkersKPMP_featureplot_orderFALSE.pdf"), width = 10, height = 10)
for (i in 1:length(KidneyMarkersKPMP)){
  plot_subtitle <- paste0("Celltype ", names(KidneyMarkersKPMP)[[i]])
  
  # Create a new page with the additional title
  plot.new()
  text(x = 0.5, y = 0.95, plot_subtitle, cex = 1.2, font = 2, adj = 0.5)
  
  # Create the feature plot with the original title
  print(FeaturePlot(seurat_merged, features = KidneyMarkersKPMP[[i]], 
                    min.cutoff='q5',
                    max.cutoff='q95',
                    order=FALSE,
                    label=TRUE, label.size=3, repel=TRUE,
                    coord.fixed=TRUE ) )
}
dev.off()

# DotPlots
pdf(paste0(projectFolder, "seurat_merged_harmony_KidneyMarkersKPMP_DotPlot.pdf"), width = 10, height = 10)
for (i in 1:length(KidneyMarkersKPMP)){
  plot_subtitle <- paste0("Celltype ", names(KidneyMarkersKPMP)[[i]])
  
  # Create a new page with the additional title
  plot.new()
  text(x = 0.5, y = 0.95, plot_subtitle, cex = 1.2, font = 2, adj = 0.5)
  
  # Create the feature plot with the original title
  print(DotPlot(seurat_merged, features = KidneyMarkersKPMP[[i]]) + RotatedAxis()) 
}
dev.off()


###### 2.8: FindAllMarkers + plotting of DEGs ###### 

## Finding all markers ##
ChosenRes <- "main_annot_HQ"
seurat_merged <- SetIdent(seurat_merged, value = "main_annot_HQ")

allMarkers <- FindAllMarkers(seurat_merged, only.pos = F, min.pct = 0.1, logfc.threshold = 0.25,
                             max.cells.per.ident = 10000)

allMarkers <- group_by(allMarkers, cluster)
write.table(allMarkers, file = paste0(projectFolder, "seurat_merged_harmony_", ChosenRes, "_FindAllMarkers_minPCT01_logFC025.csv"),
            row.names = F, col.names = T, sep = "\t", quote = F)

###### 2.9: Mapping kidney cell atlas object with SingleR ###### 

# SOURCE: https://www.kpmp.org/doi-collection/10-48698-yyvc-ak78
# FILENAME: WashU-UCSD_HuBMAP_KPMP-Biopsy_10X-R_12032021.h5Seurat
# NAME: ATLAS EXPLORER V1.3 SINGLE-NUCLEUS RNA-SEQ DATA
# DESCRIPTION: Aggregated, clustered single-nucleus RNA-seq data used in the KPMP Atlas Explorer v1.3 - Dataset published 2021 via Kidney Precision Medicine Project
# CREATOR(S): Kidney Precision Medicine Project
# DATE: December 8, 2021
# LICENSE: Creative Commons Attribution 4.0 International (CC BY 4.0)
# DOI: https://doi.org/10.48698/yyvc-ak78
# HOW TO CITE THIS DATA:Kidney Precision Medicine Project. (2021). Aggregated, clustered single-nucleus RNA-seq data used in the KPMP Atlas Explorer v1.3. Kidney Precision Medicine Project. https://doi.org/10.48698/yyvc-ak78

kidneyAtlas <- LoadH5Seurat(file = paste0(getwd(), "/data/KPMP_kidney_atlas/kpmp_kidney_atlas.h5Seurat"))

seurat_matrix <- GetAssayData(seurat_merged, assay = "RNA", slot = "data")

kidneyAtlas <- SetIdent(kidneyAtlas, value = kidneyAtlas@meta.data$subclass.l1)
kidneyAtlas_matrix <- GetAssayData(kidneyAtlas, assay = "RNA", slot = "data")

kidneyAtlas_prediction_l1 <- SingleR(test = seurat_matrix, ref = kidneyAtlas_matrix,
                                     assay.type.test=1, labels = kidneyAtlas$subclass.l1)

saveRDS(kidneyAtlas_prediction_l1, file = paste0(projectFolder, "kidneyAtlas_predictions_l1.Rds"))

kidneyAtlas <- SetIdent(kidneyAtlas, value = kidneyAtlas@meta.data$subclass.l2)
kidneyAtlas_matrix <- GetAssayData(kidneyAtlas, assay = "RNA", slot = "data")

kidneyAtlas_prediction_l2 <- SingleR(test = seurat_matrix, ref = kidneyAtlas_matrix,
                                     assay.type.test=1, labels = kidneyAtlas$subclass.l2)

saveRDS(kidneyAtlas_prediction_l2, file = paste0(projectFolder, "kidneyAtlas_predictions_l2.Rds"))

seurat_merged[["kidneyAtlas_prediction_l1"]] <- plyr::mapvalues(x = colnames(seurat_merged), from = rownames(kidneyAtlas_prediction_l1), to = kidneyAtlas_prediction_l1$labels)
seurat_merged[["kidneyAtlas_prediction_l2"]] <- plyr::mapvalues(x = colnames(seurat_merged), from = rownames(kidneyAtlas_prediction_l2), to = kidneyAtlas_prediction_l2$labels)

# Checking kidney atlas projection l1
seurat_merged$kidneyAtlas_prediction_l1 <- factor(seurat_merged$kidneyAtlas_prediction_l1)
levels(seurat_merged$kidneyAtlas_prediction_l1)

pdf(paste0(projectFolder, "seurat_merged_harmony_dimplot_KidneyAtlas_prediction_l1.pdf"), width = 10, height = 10)
DimPlot(seurat_merged, group.by = "kidneyAtlas_prediction_l1", label = T, raster=FALSE) + NoLegend()
dev.off()

# Checking kidney atlas projection l2
seurat_merged$kidneyAtlas_prediction_l2 <- factor(seurat_merged$kidneyAtlas_prediction_l2)
levels(seurat_merged$kidneyAtlas_prediction_l2)

pdf(paste0(projectFolder, "seurat_merged_harmony_dimplot_KidneyAtlas_prediction_l2.pdf"), width = 10, height = 10)
DimPlot(seurat_merged, group.by = "kidneyAtlas_prediction_l2", label = T, raster=FALSE) + NoLegend()
dev.off()

###### 2.10: Saving cells for sub-clustering ###### 

seurat_merged <- SetIdent(seurat_merged, value = "main_annot_HQ")

# LOH, DCT, CNT
seurat_LOH_DCT_CNT <- subset(seurat_merged, ident = c("CNT", "DCT", "TAL")) 
saveRDS(seurat_LOH_DCT_CNT, file = paste0(projectFolder, "seurat_LOH_DCT_CNT.Rds"))
#seurat_LOH_DCT_CNT <- readRDS(paste0(projectFolder, "seurat_LOH_DCT_CNT.Rds"))

# PT (and DTL)
seurat_PT <- subset(seurat_merged, ident = c("PT", "DTL")) 
saveRDS(seurat_PT, file = paste0(projectFolder, "seurat_PT.Rds"))
#seurat_PT <- readRDS(paste0(projectFolder, "seurat_PT.Rds"))

# PEC
seurat_PEC <- subset(seurat_merged, ident = c("PEC")) 
saveRDS(seurat_PEC, file = paste0(projectFolder, "seurat_PEC.Rds"))
#seurat_PEC <- readRDS(paste0(projectFolder, "seurat_PEC.Rds"))

# IC_A
seurat_IC_A <- subset(seurat_merged, ident = c("IC_A")) 
saveRDS(seurat_IC_A, file = paste0(projectFolder, "seurat_IC_A.Rds"))
#seurat_IC_A <- readRDS(paste0(projectFolder, "seurat_IC_A.Rds"))

# IC_B
seurat_IC_B <- subset(seurat_merged, ident = c("IC_B")) 
saveRDS(seurat_IC_B, file = paste0(projectFolder, "seurat_IC_B.Rds"))
#seurat_IC_B <- readRDS(paste0(projectFolder, "seurat_IC_B.Rds"))

# FIB, MES, IMM
seurat_FIB_MES_IMM <- subset(seurat_merged, ident = c("MES", "VSMC_P", "FIB", "T_NKT", "B_MYEL")) 
saveRDS(seurat_FIB_MES_IMM, file = paste0(projectFolder, "seurat_FIB_MES_IMM.Rds"))
#seurat_FIB_MES_IMM <- readRDS(paste0(projectFolder, "seurat_FIB_MES_IMM.Rds"))

# EC
seurat_EC <- subset(seurat_merged, ident = c("EC")) 
saveRDS(seurat_EC, file = paste0(projectFolder, "seurat_EC.Rds"))
#seurat_EC <- readRDS(paste0(projectFolder, "seurat_EC.Rds"))

# POD
seurat_POD <- subset(seurat_merged, ident = c("POD")) 
saveRDS(seurat_POD, file = paste0(projectFolder, "seurat_POD.Rds"))
#seurat_POD <- readRDS(paste0(projectFolder, "seurat_POD.Rds"))

## This concludes part 2
## Subsequently, all subclusters are analyzed individually with additional rounds of QC, to generate highQ 'clean' Seurat-objects, required for part 3.



##### 3. PART 3 - Re-map clean subclusters #####
projectFolder <- "path/HQ2/"

seurat_merged <- SetIdent(seurat_merged, value = "main_annot_HQ")

###### 3.1: Annotating 'unknown' cluster ###### 

# high ribosomal RNA when checking marker-gene list
# when checking dotplots: degenerative, some thin limb, but also contamination in general -> remove cluster "unknown".

###### 3.2: Remapping clean/HQ POD cluster (seurat_POD_HQ2) ###### 

# First create a separate Df
combinedSubAnnot <- data.frame(cellID <- colnames(seurat_merged),
                               mainAnnot <- as.character(seurat_merged$main_annot_HQ),
                               subcluster_annot_V2 <- as.character(seurat_merged$main_annot_HQ)) # first, we give subcluster_annot_V2 the same annotation as the main annot, which will be changed later
colnames(combinedSubAnnot) <- c("cellID", "mainAnnot",  "subcluster_annot_V2")
head(combinedSubAnnot)

levels(seurat_merged$main_annot_HQ)
# [1] "CNT"     "EC"      "PT"      "unknown" "IC_A"    "FIB"     "MES"     "VSMC_P"  "PEC"     "T_NKT"   "IC_B"    "B_MYEL"  "POD"     "TAL"     "DCT"     "DTL"    

## Mapping PODOCYTE sub-annotation onto UMAP
seurat_POD_HQ2 <- readRDS(file = "path/seurat_POD_HQ2_harmony_annot.Rds") # clean/HQ subcluster
seurat_POD_HQ2 <- SetIdent(seurat_POD_HQ2, value = "POD_HQ2_annot")
DimPlot(seurat_POD_HQ2, label = T)
head(seurat_POD_HQ2@meta.data)
levels(seurat_POD_HQ2$POD_HQ2_annot)

## Store annotation in ClusterAnnotations
clusterAnnotations <- Idents(seurat_POD_HQ2)
head(clusterAnnotations)

## check if cell IDs are in the same order:
podoCells <- combinedSubAnnot[combinedSubAnnot$cellID %in% colnames(seurat_POD_HQ2),] #subset Df with CellIDs that are present in seurat_POD_HQ2-object (i.e. the high quality POD)
head(podoCells)
all(names(clusterAnnotations) == podoCells$cellID) # TRUE

combinedSubAnnot$subcluster_annot_V2[combinedSubAnnot$mainAnnot == "POD"] <- "lowQ" #first call every podocyte lowQ
combinedSubAnnot$subcluster_annot_V2[combinedSubAnnot$cellID %in% podoCells$cellID] <- as.character(clusterAnnotations) #now add the HQ PODs

## Mapping sub-annotations to previous UMAP
seurat_merged$subcluster_annot_V2 <- plyr::mapvalues(x = as.character(colnames(seurat_merged)), from = colnames(seurat_merged), to = combinedSubAnnot$subcluster_annot_V2)
seurat_merged$subcluster_annot_V2 <- factor(seurat_merged$subcluster_annot_V2)
levels(seurat_merged$subcluster_annot_V2)
head(seurat_merged@meta.data)

## Plotting
seurat_merged <- SetIdent(seurat_merged, value = "subcluster_annot_V2")

pdf(paste0(projectFolder, "seurat_merged_DimPlot_POD_HQ2_mapped.pdf"), width = 11, height = 10)
DimPlot(seurat_merged, group.by = "subcluster_annot_V2", label = T)
dev.off()


###### 3.3: Remapping clean/HQ PEC cluster (seurat_PEC_HQ) ###### 

# Identical remapping procedure done for PECs


###### 3.4: Remapping clean/HQ EC cluster (seurat_EC_HQ) ###### 

# Identical remapping procedure done for ECs


###### 3.5: Remapping clean/HQ PT cluster (seurat_PT_HQ2) ###### 

# Identical remapping procedure done for PTs


###### 3.6: Remapping clean/HQ LOH_DCT_CNT cluster (seurat_LOH_DCT_CNT_HQ) ###### 

# Identical remapping procedure done for LOH, DCT and CNT


###### 3.7: Remapping clean/HQ IC_A cluster (seurat_IC_A_HQ) ######

# Identical remapping procedure done for IC A


###### 3.8: Remapping clean/HQ IC_B cluster (seurat_IC_B_HQ) ###### 

# Identical remapping procedure done for IC B


###### 3.9: Remapping clean/HQ STROMAL and IMMUNE cluster (seurat_stromal_HQ) ###### 

# Identical remapping procedure done for stromal cells and immune cells

## Plotting the result
seurat_merged <- SetIdent(seurat_merged, value = "subcluster_annot_V2")

pdf(paste0(projectFolder, "seurat_merged_DimPlot_all_mapped.pdf"), width = 11, height = 10)
DimPlot(seurat_merged, group.by = "subcluster_annot_V2", label = T)
dev.off()

###### 3.10: Remove UNKNOWN cluster (seurat_stromal_HQ) ###### 
seurat_merged <- SetIdent(seurat_merged, value = "subcluster_annot_V2")
seurat_merged <- subset(seurat_merged, ident = "unknown", invert = T)
seurat_merged$subcluster_annot_V2 <- factor(seurat_merged$subcluster_annot_V2, levels = c("POD1","POD2","POD3",
                                                                                          "PEC1","PEC2","PEC3","PEC4",
                                                                                          "EC_GC","EC_AEA_DVR","EC_PTC_AVR", "EC_others","EC_cyc",
                                                                                          "PT_S1-S2","PT_S2-S3","cPT","aPT","DTL",
                                                                                          "TAL","DCT","CNT",
                                                                                          "IC_A","IC_B",
                                                                                          "MES","REN","VSMC_P","FIB","mFIB","FAP-FIB","aFIB",
                                                                                          "B","PL","T","NKT","MDC","MAC_M2","ncMON", "cDC", "pDC",
                                                                                          "lowQ"))
seurat_merged <- SetIdent(seurat_merged, value = "subcluster_annot_V2")

pdf(paste0(projectFolder, "seurat_merged_DimPlot_all_mapped_refactorized.pdf"), width = 11, height = 10)
DimPlot(seurat_merged, group.by = "subcluster_annot_V2", label = T)
dev.off()


###### 3.11: Make annotation 'main_annot_HQ2', simplifying subcluster_annot_V2 ########
seurat_merged <- SetIdent(seurat_merged, value = "subcluster_annot_V2")
table(seurat_merged$subcluster_annot_V2)

seurat_merged <- RenameIdents(seurat_merged,'POD1'='POD') 
seurat_merged <- RenameIdents(seurat_merged,'POD2'='POD') 
seurat_merged <- RenameIdents(seurat_merged,'POD3'='POD') 
seurat_merged <- RenameIdents(seurat_merged,'PEC1'='PEC') 
seurat_merged <- RenameIdents(seurat_merged,'PEC2'='PEC') 
seurat_merged <- RenameIdents(seurat_merged,'PEC3'='PEC') 
seurat_merged <- RenameIdents(seurat_merged,'PEC4'='PEC') 
seurat_merged <- RenameIdents(seurat_merged,'EC_GC'='EC') 
seurat_merged <- RenameIdents(seurat_merged,'EC_AEA_DVR'='EC') 
seurat_merged <- RenameIdents(seurat_merged,'EC_PTC_AVR'='EC') 
seurat_merged <- RenameIdents(seurat_merged,'EC_others'='EC') 
seurat_merged <- RenameIdents(seurat_merged,'EC_cyc'='EC') 
seurat_merged <- RenameIdents(seurat_merged,'PT_S1-S2'='PT') 
seurat_merged <- RenameIdents(seurat_merged,'PT_S2-S3'='PT') 
seurat_merged <- RenameIdents(seurat_merged,'cPT'='PT') 
seurat_merged <- RenameIdents(seurat_merged,'aPT'='PT') 
seurat_merged <- RenameIdents(seurat_merged,'DTL'='DTL') 
seurat_merged <- RenameIdents(seurat_merged,'mFIB'='FIB') 
seurat_merged <- RenameIdents(seurat_merged,'FAP-FIB'='FIB') 
seurat_merged <- RenameIdents(seurat_merged,'aFIB'='FIB') 
seurat_merged <- RenameIdents(seurat_merged,'PL'='B') 
seurat_merged <- RenameIdents(seurat_merged,'T'='T_NKT') 
seurat_merged <- RenameIdents(seurat_merged,'NKT'='T_NKT') 
seurat_merged <- RenameIdents(seurat_merged,'MDC'='MYEL') 
seurat_merged <- RenameIdents(seurat_merged,'MAC_M2'='MYEL') 
seurat_merged <- RenameIdents(seurat_merged,'ncMON'='MYEL') 
seurat_merged <- RenameIdents(seurat_merged,'cDC'='MYEL') 
seurat_merged <- RenameIdents(seurat_merged,'pDC'='MYEL') 

DimPlot(seurat_merged, label=T)

## Cluster IDs in 'main_annot_HQ2'-column
seurat_merged$main_annot_HQ2 <- Idents(seurat_merged)
seurat_merged$main_annot_HQ2 <- factor(seurat_merged$main_annot_HQ2, levels = c("POD", "PEC", 
                                                                                "EC",
                                                                                "PT", "DTL",
                                                                                "TAL", "DCT", "CNT", "IC_A", "IC_B", 
                                                                                "MES", "REN", "VSMC_P", "FIB", 
                                                                                "T_NKT","B", "MYEL", 
                                                                                "lowQ"))
seurat_merged <- SetIdent(seurat_merged, value = "main_annot_HQ2")

pdf(paste0(projectFolder, "seurat_merged_main_annot_HQ2_dimPlot_umap.pdf"), width = 10, height = 8)
DimPlot(seurat_merged, group.by = "main_annot_HQ2", reduction = "umap", label = TRUE, pt.size = 0.7, raster.dpi = c(1024,1024)) + ggtitle("main_annot_HQ2")
dev.off()

library("RColorBrewer")
my_col <- c("#CAB2D6", "#FB9A99", 
            "#E31A1C",
            "#1F78B4", "#68AFDE",
            "#A6CEE3", "#B2DF8A", "#33A02C", "#FF7F00", "#FF7F00",
            "#509AA6", "#509AA6", "#509AA6", "#509AA6",
            "#FFFB74", "#FFFB74", "#FFFB74",
            "grey")

pdf(paste0(projectFolder, "seurat_merged_main_annot_HQ2_dimPlot_umap_coloured_including_lowQ.pdf"), width = 10, height = 10)
DimPlot(seurat_merged, group.by = "main_annot_HQ2", label = T, raster = F, repel = T, cols = my_col) + NoLegend()
dev.off()


###### 3.12: Deleting lowQ cells and creating object 'seurat_merged_HQ2' where lowQ is excluded ###### 
seurat_merged <- SetIdent(seurat_merged, value = "main_annot_HQ2")
DimPlot(seurat_merged, label = T)
seurat_merged_HQ2 <- subset(seurat_merged, ident = "lowQ", invert = T)

pdf(paste0(projectFolder, "seurat_merged_HQ2_main_annot_HQ2_dimPlot_umap.pdf"), width = 10, height = 10)
DimPlot(seurat_merged_HQ2, group.by = "main_annot_HQ2", label = T, raster = F, repel = T, cols = my_col) + NoLegend()
dev.off() 

## This concludes part 3


##### 4. PART 4 - Re-run Seurat pipeline on seurat_merged_HQ2 (no additional removal of cells) #####

projectFolder <- "path/HQ3/"

seurat_merged_HQ2
# An object of class Seurat 
# 36780 features across 120751 samples within 1 assay 
# Active assay: RNA (36780 features, 2000 variable features)
# 3 dimensional reductions calculated: pca, umap, harmony

table(seurat_merged_HQ2$main_annot_HQ2) # main annotation
table(seurat_merged_HQ2$subcluster_annot_V2) # subcluster annotation

###### 4.1: Re-running RunUMAP ######

numPCs <- 30 # identical PCs used as in part 2
seurat_merged_HQ2 <- RunUMAP(seurat_merged_HQ2, reduction = "harmony", dims = 1:numPCs)

###### 4.2: Checking previous annotation ###### 

# subcluster_annot_V2
seurat_merged_HQ2 <- SetIdent(seurat_merged_HQ2, value = "subcluster_annot_V2")

pdf(paste0(projectFolder, "seurat_merged_HQ2_DimPlot_subcluster_annot_V2.pdf"), width = 11, height = 10)
DimPlot(seurat_merged_HQ2, group.by = "subcluster_annot_V2", label = T, pt.size = 0.7, raster.dpi = c(1024,1024))
dev.off()

# main_annot_HQ2
seurat_merged_HQ2$main_annot_HQ2 <- factor(seurat_merged_HQ2$main_annot_HQ2, levels = c("POD", "PEC",
                                                                                "PT", "DTL",
                                                                                "TAL", "DCT", "CNT", "IC_A", "IC_B",
                                                                                "EC",
                                                                                "MES", "REN", "VSMC_P", "FIB", 
                                                                                "T_NKT","B", "MYEL"))
levels(seurat_merged_HQ2$main_annot_HQ2)
seurat_merged_HQ2 <- SetIdent(seurat_merged_HQ2, value = "main_annot_HQ2")

pdf(paste0(projectFolder, "seurat_merged_HQ2_DimPlot_main_annot_HQ2.pdf"), width = 11, height = 10)
DimPlot(seurat_merged_HQ2, group.by = "main_annot_HQ2", label = T, pt.size = 0.7, raster.dpi = c(1024,1024))
dev.off()

###### 4.3: Checking kidney cell atlas annotation (SingleR) ###### 

# Checking kidney atlas projection l1
levels(seurat_merged_HQ2$kidneyAtlas_prediction_l1)

pdf(paste0(projectFolder, "seurat_merged_HQ2_harmony_dimplot_KidneyAtlas_prediction_l1.pdf"), width = 10, height = 10)
DimPlot(seurat_merged_HQ2, group.by = "kidneyAtlas_prediction_l1", label = T, raster=FALSE)
dev.off()

# Checking kidney atlas projection l2
levels(seurat_merged_HQ2$kidneyAtlas_prediction_l2)

pdf(paste0(projectFolder, "seurat_merged_HQ2_harmony_dimplot_KidneyAtlas_prediction_l2.pdf"), width = 14, height = 10)
DimPlot(seurat_merged_HQ2, group.by = "kidneyAtlas_prediction_l2", label = T, raster=FALSE)
dev.off()

###### 4.4: Save object after RunUMAP was re-ran ###### 

seurat_merged_HQ2 <- SetIdent(seurat_merged_HQ2, value = "main_annot_HQ2")
saveRDS(seurat_merged_HQ2, file = paste0(projectFolder, "seurat_merged_HQ2_after_ReRunUMAP.Rds")) # final seurat object at main annotation level
# seurat_merged_HQ2 <- readRDS(file = "path/HQ3/seurat_merged_HQ2_after_ReRunUMAP.Rds")


###### 4.5: QC: Re-run QC on CLUSTERS (output included in manuscript) ###### 
seurat_merged_HQ2 <- SetIdent(seurat_merged_HQ2, value = "main_annot_HQ2")

# Re-run per QC per cluster
ChosenRes <- "main_annot_HQ2" 
levels(seurat_merged_HQ2$main_annot_HQ2)
seurat_merged_HQ2 <- SetIdent(seurat_merged_HQ2, value = "main_annot_HQ2")
DimPlot(seurat_merged_HQ2, label=T)

# QC per cluster - VlnPlot
pdf(paste0(projectFolder, "seurat_merged_HQ2_QC_perCluster_VlnPlot_", ChosenRes, ".pdf"), width = 15, height = 7)
VlnPlot(seurat_merged_HQ2, features = c("nCount_RNA"),
        pt.size = 0, 
        cols = rep("#1F78B4", times = 17)) &
  geom_boxplot(width=0.1, fill="white", outlier.size = 0.5) & NoLegend() &
  stat_summary(fun.y = median, geom='point', size = 2, colour = "black")

VlnPlot(seurat_merged_HQ2, features = c("nFeature_RNA"), 
        pt.size = 0, 
        cols = rep("#1F78B4", times = 17)) &
  geom_boxplot(width=0.1, fill="white", outlier.size = 0.5) & NoLegend() &
  stat_summary(fun.y = median, geom='point', size = 2, colour = "black")

VlnPlot(seurat_merged_HQ2, features = c("percent.mt"), 
        pt.size = 0, 
        cols = rep("#1F78B4", times = 17)) &
  geom_boxplot(width=0.1, fill="white", outlier.size = 0.5) & NoLegend() &
  stat_summary(fun.y = median, geom='point', size = 2, colour = "black")

dev.off()


## Number of nuclei per cluster per sample
table1 <- table(seurat_merged_HQ2$main_annot_HQ2, seurat_merged_HQ2$orig.ident)
total <- colSums(table1)
table1_with_total <- rbind(table1, Total = total)
write.table(table1_with_total, 
            file = paste0(projectFolder, "seurat_merged_HQ2_cells_perCluster_table.csv"),
            quote = FALSE, sep = ";", col.names = TRUE, row.names = TRUE)


# QC per cluster - table - without MALAT1 and without mean values
seuratDF <- data.frame(Median_nCounts = round(tapply(seurat_merged_HQ2$nCount_RNA, seurat_merged_HQ2$main_annot_HQ2, median)),
                       Median_nFeatures = round(tapply(seurat_merged_HQ2$nFeature_RNA, seurat_merged_HQ2$main_annot_HQ2, median)),
                       Median_percent.mt = round(tapply(seurat_merged_HQ2$percent.mt, seurat_merged_HQ2$main_annot_HQ2, median), 2),
                       nCells = table(seurat_merged_HQ2$main_annot_HQ2))
rownames(seuratDF) <- seuratDF$nCells.Var1
seuratDF$nCells.Var1 <- NULL

colnames(seuratDF) <- c("Median\nnCounts", "Median\nnFeatures", "Median\npercent.mt", "Nuclei")

pdf(paste0(projectFolder, "seurat_merged_HQ2_QC_perCluster_table_", ChosenRes, "_numberOfCellsCountsFeaturesMito_withoutMALAT1.pdf"), height=10, width=12)
p <- tableGrob(seuratDF)
grid.arrange(p)
dev.off()

colnames(seuratDF) <- c("Median nCounts", "Median nFeatures", "Median percent.mt", "Nuclei")
seuratDF

write.table(seuratDF, file = paste0(projectFolder, "seurat_merged_HQ2_QC_perCluster_table_", ChosenRes, "_numberOfCellsCountsFeaturesMito_withoutMALAT1.csv"),
            quote = F, sep = ";", col.names = T, row.names = T)


###### 4.6: QC: Re-run QC on SAMPLES (output included in manuscript) ###### 

# QC Per sample - VlnPlot
levels(seurat_merged_HQ2$orig.ident)
seurat_merged_HQ2 <- SetIdent(seurat_merged_HQ2, value = "orig.ident")

pdf(paste0(projectFolder, "seurat_merged_HQ2_QC_perSample_VlnPlot.pdf"), width = 14, height = 5)
VlnPlot(seurat_merged_HQ2, features = c("nCount_RNA"),
        pt.size = 0, 
        cols = rep("#1F78B4", times = 25)) &
  geom_boxplot(width=0.1, fill="white", outlier.size = 0.5) & 
  NoLegend() &
  stat_summary(fun.y = median, geom='point', size = 2, colour = "black")

VlnPlot(seurat_merged_HQ2, features = c("nFeature_RNA"), 
        pt.size = 0, 
        cols = rep("#1F78B4", times = 25)) &
  geom_boxplot(width=0.1, fill="white", outlier.size = 0.5) & 
  NoLegend() &
  stat_summary(fun.y = median, geom='point', size = 2, colour = "black")

VlnPlot(seurat_merged_HQ2, features = c("percent.mt"), 
        pt.size = 0, 
        cols = rep("#1F78B4", times = 25)) &
  geom_boxplot(width=0.1, fill="white", outlier.size = 0.5) & 
  NoLegend() &
  stat_summary(fun.y = median, geom='point', size = 2, colour = "black")

dev.off()

## QC per sample
df <- data.frame(matrix(ncol = 5, nrow = length(seurat_merged_HQ2$orig.ident)))
x <- c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "Diagnosis")
colnames(df) <- x

df$orig.ident <- seurat_merged_HQ2$orig.ident
df$nCount_RNA <- seurat_merged_HQ2$nCount_RNA
df$nFeature_RNA <- seurat_merged_HQ2$nFeature_RNA
df$percent.mt <- seurat_merged_HQ2$percent.mt
df$percent.MALAT1 <- seurat_merged_HQ2$percent.MALAT1
df$Diagnosis <- seurat_merged_HQ2$Diagnosis

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

pdf(paste0(projectFolder, "seurat_merged_HQ2_QC_perSample_table_numberOfCellsCountsFeaturesMito_withoutMALAT1.pdf"), height=8, width=13)
p <- tableGrob(seuratDF)
grid.arrange(p)
dev.off()

colnames(seuratDF) <- c("orig.ident","Diagnosis","Median nCounts", "Median nFeatures", "Median percent.mt", "Nuclei")
seuratDF

write.table(seuratDF, file = paste0(projectFolder, "seurat_merged_HQ2_QC_perSample_table_numberOfCellsCountsFeaturesMito_withoutMALAT1.csv"),
            quote = F, sep = ";", col.names = T, row.names = T)


###### 4.7: QC: Re-run QC on DIAGNOSES (output included in manuscript) ###### 

# QC per Diagnosis - VlnPlot
seurat_merged_HQ2 <- SetIdent(seurat_merged_HQ2, value = "Diagnosis")

pdf(paste0(projectFolder, "seurat_merged_HQ2_QC_perDiagnosis_VlnPlot.pdf"), width = 14, height = 5)
VlnPlot(seurat_merged_HQ2, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.MALAT1"),
        pt.size = 0, ncol=5,
        cols = rep("#1F78B4", times = 4)) &
  geom_boxplot(width=0.1, fill="white", outlier.size = 0.5) & 
  NoLegend() &
  stat_summary(fun.y = median, geom='point', size = 2, colour = "black")

dev.off()

## QC per sample - table
seuratDF <- data.frame(Median_nCounts = round(tapply(seurat_merged_HQ2$nCount_RNA, seurat_merged_HQ2$Diagnosis, median)),
                       Median_nFeatures = round(tapply(seurat_merged_HQ2$nFeature_RNA, seurat_merged_HQ2$Diagnosis, median)),
                       Median_percent.mt = round(tapply(seurat_merged_HQ2$percent.mt, seurat_merged_HQ2$Diagnosis, median), 2),
                       nCells = table(seurat_merged_HQ2$Diagnosis))
rownames(seuratDF) <- seuratDF$nCells.Var1
seuratDF$nCells.Var1 <- NULL
colnames(seuratDF) <- c("Median\nnCounts", "Median\nnFeatures", "Median\npercent.mt", "Nuclei")

pdf(paste0(projectFolder, "seurat_merged_HQ2_QC_perDiagnosis_table_numberOfCellsCountsFeaturesMito_withoutMALAT1.pdf"), height=4, width=12)
p <- tableGrob(seuratDF)
grid.arrange(p)
dev.off()

colnames(seuratDF) <- c("Median nCounts", "Median nFeatures", "Median percent.mt", "Nuclei")
seuratDF

write.table(seuratDF, file = paste0(projectFolder, "seurat_merged_HQ2_QC_perDiagnosis_table_numberOfCellsCountsFeaturesMito_withoutMALAT1.csv"),
            quote = F, sep = ";", col.names = T, row.names = T)


###### 4.8: FindAllMarkers + plotting of DEGs (output included in manuscript) ###### 

## Finding all markers ##
ChosenRes <- "main_annot_HQ2"
seurat_merged_HQ2 <- SetIdent(seurat_merged_HQ2, value = "main_annot_HQ2")

allMarkers <- FindAllMarkers(seurat_merged_HQ2, only.pos = F, min.pct = 0.1, logfc.threshold = 0.25,
                             max.cells.per.ident = 10000)
allMarkers <- group_by(allMarkers, cluster)
write.table(allMarkers, file = paste0(projectFolder, "seurat_merged_HQ2_harmony_", ChosenRes, "_FindAllMarkers_minPCT01_logFC025.csv"),
            row.names = F, col.names = T, sep = "\t", quote = F)


## Heatmap of TOP markers per cluster ##
## Reading in marker genes - TOP10
allSignMarkers <- allMarkers[allMarkers$p_val_adj < 0.05,]
top_10 <- top_n(allSignMarkers, n = 10, wt = avg_log2FC)
top_10 <- top_10[order(top_10$cluster),]
top10Expr <- as.data.frame(AverageExpression(seurat_merged_HQ2, assays = "RNA", features = unique(top_10$gene), slot = "data")) #calculates average expression in clusters

top10Zscore <- t(scale(t(top10Expr)))
top10Zscore <- top10Zscore[!is.na(top10Zscore[,1]),]

# Heatmap of z-scores for top 10 marker genes
breaks <- c(seq(min(top10Zscore), 0, -min(top10Zscore)/20),
            seq(max(top10Zscore)/21, max(top10Zscore), max(top10Zscore)/20))
myColor <- rev(colorRampPalette(brewer.pal(n = 9, name = "RdBu"))(40))

pdf(paste0(projectFolder, "seurat_merged_HQ2_harmony_", ChosenRes, "_top10_heatmap_zscore.pdf"), width = 4, height = 10)
pheatmap(top10Zscore, cluster_rows = F, cluster_cols = F, color = myColor, fontsize_row = 5,
         breaks = breaks) #, scale = "column"
dev.off()

pheatmap(top10Zscore, cluster_rows = F, cluster_cols = F, color = myColor, fontsize_row = 5,
         breaks = breaks) #, scale = "column"
dev.off()

# dotplot of top10 genes per cluster
pdf(paste0(projectFolder, "seurat_merged_HQ2_harmony_", ChosenRes, "_top10_dotplot.pdf"), width = 10, height = 22)
DotPlot(seurat_merged_HQ2, features = rev(unique(top_10$gene))) + coord_flip()
dev.off()


###### 4.9: Clean UMAP and annotation DotPlot (output included in manuscript) ###### 
seurat_merged_HQ2 <- SetIdent(seurat_merged_HQ2, value = "main_annot_HQ2")

pdf(paste0(projectFolder, "seurat_merged_HQ2_DimPlot_main_annot_HQ2_myCol.pdf"), width = 11, height = 10)
DimPlot(seurat_merged_HQ2, group.by = "main_annot_HQ2", label = T, raster = F, repel = T)
dev.off()

## New annotation (clean)
greyscale <- c("#F0F0F0", "#000000") # setting colors for DotPlot

AnnotFeatures <- c("NPHS1", "NPHS2", # POD
                   "VCAM1", "CFH", # PEC
                   "CUBN", "SLC13A1", # PT
                   "SATB2", "UNC5D", # DTL
                   "SLC12A1", "UMOD", # TAL
                   "SLC12A3", # DCT 
                   "AQP2", "SCNN1G", # CNT 
                   "SLC4A1", "SLC26A4", # IC_A and IC_B
                   "PECAM1", "PLVAP", # EC
                   "PDGFRB", "GATA3", "REN", "ACTA2", "RGS5", "COL1A2", "C7", # stromal
                   "PTPRC", "IL7R", "CCR7", "CD3D", "GNLY", "MS4A1", "IGHM", "CD79B", "ITGAM", "CD163") # immune

pdf(paste0(projectFolder, "seurat_merged_HQ2_DotPlot_main_annot_HQ2.pdf"), width = 12, height = 7)
DotPlot(seurat_merged_HQ2, features = AnnotFeatures, cols = greyscale) + RotatedAxis()
dev.off()

table <- table(seurat_merged_HQ2$main_annot_HQ2)
table <- as.data.frame(table)
colnames(table) <- c("CellType","Nuclei")

total_nuclei <- sum(table$Nuclei) # Calculate the total number of nuclei
total_nuclei # [1] 120751

table$Proportion <- round((table$Nuclei / total_nuclei *100),1)

colnames(table) <- c("CellType","Nuclei (N)","Proportion (%)")

write.table(table, file = paste0(projectFolder, "seurat_merged_main_annot_HQ2_table_cells_and_proportions_perCluster.csv"),
            row.names = T, col.names = T, sep = "\t")

# This concludes part 4


##### 5. PART 5 - Cell cluster frequencies and abundancy analysis #####
projectFolder <- "path/Abundance_analysis/"

# Read in the dataset
seurat_merged_HQ2 <- readRDS(file = "path/HQ3/seurat_merged_HQ2_after_ReRunUMAP.Rds")

###### 5.1: Calculating proportion table (per sample) ###### 

# Proportion of clusters, per sample
library(data.table)
library(rstatix) # t_test and wlicox_test
library(ggpubr)

test <- seurat_merged_HQ2@meta.data %>% group_by(orig.ident) %>% summarize(sample_total = n()) #table with the sample ID and the number of nuclei per sample
seurat_merged_HQ2@meta.data <- left_join(seurat_merged_HQ2@meta.data, test, by = "orig.ident")
table(seurat_merged_HQ2$orig.ident, seurat_merged_HQ2$sample_total)

Seurat_DF <- seurat_merged_HQ2@meta.data[which(!is.na(seurat_merged_HQ2$orig.ident)), ] # The first line subsets the metadata table by extracting the rows where the "orig.ident" column is not missing
Seurat_DF <- Seurat_DF[which(!is.na(Seurat_DF$main_annot_HQ2)), ] # The second line further subsets "Seurat_DF" by extracting the rows where the "seurat_merged_HQ2" column is not missing

table(Seurat_DF$orig.ident, Seurat_DF$main_annot_HQ2)
DF <- Seurat_DF
#DF$orig.ident <- as.factor(DF$orig.ident) # already factorized before

# Create a probability table of each celltype per sample (in essence the fraction of cells assigned to a given type per sample, the sum of which equals 1)
L <-tapply (DF$main_annot_HQ2, DF$orig.ident, table)
for(i in 1:length(L)) {
  L[[i]] <- prop.table(L[[i]])
}

i<- 1
K <-  as.data.frame((L[[i]])) %>%tibble::rownames_to_column(var ="Cell")
K$Pat <-rep(names(L)[i], dim(K)[1])
F <- K


for(i in 2:length(L)) {
  K <-  as.data.frame((L[[i]])) %>%tibble::rownames_to_column(var ="Cell")
  K$Pat <-rep(names(L)[i], dim(K)[1])
  F <-rbind(F,K)
}

F

#### Important that the names are characters, not factors
ANNO <- unique(as.data.frame(cbind(as.character(Seurat_DF$orig.ident), as.character(Seurat_DF$Diagnosis), Seurat_DF$sample_total)))

F$Diagnosis  <- plyr::mapvalues(x=F$Pat, from= ANNO$V1, to= as.character(ANNO$V2))
F$SampleTotal  <- plyr::mapvalues(x=F$Pat, from= ANNO$V1, to= as.character(ANNO$V3))

colnames(F) <- c("Cell", "CellType", "ratio", "SampleID", "Diagnosis", "Totals")
F$Diagnosis <- factor(F$Diagnosis, level=c("Primary FSGS", "Maladaptive FSGS", "PLA2R+ MN", "Healthy control"))
head(F)

F$Totals <- as.double(F$Totals)
F$Size <- round(F$ratio * F$Totals, 0)
#F$Diagnosis <- as.factor(F$Diagnosis) # already factorized before
F

F$Percentage <- round(F$ratio*100, 1)
head(F)

# Create a stacked bar chart
# Re-factor orig.ident according to diagnoses
levels(F$SampleID)
F$SampleID <- factor(F$SampleID, levels=rev(c("snrFSGS001", 
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
                                              "snrFSGS017",
                                              "snrNEPH029",
                                              "snrNEPH031",
                                              "snrNEPH033",
                                              "snrNEPH027",
                                              "snrNEPH028",
                                              "snrNEPH030",
                                              "snrNEPH032")))
levels(F$SampleID)

###### 5.2: Horizontal bar chart (per sample) ###### 

# Create the stacked bar chart with specified colors
barchart_plot <- ggplot(F, aes(x = ratio, y = SampleID, fill = CellType)) +
  geom_bar(stat = "identity") +
  labs(title = "Proportion of cell types per sample",
       x = "Proportion",
       y = "Sample ID") +
  theme_minimal()

pdf(paste0(projectFolder, "seurat_merged_HQ2_stackedbarplot_allCelltypes.pdf"), width = 7, height = 10)
barchart_plot
dev.off()

###### 5.3: Fraction plot - main_annot_HQ2 - all cell types (output included in manuscript) ###### 

#Fraction plots
fractions_plot <- ggplot(F,aes(x =CellType, y = Percentage, fill = Diagnosis , size = Size)) + 
  geom_boxplot(outlier.shape =NA, position=position_dodge(width=1), size = 0.2, colour = "black") +
  facet_grid(cols = vars(CellType), scales = "free", switch = "both") +
  theme_classic() +
  theme(plot.title = element_text(size = 15, colour = "black"),
        axis.text.x =element_text(angle =45,hjust = 1, size =15, colour = "black"),
        axis.title.x = element_text(size = 15, colour = "black"),
        axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),
        strip.text.x = element_blank()) +
  ggtitle("Main annotation", subtitle = "According to diagnoses") +
  geom_point(position =position_jitterdodge(), shape = 21, color = "black") + 
  scale_fill_manual(values = c("#1F7DC5","#76C6FB","#29A3FF", "grey")) +
  guides(fill = guide_legend(reverse = FALSE)) +
  scale_y_continuous(labels = scales::percent_format(scale = 1, suffix = "%"),
                     name = "% of total cells, per sample")

pdf(paste0(projectFolder, "seurat_merged_HQ2_fractionplot_allCelltypes.pdf"), width = 20, height = 10)
fractions_plot
dev.off()

###### 5.4: Statistical testing of abundance - main_annot_HQ2 - all cell types ###### 

# Let's perform Kruskal-Wallis test on these CellTypes: checking for statistically significant differences
library(FSA) # for post-hoc pairwise comparisons Dunn-test

head(F)
F$Percentage <- F$ratio*100 # unrounded values for statistical testing
head(F)

# POD
F_filtered <- F %>% filter(CellType == "POD")

kruskal_result <- kruskal.test(Percentage ~ Diagnosis, data = F_filtered)
print(kruskal_result)

# Kruskal-Wallis rank sum test
# 
# data:  Percentage by Diagnosis
# Kruskal-Wallis chi-squared = 6.1005, df = 3, p-value = 0.1068

# Check if the Kruskal-Wallis test is significant
if (kruskal_result$p.value < 0.05) {
  
  # Perform Dunn test as follow-up pairwise test
  dunn_result <- dunnTest(Percentage ~ Diagnosis, data = F_filtered, method = "bonferroni")
  print(dunn_result)
} else {
  # If the Kruskal-Wallis test is not significant, no follow-up pairwise test is needed
  print("Kruskal-Wallis test not significant.")
}

output <- list(
  statistic = kruskal_result$statistic,
  p.value = kruskal_result$p.value,
  method = kruskal_result$method,
  data = kruskal_result$data.name)

# Save the output of Kruskal Wallis and follow-up Dunn test to a CSV file
write.csv2(output, file = paste0(projectFolder, "seurat_merged_HQ2_abundance_POD_KruskalTest.csv"))

# The same is done for all other cell subclusters: PEC, PT, DTL, TAL, DCT, CNT, IC_A, IC_B, EC, MES, REN, VSMC_P, FIB, T_NKT, B, MYEL, 
# Additional fraction plots are made for immune cells only

###### 5.5: Subclustering of glomerular cells (for abundancy analysis and downstream MultiNicheNet anaysis) ###### 

# Read in the dataset
seurat_merged_HQ2 <- readRDS(file = "path/HQ3/seurat_merged_HQ2_after_ReRunUMAP.Rds")

# Subcluster glomerular cells
seurat_merged_HQ2 <- SetIdent(seurat_merged_HQ2, value = "subcluster_annot_V2")
seurat_merged_HQ2_glomerular <- subset(seurat_merged_HQ2, ident = c("POD1","POD2","POD3",
                                                                    "PEC1","PEC2","PEC3","PEC4",
                                                                    "EC_GC",
                                                                    "MES"))

seurat_merged_HQ2_glomerular <- RenameIdents(seurat_merged_HQ2_glomerular, "POD1" = "POD")
seurat_merged_HQ2_glomerular <- RenameIdents(seurat_merged_HQ2_glomerular, "POD2" = "POD")
seurat_merged_HQ2_glomerular <- RenameIdents(seurat_merged_HQ2_glomerular, "POD3" = "POD")
seurat_merged_HQ2_glomerular <- RenameIdents(seurat_merged_HQ2_glomerular, "PEC1" = "PEC")
seurat_merged_HQ2_glomerular <- RenameIdents(seurat_merged_HQ2_glomerular, "PEC2" = "PEC")
seurat_merged_HQ2_glomerular <- RenameIdents(seurat_merged_HQ2_glomerular, "PEC3" = "PEC")
seurat_merged_HQ2_glomerular <- RenameIdents(seurat_merged_HQ2_glomerular, "PEC4" = "PEC")
seurat_merged_HQ2_glomerular <- RenameIdents(seurat_merged_HQ2_glomerular, "EC_GC" = "EC_GC")
seurat_merged_HQ2_glomerular <- RenameIdents(seurat_merged_HQ2_glomerular, "MES" = "MES")

seurat_merged_HQ2_glomerular$glomerular_annot <- Idents(seurat_merged_HQ2_glomerular)
seurat_merged_HQ2_glomerular$glomerular_annot <- factor(seurat_merged_HQ2_glomerular$glomerular_annot, levels = c("POD", 
                                                                                                                  "PEC", 
                                                                                                                  "MES",
                                                                                                                  "EC_GC"))

seurat_merged_HQ2_glomerular <- SetIdent(seurat_merged_HQ2_glomerular, value = "glomerular_annot")

table(seurat_merged_HQ2$main_annot_HQ2)
# POD    PEC     PT    DTL    TAL    DCT    CNT   IC_A   IC_B     EC    MES    REN VSMC_P    FIB  T_NKT      B   MYEL 
# 2471   1574  34217   1256  26890   9125  18724   4215   1922  12485    839    267    826   2417   1536    908   1079 

table(seurat_merged_HQ2_glomerular$glomerular_annot)
# POD   PEC   MES EC_GC 
# 2471  1574   839  4280 

seurat_merged_HQ2_glomerular
# An object of class Seurat 
# 36780 features across 9164 samples within 1 assay 
# Active assay: RNA (36780 features, 2000 variable features)
# 3 dimensional reductions calculated: pca, umap, harmony

# Next, we make additional fraction plots and statistical testing of difference in cell type abundance across diagnoses in only glomerular cells

# This concludes part 5