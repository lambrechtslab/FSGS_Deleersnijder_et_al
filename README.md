# FSGS_Deleersnijder_et_al

# Introduction

This repository contains code related to data processing and downstream analysis associated with the study "Single-Nucleus RNA-sequencing Identifies a Differential Profibrotic Response in Parietal Epithelial Cells in Primary versus Maladaptive FSGS”, under submission at Kidney International Reports (KI Reports).

# Data Availability

## Processed & filtered single-nucleus RNA-sequencing (snRNA-seq) data
Processed snRNA-seq data (i.e., annotated Seurat object with all cell types included) are available as Seurat-object entitled “seurat_merged_HQ2_after_ReRunUMAP.Rds” on EGA.

## Gene sets used for quality control
Genes related to cell stress can be found in the “stress.genes.csv”-file.
Genes related to cell cycle can be found in “regev_lab_cell_cycle_genes.txt”-file.

## Gene sets used for gene set enrichment analysis (GSEA)
We used fgsea package for GSEA. Used gene sets include Reactome, Hallmark, KEGG Medicus and Gene Ontology Biological Process (GOBP). Gene set .gmt-files are included in this repository.

## List and description of scripts

## SCRIPT 1: Creation of main Seurat object, QC, main clustering, annotation and abundancy analysis
This script starts from loading in Cellbender-corrected indivual samples, creating a merged object, and performing quality control (QC).
-	Part 1: We run Seurat pipeline for the first time with all low-quality cells still included. We perform first QC checks with first removal of low-quality cells.
-	Part 2: We re-run Seurat pipeline for the 2nd time after first removal of low-quality cells in part 1. We perform additional QC checks. We also use the SingleR-package to import the cell annotation of a publicly available snRNA-seq kidney dataset from Kidney Precision Medicine Project (KPMP), which aids in annotation of the Seurat object. We perform QC checks on every individual subcluster level (separate scripts) and identify additional low-quality cells on subcluster level.
-	Part 3: We re-map all "clean" (high-quality) subclusters back to the main annotation level. All low-quality cells that were identified on subcluster level (in part 2 and in separate analyses) are now also visualized on main annotation level and subsequently removed. This concludes QC of the main Seurat object, yielding 120,751 high-quality nuclei.
-	Part 4: We re-run Seurat pipeline for a final 3rd time to recreate a clean and reclustered UMAP; no additional cells are removed in this step. The final high-quality Seurat-object containing all celltypes is called "seurat_merged_HQ2_after_ReRunUMAP.Rds" (120,751 high-quality nuclei). The final main annotation can be found in meta.data "seurat_merged_HQ2$main_annot_HQ2". The final subcluster annotation can be found in meta.data "seurat_merged_HQ2$subcluster_annot_V2".
-	Part 5: We calculate cell cluster frequencies and abundancy analysis on the final main Seurat object.

## SCRIPT 2: Podocyte subclustering and analysis 
This script starts from loading in the podocyte (POD) subcluster, extracted from the main Seurat object in part 2 of script 1.
-	Part 1: We re-run Seurat pipeline for the first time on the podocyte subcluster with low-quality cells still included. We perform first QC checks with first removal of low-quality cells, based on visualization of QC parameters, doublet visualization, marker gene expression and KPMP annotation.
-	Part 2: We re-run Seurat pipeline for a 2nd time after first removal of low-quality cells in part 1. We perform additional QC checks and remove additional low-quality cells.
-	Part 3: We re-run Seurat pipeline for a final 3rd time after removal of low-quality cells in part 2. We re-annotate cell subclusters, but do not remove additional cells. Result of part 3 is the final podocyte-object: "seurat_POD_HQ2"-object with annotation "POD_HQ2_annot".
-	Part 4: We perform differential gene expression (DGE) analysis (using EdgeR package) of one diagnostic group vs. all other three groups. Output is visualized using heatmaps, dotplots and volcanoplots. Subsequently, the output from EdgeR DGE analysis is used for gene set enrichment analysis (GSEA) using the fgsea package.
-	Part 5: We perform DGE analysis (using EdgeR package) of one diagnostic group vs. only the two control groups.
-	Part 6: We perform DGE analysis (using EdgeR package) of primary FSGS vs. maladaptive FSGS.

## SCRIPT 3: Cell-cell interaction analysis with MultiNicheNet
In this script, we  used the dev-version of MultiNicheNet (i.e. V2 in development at time of writing the manuscript). 
Documentation for the dev-branch: https://github.com/saeyslab/multinichenetr/blob/dev-branch/vignettes/basic_analysis_steps_MISC.Rmd.
-	Part 0: Preparation of the analysis: load packages, NicheNet ligand-receptor (LR) network & ligand-target matrix, single-cell expression data.
-	Part 1: First steps of MultiNicheNet analysis.
-	Part 2: Perform genome-wide DGE analysis of receiver and sender cell types to define DE genes between the conditions of interest.
-	Part 3: Predict NicheNet ligand activities and NicheNet ligand-target links based on these differential expression results.
-	Part 4: Use the information collected above to prioritize all sender-ligand—receiver-receptor pairs.
-	Part 5: Add information on prior knowledge and expression correlation between LR and target expression.
-	Part 6: Visualizing all results.

# Contacts
All other relevant data and analysis are available from the authors upon request. For further enquires, please either raise an issue via GitHub or email Tom Venken (tom.venken[at]kuleuven.be), Diether Lambrechts (diether.lambrechts[at]kuleuven.be) or Amaryllis Van Craenenbroeck (amaryllis.vancraenenbroeck[at]kuleuven.be).

