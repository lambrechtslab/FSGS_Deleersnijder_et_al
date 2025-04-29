#### CELL-CELL INTERACTION ANALYSIS WITH MULTINICHENET ####

##### Script info #####

# In this script, we  used the dev-version of MultiNicheNet (i.e. V2 in development)
# Documentation from V1.0.1:  https://github.com/saeyslab/multinichenetr/blob/main/vignettes/basic_analysis_steps_MISC.md
# Documentation for the dev-branch: https://github.com/saeyslab/multinichenetr/blob/dev-branch/vignettes/basic_analysis_steps_MISC.Rmd

# Install the dev-branch:
devtools::install_github("saeyslab/multinichenetr", "dev-branch")

##### 0. PART 0 - Preparation of the analysis: load packages, NicheNet LR network & ligand-target matrix, single-cell expression data #####
library(nichenetr) 
library(multinichenetr) # note: install dev-branch
library(Seurat)
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)

## Setting wd and projectfolder:
setwd("path")
projectFolder <- "path/"

# Read in the dataset with the glomerular cells 
seurat_obj <- readRDS(file = "path/seurat_merged_HQ2_after_ReRunUMAP_glomerular_subset.Rds")

seurat_obj
# An object of class Seurat 
# 36780 features across 9164 samples within 1 assay 
# Active assay: RNA (36780 features, 2000 variable features)
# 3 dimensional reductions calculated: pca, umap, harmony

head(seurat_obj@meta.data)
DimPlot(seurat_obj, group.by = "glomerular_annot", label = T)

## Re-factor orig.ident according to diagnoses
levels(seurat_obj$orig.ident) # ok

## Re-factor diagnoses
levels(seurat_obj$Diagnosis) # ok

## Re-factor cells
levels(seurat_obj$glomerular_annot)
# "POD"   "PEC"   "MES"   "EC_GC" # ok

# Load NicheNet V2 networks and matrices from Zenodo

organism = "human"
if(organism == "human"){
  lr_network_all = readRDS(url("https://zenodo.org/record/10229222/files/lr_network_human_allInfo_30112033.rds")) %>% mutate(ligand = convert_alias_to_symbols(ligand, organism = organism), receptor = convert_alias_to_symbols(receptor, organism = organism))
  lr_network = lr_network_all  %>% mutate(ligand = make.names(ligand), receptor = make.names(receptor)) %>% distinct(ligand, receptor)
  ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
  colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% convert_alias_to_symbols(organism = organism) %>% make.names()
  rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% convert_alias_to_symbols(organism = organism) %>% make.names()
  lr_network = lr_network %>% filter(ligand %in% colnames(ligand_target_matrix))
  ligand_target_matrix = ligand_target_matrix[, lr_network$ligand %>% unique()]
} else if(organism == "mouse"){
  lr_network_all = readRDS(url("https://zenodo.org/record/10229222/files/lr_network_mouse_allInfo_30112033.rds")) %>% mutate(ligand = convert_alias_to_symbols(ligand, organism = organism), receptor = convert_alias_to_symbols(receptor, organism = organism))
  lr_network = lr_network_all  %>% mutate(ligand = make.names(ligand), receptor = make.names(receptor)) %>% distinct(ligand, receptor)
  ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
  colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% convert_alias_to_symbols(organism = organism) %>% make.names()
  rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% convert_alias_to_symbols(organism = organism) %>% make.names()
  lr_network = lr_network %>% filter(ligand %in% colnames(ligand_target_matrix))
  ligand_target_matrix = ligand_target_matrix[, lr_network$ligand %>% unique()]
}

# Make and prepare SingleCellExperiment object from Seurat object
sce = Seurat::as.SingleCellExperiment(seurat_obj, assay = "RNA")

# Because the NicheNet 2.0. networks are in the most recent version of the official gene symbols, 
# we will make sure that the gene symbols used in the expression data are also updated (= converted from their “aliases” to official gene symbols). 
# Afterwards, we will make them again syntactically valid.
sce =  alias_to_symbol_SCE(sce, "human") %>% makenames_SCE()

dim(sce)
# [1] 36780  9164
sce <- sce[rowSums(counts(sce) > 0) > 0, ] # this code removes rows from the SCE object where all the counts are zero, ensuring that only rows with at least one non-zero count value are retained.
dim(sce) 
# [1] 31474  9164 # rows (i.e. genes) with zero counts across all cells are filtered out

table(SummarizedExperiment::colData(sce)$glomerular_annot, SummarizedExperiment::colData(sce)$orig.ident) # cell types vs samples
table(SummarizedExperiment::colData(sce)$glomerular_annot, SummarizedExperiment::colData(sce)$Diagnosis) # cell types vs conditions

# Prepare settings of the MultiNicheNet cell-cell communication analysis

## Define the sampleID
sample_id = "orig.ident"
SummarizedExperiment::colData(sce)$orig.ident = SummarizedExperiment::colData(sce)$orig.ident %>% make.names()

## Define the condition/diagnosis
group_id = "Diagnosis"
SummarizedExperiment::colData(sce)$Diagnosis = SummarizedExperiment::colData(sce)$Diagnosis %>% make.names() # make it syntactically valid

## Define the annotation/cell ID (each sample should have sufficient cells per subtype)
celltype_id = "glomerular_annot"
SummarizedExperiment::colData(sce)$glomerular_annot = SummarizedExperiment::colData(sce)$glomerular_annot %>% make.names()

covariates = NA # we do not correct for additional covariates 
batches = NA # we do not correct for batch effects

# Sender and receiver cell types also need to be defined. 
# Both include all cell types in the dataset because we are interested in an All-vs-All analysis.

celltypes = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
senders_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
receivers_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()

##### 1. PART 1 - First steps of MultiNicheNet analysis #####

min_cells = 10 # recommended in the vignette

abundance_info = get_abundance_info(sce = sce, 
                                    sample_id = sample_id, 
                                    group_id = group_id, 
                                    celltype_id = celltype_id, 
                                    min_cells = min_cells, 
                                    senders_oi = senders_oi, 
                                    receivers_oi = receivers_oi, 
                                    batches = batches)

pdf(paste0(projectFolder, "MultiNicheNet_PART1_CellType_AbsoluteNumbers.pdf"), width = 10, height = 10)
abundance_info$abund_plot_sample
dev.off()

# For the DE analysis in the next step, only cell types will be considered if there are at least two samples per group with a sufficient number of cells.

saveRDS(abundance_info, file="MultiNicheNet_PART1_abundance_info.rds")


##### 2. PART 2 - Perform genome-wide DGE analysis of receiver and sender cell types to define DE genes between the conditions of interest #####

# Now we will go over to the multi-group, multi-sample differential expression (DE) analysis (also called ‘differential state’ analysis by the developers of Muscat).

# Define the contrasts and covariates of interest for the DE analysis.
table(sce$Diagnosis)
# Healthy.control Maladaptive.FSGS        PLA2R..MN     Primary.FSGS 
# 1413             2216             1497             4038 

# #Remove spaces from group_id
# SummarizedExperiment::colData(sce)[,group_id] <- gsub(" ",".",SummarizedExperiment::colData(sce)[,group_id])
# SummarizedExperiment::colData(sce)[,group_id] <- gsub("\\+",".",SummarizedExperiment::colData(sce)[,group_id]) 

contrasts_oi = c("'Primary.FSGS-(Healthy.control+Maladaptive.FSGS+PLA2R..MN)/3','Healthy.control-(Maladaptive.FSGS+PLA2R..MN+Primary.FSGS)/3','Maladaptive.FSGS-(Healthy.control+PLA2R..MN+Primary.FSGS)/3', 'PLA2R..MN-(Healthy.control+Maladaptive.FSGS+Primary.FSGS)/3'")
contrast_tbl = tibble(contrast =
                        c("Primary.FSGS-(Healthy.control+Maladaptive.FSGS+PLA2R..MN)/3","Healthy.control-(Maladaptive.FSGS+PLA2R..MN+Primary.FSGS)/3","Maladaptive.FSGS-(Healthy.control+PLA2R..MN+Primary.FSGS)/3","PLA2R..MN-(Healthy.control+Maladaptive.FSGS+Primary.FSGS)/3"),
                      group = c("Primary.FSGS","Healthy.control","Maladaptive.FSGS","PLA2R..MN"))

# First define expressed genes (this is new V2 functionality)
fraction_cutoff = 0.05
min_sample_prop = 0.50
frq_list = get_frac_exprs(sce = sce, 
                          sample_id = sample_id, 
                          celltype_id =  celltype_id, 
                          group_id = group_id, 
                          batches = batches, 
                          min_cells = min_cells, 
                          fraction_cutoff = fraction_cutoff, 
                          min_sample_prop = min_sample_prop)

# Perform the DE analysis for each cell type
DE_info = get_DE_info(sce = sce, 
                      sample_id = sample_id, 
                      group_id = group_id, 
                      celltype_id = celltype_id, 
                      batches = batches, 
                      covariates = covariates, 
                      contrasts_oi = contrasts_oi, 
                      min_cells = min_cells, 
                      expressed_df = frq_list$expressed_df # new V2 functionality
)

# Table with logFC and p-values for each gene-celltype-contrast:

DE_info$celltype_de$de_output_tidy %>% arrange(p_adj) %>% head()

celltype_de <- DE_info$celltype_de$de_output_tidy %>% arrange(p_adj)

write.csv2(celltype_de, file = paste0(projectFolder, "MultiNicheNet_PART2_DE_info.csv"),
           row.names = T, col.names = T, quote = F)

# It is always a good idea to check distribution of the p-values resulting from this DE expression analysis. 
# In order to trust the p-values, the p-value distributions should be uniform distributions, 
# with a peak allowed between 0 and 0.05 if there would be a clear biological effect in the data.

pdf(paste0(projectFolder, "MultiNicheNet_PART2_P-value_histograms.pdf"), width = 10, height = 8)
DE_info$hist_pvals
dev.off()

## Combine DE information for ligand-senders and receptors-receivers (similar to step1 - abundance_expression_info$sender_receiver_info)

sender_receiver_de = combine_sender_receiver_de(
  sender_de = celltype_de,
  receiver_de = celltype_de,
  senders_oi = senders_oi,
  receivers_oi = receivers_oi,
  lr_network = lr_network
)

sender_receiver_de %>% head(20)

write.csv2(sender_receiver_de, file = paste0(projectFolder, "MultiNicheNet_PART2_sender_receiver_de.csv"))

## Calculate (pseudobulk expression averages) = New V2 functionality!

abundance_expression_info = process_abundance_expression_info(sce = sce, 
                                                              sample_id = sample_id, 
                                                              group_id = group_id, 
                                                              celltype_id = celltype_id, 
                                                              min_cells = min_cells, 
                                                              senders_oi = senders_oi, 
                                                              receivers_oi = receivers_oi, 
                                                              lr_network = lr_network, 
                                                              batches = batches, 
                                                              frq_list = frq_list, 
                                                              abundance_info = abundance_info)

##### 3. PART 3 - Predict NicheNet ligand activities and NicheNet ligand-target links based on these differential expression results #####

# Define the parameters for the NicheNet ligand activity analysis

# Here we choose for a minimum logFC of 0.50, maximum p-value of 0.05, and minimum fraction of expression of 0.05.
logFC_threshold = 0.50
p_val_threshold = 0.05
#fraction_cutoff = 0.05  # this was already upstream defined

p_val_adj = FALSE # In case of more samples per group, and a sufficient high number of DE genes per group-celltype (> 20-50), the vignette recommends using the adjusted p-values. In this case, we used unadjusted P-values, as recommended.

# For the NicheNet ligand-target inference, 
# we also need to select which top n of the predicted target genes will be considered (here: top 250 targets per ligand).
top_n_target = 250 # default

cores_system = 18 # adapt according to machine
n.cores = min(cores_system, sender_receiver_de$receiver %>% unique() %>% length()) # use one core per receiver cell type

# Run the NicheNet ligand activity analysis, this takes some time
ligand_activities_targets_DEgenes = suppressMessages(suppressWarnings(get_ligand_activities_targets_DEgenes(
  receiver_de = celltype_de,
  receivers_oi = receivers_oi,
  ligand_target_matrix = ligand_target_matrix,
  logFC_threshold = logFC_threshold,
  p_val_threshold = p_val_threshold,
  p_val_adj = p_val_adj,
  top_n_target = top_n_target,
  n.cores = n.cores
)))

saveRDS(ligand_activities_targets_DEgenes, paste0(projectFolder, "MultiNicheNet_PART3_ligand_activities_targets_DEgenes.rds"))
#ligand_activities_targets_DEgenes <- readRDS(file = paste0(projectFolder, "MultiNicheNet_ligand_activities_targets_DEgenes.rds"))

ligand_activities_targets_DEgenes$de_genes_df %>% head(20)
write.csv2(ligand_activities_targets_DEgenes$de_genes_df, file = paste0(projectFolder, "MultiNicheNet_PART3_ligand_activities_targets_DEgenes_genes.csv"))

ligand_activities_targets_DEgenes$ligand_activities %>% head(20)
write.csv2(ligand_activities_targets_DEgenes$ligand_activities, file = paste0(projectFolder, "MultiNicheNet_PART3_ligand_activities_targets_DEgenes_ligandActivities.csv"))


##### 4. PART 4 - Use the information collected above to prioritize all sender-ligand—receiver-receptor pairs. #####

# In the 3 previous steps, 
# (1) we calculated expression, 
# (2) differential expression and 
# (3) NicheNet activity information. 
# Now we will combine these different types of information in one prioritization scheme.

# # In V1.0.1, you had to set the weights manually
# # In dev-branch, this is incorporated in the 'prioritization_tables'-function as 'scenario = regular'

# Make necessary grouping data frame
sender_receiver_tbl = sender_receiver_de %>% dplyr::distinct(sender, receiver)

metadata_combined = SummarizedExperiment::colData(sce) %>% tibble::as_tibble()

if(!is.na(batches)){
  grouping_tbl = metadata_combined[,c(sample_id, group_id, batches)] %>% tibble::as_tibble() %>% dplyr::distinct()
  colnames(grouping_tbl) = c("sample","group",batches)
} else {
  grouping_tbl = metadata_combined[,c(sample_id, group_id)] %>% tibble::as_tibble() %>% dplyr::distinct()
  colnames(grouping_tbl) = c("sample","group")
}

prioritization_tables = suppressMessages(generate_prioritization_tables(
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de = sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  contrast_tbl = contrast_tbl,
  sender_receiver_tbl = sender_receiver_tbl,
  grouping_tbl = grouping_tbl,
  scenario = "regular", # all prioritization criteria will be weighted equally
  fraction_cutoff = fraction_cutoff, 
  abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
  abundance_data_sender = abundance_expression_info$abundance_data_sender,
  ligand_activity_down = FALSE # by setting FALSE, we only use upregulatory ligand activities to prioritize, which is more intuitive and useful
))

# Check output tables
# First: group-based summary table
prioritization_tables$group_prioritization_tbl %>% head(20)

# Second: sample-based summary table: contains expression information of each LR pair per sample
prioritization_tables$sample_prioritization_tbl %>% head(20)

saveRDS(prioritization_tables, paste0(projectFolder, "MultiNicheNet_PART4_prioritization_tables.rds"))


##### 5. PART 5 - Add information on prior knowledge and expression correlation between LR and target expression. #####

# In multi-sample datasets, we have the opportunity to look whether expression of ligand-receptor across all samples is correlated with the expression of their by NicheNet predicted target genes. 
# This is what we will do with the following line of code:

lr_target_prior_cor = lr_target_prior_cor_inference(prioritization_tables$group_prioritization_tbl$receiver %>% unique(), 
                                                    abundance_expression_info, 
                                                    celltype_de, 
                                                    grouping_tbl, 
                                                    prioritization_tables, 
                                                    ligand_target_matrix, 
                                                    logFC_threshold = logFC_threshold, 
                                                    p_val_threshold = p_val_threshold, 
                                                    p_val_adj = p_val_adj)

# Save all the output of MultiNicheNet

multinichenet_output = list(
  celltype_info = abundance_expression_info$celltype_info,
  celltype_de = celltype_de,
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de =  sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  prioritization_tables = prioritization_tables,
  grouping_tbl = grouping_tbl,
  lr_target_prior_cor = lr_target_prior_cor) 

multinichenet_output = make_lite_output(multinichenet_output)

saveRDS(multinichenet_output, paste0(projectFolder, "Multinichenet_PART5_Multinichenet_all_output.rds"))

## Reload if necessary
multinichenet_output <- readRDS(file = "path/Multinichenet_PART5_Multinichenet_all_output.rds")

# Reload this as well, if you want to reload RDS file
contrasts_oi = c("'Primary.FSGS-(Healthy.control+Maladaptive.FSGS+PLA2R..MN)/3','Healthy.control-(Maladaptive.FSGS+PLA2R..MN+Primary.FSGS)/3','Maladaptive.FSGS-(Healthy.control+PLA2R..MN+Primary.FSGS)/3', 'PLA2R..MN-(Healthy.control+Maladaptive.FSGS+Primary.FSGS)/3'")
contrast_tbl = tibble(contrast =
                       c("Primary.FSGS-(Healthy.control+Maladaptive.FSGS+PLA2R..MN)/3","Healthy.control-(Maladaptive.FSGS+PLA2R..MN+Primary.FSGS)/3","Maladaptive.FSGS-(Healthy.control+PLA2R..MN+Primary.FSGS)/3","PLA2R..MN-(Healthy.control+Maladaptive.FSGS+Primary.FSGS)/3"),
                       group = c("Primary.FSGS","Healthy.control","Maladaptive.FSGS","PLA2R..MN"))
top_n_target = 250 # default


##### 6. PART 6 - Visualizing all results #####

# Visualization of the results of the cell-cell communication analysis

###### 6.1 - CIRCOS-plots ######

# We will look here at the top 50 predictions across all contrasts, senders, and receivers of interest.
prioritized_tbl_oi_all = get_top_n_lr_pairs(multinichenet_output$prioritization_tables, 50, rank_per_group = FALSE)

prioritized_tbl_oi = multinichenet_output$prioritization_tables$group_prioritization_tbl %>%
  filter(id %in% prioritized_tbl_oi_all$id) %>%
  distinct(id, sender, receiver, ligand, receptor, group) %>% left_join(prioritized_tbl_oi_all)
prioritized_tbl_oi$prioritization_score[is.na(prioritized_tbl_oi$prioritization_score)] = 0

senders_receivers = union(prioritized_tbl_oi$sender %>% unique(), prioritized_tbl_oi$receiver %>% unique()) %>% sort()

colors_sender = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)
colors_receiver = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)

# CircosPlot
pdf(paste0(projectFolder, "MultiNicheNet_PART6_CircosPlot_Top50.pdf"), width = 10, height = 10)
circos_list = make_circos_group_comparison(prioritized_tbl_oi, colors_sender, colors_receiver)
dev.off()


# Now you can also make a full circos plot for one group of interest, where we show the top30 per group
prioritized_tbl_oi_primary_30 = get_top_n_lr_pairs(multinichenet_output$prioritization_tables, 30, groups_oi = "Primary.FSGS")
prioritized_tbl_oi_maladaptive_30 = get_top_n_lr_pairs(multinichenet_output$prioritization_tables, 30, groups_oi = "Maladaptive.FSGS")
prioritized_tbl_oi_PLA2R_30 = get_top_n_lr_pairs(multinichenet_output$prioritization_tables, 30, groups_oi = "PLA2R..MN")
prioritized_tbl_oi_healthy_30 = get_top_n_lr_pairs(multinichenet_output$prioritization_tables, 30, groups_oi = "Healthy.control")

pdf(paste0(projectFolder, "MultiNicheNet_PART6_CircosPlot_Top30_Primary.pdf"), width = 10, height = 10)
make_circos_one_group(prioritized_tbl_oi_primary_30, colors_sender, colors_receiver)
dev.off()

pdf(paste0(projectFolder, "MultiNicheNet_PART6_CircosPlot_Top30_Maladaptive.pdf"), width = 10, height = 10)
make_circos_one_group(prioritized_tbl_oi_maladaptive_30, colors_sender, colors_receiver)
dev.off()

###### 6.2 - Scaled ligand-receptor and ligand activity ######

# Now do visualization of scaled ligand-receptor pseudobulk products and ligand activity - first default settings 

pdf(paste0(projectFolder, "MultiNicheNet_PART6_scaledLRpseudobulk_and_LA_default.pdf"), width = 20, height = 10)

group_oi = "Primary.FSGS"
prioritized_tbl_oi_M_50 = get_top_n_lr_pairs(multinichenet_output$prioritization_tables, 50, groups_oi = group_oi)
plot_oi = make_sample_lr_prod_activity_plots(multinichenet_output$prioritization_tables, prioritized_tbl_oi_M_50)
plot_oi

group_oi = "Maladaptive.FSGS"
prioritized_tbl_oi_M_50 = get_top_n_lr_pairs(multinichenet_output$prioritization_tables, 50, groups_oi = group_oi)
plot_oi = make_sample_lr_prod_activity_plots(multinichenet_output$prioritization_tables, prioritized_tbl_oi_M_50)
plot_oi


###### 6.3 - Visualization of expression-correlated target genes of ligand-receptor pairs ######

# Before, we had calculated the correlation in expression between ligand-receptor pairs and DE genes. 
# Now we will filter out correlated ligand-receptor –> target links that both show high expression correlation(spearman or activity > 0.50 in this example) and have some prior knowledge to support their link.

# With PEC as receiver

pdf(paste0(projectFolder, "MultiNicheNet_PART6_expressionCorrelated_targetGenes_PrimaryFSGS_PEC_receiver.pdf"), width = 20, height = 10)

group_oi = "Primary.FSGS"
receiver_oi = "PEC"
lr_target_prior_cor_filtered = multinichenet_output$lr_target_prior_cor %>% inner_join(multinichenet_output$ligand_activities_targets_DEgenes$ligand_activities %>% distinct(ligand, target, direction_regulation, contrast)) %>% inner_join(contrast_tbl) %>% filter(group == group_oi, receiver == receiver_oi)
lr_target_prior_cor_filtered_up = lr_target_prior_cor_filtered %>% filter(direction_regulation == "up") %>% filter( (rank_of_target < top_n_target) & (pearson > 0.50 | spearman > 0.50))
lr_target_prior_cor_filtered_down = lr_target_prior_cor_filtered %>% filter(direction_regulation == "down") %>% filter( (rank_of_target < top_n_target) & (pearson < -0.50 | spearman < -0.50)) # downregulation -- negative correlation
lr_target_prior_cor_filtered = bind_rows(lr_target_prior_cor_filtered_up, lr_target_prior_cor_filtered_down)

#Now we will visualize the top correlated target genes for the LR pairs that are also in the top 30 LR pairs discriminating the groups from each other:
prioritized_tbl_oi = get_top_n_lr_pairs(multinichenet_output$prioritization_tables, 30, groups_oi = group_oi, receivers_oi = receiver_oi)
lr_target_correlation_plot = make_lr_target_correlation_plot(multinichenet_output$prioritization_tables, prioritized_tbl_oi,  lr_target_prior_cor_filtered , multinichenet_output$grouping_tbl, multinichenet_output$celltype_info, receiver_oi,plot_legend = FALSE)
lr_target_correlation_plot$combined_plot
lr_target_correlation_plot$legends
dev.off()

pdf(paste0(projectFolder, "MultiNicheNet_PART6_expressionCorrelated_targetGenes_MaladaptiveFSGS_PEC_receiver.pdf"), width = 20, height = 10)

group_oi = "Maladaptive.FSGS"
receiver_oi = "PEC"
lr_target_prior_cor_filtered = multinichenet_output$lr_target_prior_cor %>% inner_join(multinichenet_output$ligand_activities_targets_DEgenes$ligand_activities %>% distinct(ligand, target, direction_regulation, contrast)) %>% inner_join(contrast_tbl) %>% filter(group == group_oi, receiver == receiver_oi)
lr_target_prior_cor_filtered_up = lr_target_prior_cor_filtered %>% filter(direction_regulation == "up") %>% filter( (rank_of_target < top_n_target) & (pearson > 0.50 | spearman > 0.50))
lr_target_prior_cor_filtered_down = lr_target_prior_cor_filtered %>% filter(direction_regulation == "down") %>% filter( (rank_of_target < top_n_target) & (pearson < -0.50 | spearman < -0.50)) # downregulation -- negative correlation
lr_target_prior_cor_filtered = bind_rows(lr_target_prior_cor_filtered_up, lr_target_prior_cor_filtered_down)

#Now we will visualize the top correlated target genes for the LR pairs that are also in the top 30 LR pairs discriminating the groups from each other:
prioritized_tbl_oi = get_top_n_lr_pairs(multinichenet_output$prioritization_tables, 30, groups_oi = group_oi, receivers_oi = receiver_oi)
lr_target_correlation_plot = make_lr_target_correlation_plot(multinichenet_output$prioritization_tables, prioritized_tbl_oi,  lr_target_prior_cor_filtered , multinichenet_output$grouping_tbl, multinichenet_output$celltype_info, receiver_oi,plot_legend = FALSE)
lr_target_correlation_plot$combined_plot
lr_target_correlation_plot$legends
dev.off()


###### 6.4 - Intercellular regulatory network systems view ######

prioritized_tbl_oi = get_top_n_lr_pairs(multinichenet_output$prioritization_tables, 250, rank_per_group = FALSE)

lr_target_prior_cor_filtered = multinichenet_output$prioritization_tables$group_prioritization_tbl$group %>% unique() %>% lapply(function(group_oi){
  lr_target_prior_cor_filtered = multinichenet_output$lr_target_prior_cor %>% inner_join(multinichenet_output$ligand_activities_targets_DEgenes$ligand_activities %>% distinct(ligand, target, direction_regulation, contrast)) %>% inner_join(contrast_tbl) %>% filter(group == group_oi)
  lr_target_prior_cor_filtered_up = lr_target_prior_cor_filtered %>% filter(direction_regulation == "up") %>% filter( (rank_of_target < top_n_target) & (pearson > 0.50 | spearman > 0.50))
  lr_target_prior_cor_filtered_down = lr_target_prior_cor_filtered %>% filter(direction_regulation == "down") %>% filter( (rank_of_target < top_n_target) & (pearson < -0.50 | spearman < -0.50))
  lr_target_prior_cor_filtered = bind_rows(lr_target_prior_cor_filtered_up, lr_target_prior_cor_filtered_down)
}) %>% bind_rows()

lr_target_df = lr_target_prior_cor_filtered  %>% distinct(group, sender, receiver, ligand, receptor, id, target, direction_regulation) 

network = infer_intercellular_regulatory_network(lr_target_df, prioritized_tbl_oi)
network$links %>% head()
network$nodes %>% head()

network_graph = visualize_network(network, colors_sender)

pdf(paste0(projectFolder, "MultiNicheNet_PART6_IntercellularRegulatoryNetwork_250.pdf"), width = 20, height = 10)
network_graph$plot
dev.off()

# An interesting follow-up question would be to find nodes that are at the top of the hierarchy in this network, and could thus explain the cascade of intercellular communication events we observe.
# Using PageRank on the opposite-directionality graph can help identify important nodes that act as sources in the hierarchical structure. 
# By running PageRank on the reversed graph (where edges point in the opposite direction), you can potentially pinpoint nodes that have significant influence or act as sources of information in your network.

## First for primary FSGS group

library(igraph)

# Create the normal graph (normal directionality)
normal_edges <- get.data.frame(network_graph$network_igraph) %>% filter(group == "Primary.FSGS") %>% .[, c("from", "to")]  
normal_graph <- graph_from_data_frame(normal_edges, directed = TRUE)

# Calculate PageRank on the normal graph
pagerank_scores <- page_rank(normal_graph)$vector %>% sort(decreasing = TRUE)
pr_score_tbl = tibble(node = names(pagerank_scores), score_normal = pagerank_scores) %>% arrange(-score_normal)

# Create the reversed graph (opposite directionality)
reversed_edges <- get.data.frame(network_graph$network_igraph) %>% filter(group == "Primary.FSGS") %>% .[, c("to", "from")]  # Swap 'from' and 'to' columns
reversed_graph <- graph_from_data_frame(reversed_edges, directed = TRUE)

# Calculate PageRank on the reversed graph
reversed_pagerank_scores <- page_rank(reversed_graph)$vector %>% sort(decreasing = TRUE)

pr_score_tbl = pr_score_tbl %>% inner_join(
  tibble(node = names(reversed_pagerank_scores), score_reverse = reversed_pagerank_scores) 
) %>% arrange(-score_reverse)
pr_score_tbl

# Let's zoom now in on the Primary.FSGS-group network to check this

network_subset = infer_intercellular_regulatory_network(lr_target_df %>% filter(group == "Primary.FSGS"), prioritized_tbl_oi %>% filter(group == "Primary.FSGS"))
network_subset$links %>% head()
network_subset$nodes %>% head()

network_subset_graph = visualize_network(network_subset, colors_sender)

pdf(paste0(projectFolder, "MultiNicheNet_PART6_IntercellularRegulatoryNetwork_PrimaryFSGS.pdf"), width = 8, height = 10)
network_subset_graph$plot
dev.off()

## Now for maladaptive FSGS group

# Create the normal graph (normal directionality)
normal_edges <- get.data.frame(network_graph$network_igraph) %>% filter(group == "Maladaptive.FSGS") %>% .[, c("from", "to")]  
normal_graph <- graph_from_data_frame(normal_edges, directed = TRUE)

# Calculate PageRank on the normal graph
pagerank_scores <- page_rank(normal_graph)$vector %>% sort(decreasing = TRUE)
pr_score_tbl = tibble(node = names(pagerank_scores), score_normal = pagerank_scores) %>% arrange(-score_normal)

# Create the reversed graph (opposite directionality)
reversed_edges <- get.data.frame(network_graph$network_igraph) %>% filter(group == "Maladaptive.FSGS") %>% .[, c("to", "from")]  # Swap 'from' and 'to' columns
reversed_graph <- graph_from_data_frame(reversed_edges, directed = TRUE)

# Calculate PageRank on the reversed graph
reversed_pagerank_scores <- page_rank(reversed_graph)$vector %>% sort(decreasing = TRUE)

pr_score_tbl = pr_score_tbl %>% inner_join(
  tibble(node = names(reversed_pagerank_scores), score_reverse = reversed_pagerank_scores) 
) %>% arrange(-score_reverse)
pr_score_tbl

# Let's zoom now in on the Maladaptive.FSGS-group network to check this

network_subset = infer_intercellular_regulatory_network(lr_target_df %>% filter(group == "Maladaptive.FSGS"), prioritized_tbl_oi %>% filter(group == "Maladaptive.FSGS"))
network_subset$links %>% head()
network_subset$nodes %>% head()

network_subset_graph = visualize_network(network_subset, colors_sender)

pdf(paste0(projectFolder, "MultiNicheNet_PART6_IntercellularRegulatoryNetwork_MaladaptiveFSGS.pdf"), width = 8, height = 10)
network_subset_graph$plot
dev.off()

# End of chapter 6 and script.