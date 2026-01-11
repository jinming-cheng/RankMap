
# 2026.01.11
# R codes for generating vignette data
#
# 50 cells from each of the three cell types (Tumor, Basal, LP) were sampled 
# to reduce the example dataset size

# Raw data were downloaded from 10x Genomics website
# 10x Xenium human breast cancer data were used in this analysis:
# https://www.10xgenomics.com/products/xenium-in-situ/preview-dataset-human-breast

# Processed data with annotations were obtained from
# Cheng et al. 2025. BMC Bioinformatics. (PMID: 39833693)
# Cell type annotations can be downloaded from the GitHub repository: 
# https://github.com/jinming-cheng/cheng_annotation_bmc_bioinfo

###############################################
#                                             #
# #1. prepare single-cell data (seurat object)

filename_rds = "data_seurat_cell_type_Flex_Sample1.rds"
sc_ref_seu = readRDS(filename_rds)

# subset data with three cell types
sc_ref_seu@meta.data$cell_name <- colnames(sc_ref_seu)
df <- sc_ref_seu@meta.data[,c("cell_name","cell_type_label")]
colnames(df) <- c("cell_id","cell_type")
df <- subset(df,cell_type %in% c("Tumor", "Basal","LP") )

# select 50 cells for each cell type
set.seed(42)
selected_cells <- RankMap::SampleCellsByType(df, n_total_cells = 150)
seu <- sc_ref_seu[,selected_cells$cell_id]

# single-cell data for vignette
library(Seurat)
meta_seu <- seu@meta.data[,c("cell_name","cell_type_label")]
colnames(meta_seu) <- c("cell_id","cell_type")
meta_seu$cell_type <- droplevels(meta_seu$cell_type)
seu <- CreateSeuratObject(counts = GetAssayData(seu, layer = "counts"),
                          meta.data = meta_seu)
seu <- NormalizeData(seu)
saveRDS(seu,file = "seu_sc.rds", compress = "xz")

#                                             #
###############################################


###############################################
#                                             #
# #2. prepare spatial data (seurat object)

filename_rds = "xenium_obj_cell_type_Rep1_Xenium_v2.rds"
xenium_obj = readRDS(filename_rds)

xenium_obj@meta.data$cell_name <- colnames(xenium_obj)
df <- xenium_obj@meta.data[,c("cell_name","predicted_cell_type_SingleR")]
colnames(df) <- c("cell_id","cell_type")
df <- subset(df,cell_type %in% c("Tumor", "Basal","LP") )

# select 50 cells for each cell type
set.seed(42)
selected_cells <- RankMap::SampleCellsByType(df, n_total_cells = 150)
seu_sp <- xenium_obj[,selected_cells$cell_id]

# spatial data for vignette
library(Seurat)
meta_seu_sp <- seu_sp@meta.data[,c("cell_name","predicted_cell_type_SingleR")]
colnames(meta_seu_sp) <- c("cell_id","cell_type_SingleR")
meta_seu_sp$cell_type_SingleR <- droplevels(meta_seu_sp$cell_type_SingleR)
seu_sp <- CreateSeuratObject(counts = GetAssayData(seu_sp, layer = "counts"),
                             meta.data = meta_seu_sp)
seu_sp <- NormalizeData(seu_sp)
saveRDS(seu_sp,file = "seu_xen.rds", compress = "xz")

#                                             #
###############################################


