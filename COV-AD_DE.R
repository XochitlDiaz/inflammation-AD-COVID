# AD - COVID DE 
# (aaaaa stop it with the accronyms)
# R 4.2.0
# 28 Nov 2023

library(Seurat)
library(ggplot2)
library(cowplot)
library(patchwork)
library(Azimuth)
library(gridExtra)
library(tidyr)
library(ggplot2)
library(reshape2)
library(SeuratDisk)
library(SeuratObject)
library(RColorBrewer)
library(VennDiagram)

setwd("/working/joint_projects/AD-COVID/inflammation_organoids/AD-COVtissue")

# Using data from two-way-integration.R, in which COVID labels where 
# transferred to the AD data set.

# COVID tissue data processing before transfer anchors
#setwd("/working/joint_projects/AD-COVID/inflammation_organoids")
#COVtissue_expmat = readRDS("./COVtissue_expmatV2.rds")
#COVtissue_expmat = ScaleData(COVtissue_expmat, verbose = F)
#COVtissue_expmat = RunPCA(COVtissue_expmat, npcs = 30, verbose = F)
#COVtissue_expmat <- RunUMAP(COVtissue_expmat, reduction = "pca", dims = 1:30, umap.method = "umap-learn", metric = "correlation",verbose = FALSE, return.model = T)

# AD tissue data and processing to add anchors
#ADtissue_expmat = readRDS("./ADtissue_GSE/filt_expmat.rds")
#Drop some seurat layers
#ADtissue_expmat= DietSeurat(ADtissue_expmat)
# Ran FindTransferAnchors()
# Ran TransferData()
# Ran AddMetaData()
# Ran MapQuery()
#saveRDS(ADtissue_expmat, "./two-way-integration/ADtissue_expmat_transfered.rds")

# From two-way_integration_p2.R,
# Data used was: 
#setwd("/working/joint_projects/AD-COVID/inflammation_organoids/two-way-integration")
#COVtissue_expmat = readRDS("./COVtissue_reductions.rds")
#COVorg_expmat= readRDS("./COVorg_expmat_transfered.rds")
#ADtissue_expmat = readRDS("./ADtissue_expmat_transfered.rds")

##############################################################################

# Load Data
COVtissue_expmat = readRDS("../two-way-integration/COVtissue_reductions.rds")
ADtissue_expmat = readRDS("../two-way-integration/ADtissue_expmat_transfered.rds")

merged_seurat <- merge(x = COVtissue_expmat, y = ADtissue_expmat, add.cell.ids = c("COVID_tissue_ds", "AD_tissue_ds"))

saveRDS(merged_seurat, "COV-ADv1.rds")
SaveH5Seurat(merged_seurat, filename = "COV-ADv1.h5Seurat")
Convert("COV-ADv1.h5Seurat", dest = "h5ad")

#https://satijalab.org/seurat/articles/integration_mapping

#need to find library to do this in older scripts
#seurat_write_h5(merged_seurat,"COV-ADv1.h5")
#https://rdrr.io/github/JiekaiLab/RIOH5/man/seurat_write_h5.html