#Find cell markers against all other cell types
#6 Sep 2023

library(Seurat)
#library(hdf5)
#library(SeuratData)
library(SeuratDisk)

celltype= c("Excitatory", "Oligodendrocytes", "Microglia","Astrocytes",     
            "Inhibitory", "OPCs", "Pericytes", "Endothelial")

####################################
###############Setup################
#Select Cell to find markers
#cell= celltype[1]

#Select name of the previously processed Seurat Object file
matfile= "expression_matrix.h5.h5seurat"
####################################

#Read Seurat object
expmat= LoadH5Seurat(matfile)
#explore the file without loading
#expmat=Connect(matfile)

i=0
for (cell in celltype) {
i = i+1
ThisCell_markers = FindMarkers(expmat, ident.1=cell, ident.2=NULL)
ThisCell_markers$marker=row.names(ThisCell_markers)
ThisCell_markers$cell_type= rep(i,nrow(ThisCell_markers))


#append markers to existing table
write.table(ThisCell_markers, file = 'cell_markers.csv', append = T, row.names = F,col.names = F)
}