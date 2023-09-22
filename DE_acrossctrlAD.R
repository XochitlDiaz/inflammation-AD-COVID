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

#cols= c("marker","p_val", "avg_log2FC", "pct.1","pct.2","p_val_adj","cell_type")
#cell_markers= data.frame(matrix(nrow=0,ncol = length(cols)))
#colnames(cell_markers)=cols
#write.table(cell_markers, file = "status_markers.csv")
####################################

#Read Seurat object
expmat= LoadH5Seurat(matfile)
#explore the file without loading
#expmat=Connect(matfile)

i=0
for (cell in celltype) {
i = i+1
cellmat=subset(expmat,subset= predicted.id==cell)
Idents(cellmat)="Status"

ThisCell_markers = FindMarkers(cellmat, ident.1="AD", ident.2="Ctrl")
ThisCell_markers$marker=row.names(ThisCell_markers)
ThisCell_markers$cell_type= rep(i,nrow(ThisCell_markers))


#append markers to existing table
write.table(ThisCell_markers, file = 'status_markers.csv', append = T, row.names = F,col.names = F)
}