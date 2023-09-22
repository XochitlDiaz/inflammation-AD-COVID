#Marker plots
#18 Sep 2023
#R.4.0.2
#must configure proxy before using this script
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(reshape2)
library(dplyr)
library(cowplot)
library(pheatmap)
library(Matrix)
library(edgeR)
library(MAST)
library(ggrepel)
library(gprofiler2)
library(forcats)


# library(SeuratData)
# library(SeuratObject)
# library(biomaRt)
# library(DoubletFinder)

####################################
  ### SET UP FILES AND OUTPUT ###
####################################
#Edit
comparison = "Status Markers"
matfile = "./GSE214979_expression_matrix.h5.h5seurat"
all_markers = read.table("./status_markers.csv",header = T)
outdir = "/working/lab_miguelr/xochitlY/inflammation_organoids/ADtissue_GSE"
trhlog2FC = 0.6
trhpval   = 0.05

#Don't Change
celltype= c("Excitatory", "Oligodendrocytes", "Microglia","Astrocytes",     
            "Inhibitory", "OPCs", "Pericytes", "Endothelial")
colnames(all_markers)=c("p_val","avg_log2FC","pct.1","pct2","p_val_adj","marker", "cell")
comp= gsub(" ","",comparison)
####################################



#Read Seurat object
totexpmat= LoadH5Seurat(matfile)


all_markers$de="Not DE"
all_markers$de[all_markers$avg_log2FC >  trhlog2FC & all_markers$p_val < trhpval] = "Upregulated"
all_markers$de[all_markers$avg_log2FC < -trhlog2FC & all_markers$p_val < trhpval] = "Downregulated"



# #QC check sex from the samples
# # initialize connection to mart, may take some time if the sites are
# # unresponsive.
# mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
# # fetch chromosome info plus some other annotations
# genes.table <- try(biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name",
#                                                  "description", "gene_biotype", "chromosome_name", "start_position"), mart = mart,
#                                   useCache = F))
# chrY.gene = genes.table$external_gene_name[genes.table$chromosome_name == "Y"]
# 
# data.filt$pct_chrY = colSums(data.filt@assays$RNA@counts[chrY.gene, ])/colSums(data.filt@assays$RNA@counts)
# FeatureScatter(data.filt, feature1 = "XIST", feature2 = "pct_chrY")
# VlnPlot(data.filt, features = c("XIST", "pct_chrY"))

cell_i=0
for (ccell in celltype) {
cell_i= cell_i + 1

cell_markers=subset(all_markers, cell==cell_i)
expmat=subset(totexpmat, predicted.id==ccell)
Idents(expmat)="Status"
#Ploting QC ond UMAP distribution
png(file= paste0(outdir,"/",comp,"_",ccell, "UMAP.png"), width=800, height = 800)
plot_grid(ncol=2, DimPlot(expmat,label=T,group.by = "predicted.id")+NoAxes(),
          DimPlot(expmat,label=T,group.by = "Status")+NoAxes())
dev.off()


#Plotting differentialy expressed genes
top25=cell_markers[order(cell_markers$p_val_adj, decreasing = T),]
top25= subset(top25, p_val_adj!=0)
top25= top_n(top25,-25,p_val_adj)
png(file=paste0(outdir,"/",comp,"_",ccell, "TopDEG.png"), width=660, height=850)
par(mar=c(5.1,8,4.1,2.1))
barplot(sort(setNames(top25$avg_log2FC, top25$marker), F),
        horiz = T, las = 1, main = paste0(comparison," - ",ccell), border = "white", yaxs = "i",
        xlab = "Log2FC",col = "#4c6673")
abline(v = c(0, 0.25), lty = c(1, 2))
abline(v = c(0,-0.25), lty = c(1, 2))
dev.off()

#Pick top 25 up and down DE genes
upreg_markers = head(cell_markers[order(-cell_markers$avg_log2FC), ],25)
downreg_markers = head(cell_markers[order(cell_markers$avg_log2FC), ],25)

markers= c(upreg_markers$marker,downreg_markers$marker)



#Heat Map
png(file = paste0(outdir,"/",comp,"_",ccell, "Heatmap.png"), width = 1000, height = 800) 
  DoHeatmap(subset(expmat,downsample=100),features = markers,size=3) +
    labs(title = "Heatmap",
         subtitle = ccell,
         caption = paste0("Data Source: ",comparison))
dev.off()



#DotPlot
upreg_markers = head(cell_markers[order(-cell_markers$avg_log2FC), ],10)
downreg_markers = head(cell_markers[order(cell_markers$avg_log2FC), ],10)
markers= c(upreg_markers$marker,downreg_markers$marker)

png(file = paste0(outdir,"/",comp,"_",ccell,"DotPlot.png"), width = 1000, height = 800) 
  DotPlot(expmat,features=markers)+RotatedAxis() +
  labs(title = "DotPlot",
       subtitle = ccell,
       caption = paste0("Data Source: ",comparison))
dev.off()



#Violin Plot
upreg_markers = head(top25[order(-top25$avg_log2FC), ],6)
png(file = paste0(outdir,"/",comp,"_",ccell,"upregViolin.png"), width = 1000, height = 800) 
  VlnPlot(expmat,features = upreg_markers$marker,split.by="Status") +
   labs(title = "Upreg Violin Plot",
        subtitle = ccell,
        caption = paste0("Data Source: ",comparison))
dev.off()
dreg_markers = head(top25[order(top25$avg_log2FC), ],6)
png(file = paste0(outdir,"/",comp,"_",ccell,"dregViolin.png"), width = 1000, height = 800) 
 VlnPlot(expmat,features = dreg_markers$marker,split.by="Status") +
  labs(title = "Downreg Violin Plot",
       subtitle = ccell,
       caption = paste0("Data Source: ",comparison))
dev.off()

#Volcano Plot
png(file = paste0(outdir,"/",comp,"_",ccell,"Volcano.png"), width = 800, height = 800) 
ggplot(data=cell_markers, aes(x=avg_log2FC, y=-log10(p_val), col=de,label=marker)) + #<<
  geom_point() +
  #scale_color_manual(values=c("blue", "black", "red")) + #change dot colour #<<
  theme_minimal() +
  geom_text_repel() +
  xlim(-4.5,4) +
  ylim(0,320) + 
  labs(title = "Volcano plot",
       subtitle = ccell,
       caption = paste0("Data Source: ",comparison))
dev.off()



#Functional Analysis
genes_universe <- cell_markers$marker


# For up-regulated genes
# subset results for genes of interest
resSig <- subset(cell_markers, p_val_adj < trhpval & avg_log2FC > trhlog2FC) # subset
resSig <- resSig[ order(resSig$avg_log2FC, decreasing = TRUE), ]

# define gene lists
DEG <- resSig$marker

# enrichment analysis
gostres = gost(query = DEG, organism = "hsapiens", significant = F, 
               correction_method = "fdr", domain_scope = "custom_annotated", 
               custom_bg = genes_universe, ordered_query = TRUE)
go <- as.data.frame(gostres$result)

png(file = paste0(outdir,"/",comp,"_",ccell,"GOupreg.png"), width = 1800, height = 800)

# Reorder following the value of another column:
go %>%
  dplyr::slice(1:15) %>%
  select(p_value, term_name) %>%
  mutate(go = fct_reorder(term_name, p_value)) %>%
  ggplot( aes(x=go, y=p_value)) +
  geom_bar(stat="identity", fill="#4c6673", alpha=.6, width=.4) +
  coord_flip() +
  labs(x = "", y = "Adjusted P-value") +
  theme_bw(base_size = 20) +
  labs(title = "Enrichment Analysis \n Upregulated",
       subtitle = ccell,
       caption = paste0("Data Source: ",comparison))
dev.off()


# For down-regulated genes
# subset results for genes of interest
resSig <- subset(cell_markers, p_val_adj < trhpval & avg_log2FC < -trhlog2FC) # subset
resSig <- resSig[ order(resSig$avg_log2FC, decreasing = TRUE), ]

# define gene lists
DEG <- resSig$marker

# enrichment analysis
gostres = gost(query = DEG, organism = "hsapiens", significant = F, 
               correction_method = "fdr", domain_scope = "custom_annotated", 
               custom_bg = genes_universe, ordered_query = TRUE)
go <- as.data.frame(gostres$result)

png(file = paste0(outdir,"/",comp,"_",ccell,"GOdownreg.png"), width = 1800, height = 800)
# Reorder following the value of another column:
go %>%
  dplyr::slice(1:15) %>%
  select(p_value, term_name) %>%
  mutate(go = fct_reorder(term_name, p_value)) %>%
  ggplot( aes(x=go, y=p_value)) +
  geom_bar(stat="identity", fill="#4c6673", alpha=.6, width=.4) +
  coord_flip() +
  labs(x = "", y = "Adjusted P-value") +
  theme_bw(base_size = 20) +
  ggtitle("Enrichment Analysis \n Downregulated")
dev.off()

}

#immune.combined.sct.average<-AverageExpression(immune.combined.sct, return.seurat = T)


