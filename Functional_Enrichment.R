# Functional enrichment 

# 1 Dec 2023
# Xochitl Diaz

library(ggplot2)
library(gridExtra)
library(car)
library(msigdbr)
library(biomaRt)
library(clusterProfiler)

# Load data
#########################################################
# Set up for AD

setwd("/Users/xdiaz/Documents/Inflammation/AD_DE")
AD_all_markers= read.table("status_markers.csv", header = T)
oldcol = colnames(AD_all_markers)
colnames(AD_all_markers)= c(oldcol[2:6],oldcol[1],oldcol[7]) 

AD_celltypes = c("Excitatory", "Oligodendrocytes", "Microglia","Astrocytes",     
                 "Inhibitory", "OPCs", "Pericytes", "Endothelial")

# Keep only those cell types of interest
AD_celltypes = c(3,4,6)
# Oligodendrocyte progenitor cells
keep = which(AD_all_markers$cell_type %in% AD_celltypes)
AD_all_markers = AD_all_markers[keep,]

AD_all_markers$cell_type[AD_all_markers$cell_type==3] = "Microglia"
AD_all_markers$cell_type[AD_all_markers$cell_type==4] = "Astrocytes"
AD_all_markers$cell_type[AD_all_markers$cell_type==6] = "OPCs"
#########################################################
# Set up for COVID

setwd("/Users/xdiaz/Documents/Inflammation/COVID_DE")
COVID_all_markers= read.table("status_markers.csv", header = T)

colnames(COVID_all_markers)= c(oldcol[2:6],oldcol[1],oldcol[7]) 

COV_celltypes= c("Ast1", "Opc", "Mic")

keep = which(COVID_all_markers$cell_type %in% COV_celltypes)
COVID_all_markers = COVID_all_markers[keep,]

COVID_all_markers$cell_type[COVID_all_markers$cell_type=="Ast1"] = "Astrocytes"
COVID_all_markers$cell_type[COVID_all_markers$cell_type=="Opc"] = "OPCs"
COVID_all_markers$cell_type[COVID_all_markers$cell_type=="Mic"] = "Microglia"

########################################################

AD_filt_markers = subset(AD_all_markers, p_val_adj < 0.05 & (avg_log2FC < -0.6 | avg_log2FC > 0.6))
table(AD_filt_markers$cell_type)
#. Astrocytes  Microglia   OPCs 
#  53          111         16 
COV_filt_markers = subset(COVID_all_markers, p_val_adj < 0.05 & (avg_log2FC < -0.6 | avg_log2FC > 0.6))
table(COV_filt_markers$cell_type)
# Astrocytes  Microglia   OPCs 
# 1090        237         131


# Plotting overlapped Genes
olap = which(AD_all_markers$marker %in% COVID_all_markers$marker)
overlapping_ADtab = AD_all_markers[olap,]

olap = which(COVID_all_markers$marker %in% AD_all_markers$marker)
overlapping_COVIDtab = COVID_all_markers[olap,]

table(overlapping_ADtab$cell_type)
# Astrocytes  Microglia       OPCs
# 392         694             205

ast_AD= overlapping_ADtab[overlapping_ADtab$cell_type=="Astrocytes",]
ast_AD$condition = rep("AD",nrow(ast_AD))
ast_COV= overlapping_COVIDtab[overlapping_COVIDtab$cell_type=="Astrocytes",]
ast_COV$condition = rep("COVID",nrow(ast_COV))

astT = rbind(ast_COV,ast_AD)
astT = subset(astT, p_val_adj < 0.05 & (avg_log2FC < -0.6 | avg_log2FC > 0.6))
table(astT$cell_type,astT$condition)
#            AD    COVID
# Astrocytes 48    117
ast_olap= astT$marker[duplicated(astT$marker)]
#"GFAP"       "COLEC12"    "AC002429.2" "ARHGEF3"    "MAOB"     "LRMDA"     
# "ARL17B"     "PDZRN4"     "PTGDS"


mic_AD= overlapping_ADtab[overlapping_ADtab$cell_type=="Microglia",]
mic_AD$condition = rep("AD",nrow(mic_AD))
mic_COV= overlapping_COVIDtab[overlapping_COVIDtab$cell_type=="Microglia",]
mic_COV$condition = rep("COVID",nrow(mic_COV))

micT = rbind(mic_COV,mic_AD)
micT = subset(micT, p_val_adj < 0.05 & (avg_log2FC < -0.6 | avg_log2FC > 0.6))
table(micT$cell_type,micT$condition)
#            AD    COVID
# Microglia  98    88
mic_olap= micT$marker[duplicated(micT$marker)]
#[1] "AC008691.1" "FMN1"       "FKBP5"      "CCDC26"     "ACSL1"      "MS4A6A"    
#[7] "SIPA1L1"    "NAV2"       "PTPRG"      "TANC2"      "TMEM163"    "SYNDIG1"   
#[13] "RASGEF1C"   "C22orf34"   "LRRK2"      "MIR181A1HG" "CD83"       "TBC1D14"   
#[19] "XACT"       "APOE"       "TXNIP"      "ZBTB16"     "NHSL1"      "IPCEF1"    
#[25] "CPM"        "ARID5B"     "IQGAP2"     "TPRG1"      "CX3CR1"     "STARD13"   
#[31] "DIRC3"      "F13A1"      "LINC01141"  "NR4A2"      "P2RY12"     "ASTN1"     
#[37] "SPARCL1"    "PLCXD3"     "SGK1"       "ESRRG"      "SRGN"       "FCGBP"     
#[43] "MAP1B"      "ITGB5"      "MB21D2"     "PPP2R2B"    "ZNF331"     "PDK4"      
#[49] "CREM"

opc_AD= overlapping_ADtab[overlapping_ADtab$cell_type=="OPCs",]
opc_AD$condition = rep("AD",nrow(opc_AD))
opc_COV= overlapping_COVIDtab[overlapping_COVIDtab$cell_type=="OPCs",]
opc_COV$condition = rep("COVID",nrow(opc_COV))

opcT = rbind(opc_COV,opc_AD)
opcT = subset(opcT, p_val_adj < 0.05 & (avg_log2FC < -0.6 | avg_log2FC > 0.6))
table(opcT$cell_type,opcT$condition)
#            AD    COVID
#     OPCs   12    74
opc_olap= opcT$marker[duplicated(opcT$marker)]
#"FKBP5"  "ARL17B" "LIFR" 

write.csv(astT, "olapAD-COVID_astrocytes.csv")
write.csv(micT, "olapAD-COVID_microglia.csv")
write.csv(opcT, "olapAD-COVID_opc.csv")

#######################################################
Trestab = subset(resDA_Ctrl, p_val < 0.05 & (avg_log2FC < -0.6 | avg_log2FC > 0.6))
restab = Trestab[Trestab$cell_type=="Microglia",]
restab = Trestab[Trestab$cell_type=="Astrocytes",]
restab = Trestab[Trestab$cell_type=="OPCs",]
#write.csv(restab, paste0("res_filt_prot",nam_model,cellsite,".csv"))

#Table with genes that passed the test using adj.pval and logFC
Trestab2 = subset(resDA_Ctrl, p_val_adj < 0.05 & (avg_log2FC < -0.6 | avg_log2FC > 0.6))
restab2 = Trestab[Trestab$cell_type=="Microglia",]
col= "#0078D5"
restab2 = Trestab[Trestab$cell_type=="Astrocytes",]
col= "#FF7F7F"
restab2 = Trestab[Trestab$cell_type=="OPCs",]
col="#E88DFF"
#write.csv(restab2, paste0("res_filt_protPRO",nam_model,cellsite,".csv"))
cellsite = "AD n COVID - "  
cellsite = "AD" 
cellsite = "COVID" 

### ENSEMBL for GO set up ###
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)
###

geneNames <- restab2$marker
geneIDs <- bitr(geneNames, fromType = "SYMBOL", toType = "ENTREZID",OrgDb = "org.Hs.eg.db")


# Cellular Component (CC) sub-ontology  
nam_model ="Cellular Component GO"
goResults_cc <- enrichGO(gene = geneIDs$ENTREZID,
                         OrgDb = "org.Hs.eg.db",
                         ont = "CC",  
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.05,
                         readable = TRUE)

topTerms_GO_cc <- head(goResults_cc, n = 25)
print(topTerms_GO_cc)

png(filename = paste0("./GO_CC",nam_model,cellsite,".png"), width= 610, height = 640)
par(mar = c(6, 15, 1, 1.9))
# Plot the GO terms as a bar plot
bar_cc<-barplot(topTerms_GO_cc$Count, names.arg = topTerms_GO_cc$Description, horiz = TRUE, las = 1,
                sub =  paste0(nam_model, "-",cellsite),
                xlab = "Number of Genes", col = col)
dev.off()

GOterm_plot_CC = dotplot(goResults_cc, showCategory = 10)
png(filename = paste0("./GO_CCgr",nam_model,cellsite,".png"), width= 720, height = 570)
par(mar = c(5.1, 4.1, 4.1, 2.1))
GOterm_plot_CC + labs(title = "Enriched GO Cellular Component", subtitle = paste0(nam_model, "-",cellsite))
dev.off()


# Molecular Function (MF) sub-ontology
nam_model ="Molecular Function GO"
goResults_mf <- enrichGO(gene = geneIDs$ENTREZID,
                         OrgDb = "org.Hs.eg.db",
                         ont = "MF",  
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.05,
                         readable = TRUE)

topTerms_GO_mf <- head(goResults_mf, n = 25)
print(topTerms_GO_mf)

png(filename = paste0("./GO_MF",nam_model,cellsite,".png"), width= 610, height = 640)
par(mar = c(6, 20, 1, 1.9))
# Plot the GO terms as a bar plot
bar_mf<-barplot(topTerms_GO_mf$Count, names.arg = topTerms_GO_mf$Description, horiz = TRUE, las = 1,
                main = "Enriched GO Molecular Function",
                sub =  paste0(nam_model, "-",cellsite),
                xlab = "Number of Genes", col = col)
dev.off()

GOterm_plot_MF = dotplot(goResults_mf, showCategory = 10)
png(filename = paste0("./GO_MF_gene_ratio",nam_model,cellsite,".png"), width= 720, height = 570)
GOterm_plot_MF + labs(title = "Enriched GO Molecular Function", subtitle = paste0(nam_model, "-",cellsite))
dev.off()


# Biological Process (BP) sub-ontology
nam_model ="Molecular Function GO"
goResults_bp <- enrichGO(gene = geneIDs$ENTREZID,
                         OrgDb = "org.Hs.eg.db",
                         ont = "BP",  # 
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.05,
                         readable = TRUE)

topTerms_GO_bp <- head(goResults_bp, n = 25)
print(topTerms_GO_bp)

png(filename = paste0("./GO_BP",nam_model,cellsite,".png"), width= 610, height = 640)
par(mar = c(6, 15, 1, 1.9))
# Plot the GO terms as a bar plot
bar_bp<-barplot(topTerms_GO_bp$Count, names.arg = topTerms_GO_bp$Description, horiz = TRUE, las = 1,
                main = "Enriched GO Biological Process",
                sub =  paste0(nam_model, "-",cellsite),
                xlab = "Number of Genes", col = col)
dev.off()

GOterm_plot_BP = dotplot(goResults_bp, showCategory = 10)
png(filename = paste0("./GO_BP_gene_ratio",nam_model,cellsite,".png"), width= 720, height = 570)
GOterm_plot_BP + labs(title = "Enriched GO Biological Process", subtitle = paste0(nam_model, "-",cellsite))
dev.off()




