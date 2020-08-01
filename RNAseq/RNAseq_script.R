# load packages
library(pheatmap)
library(RColorBrewer)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(DOSE)
library(clusterProfiler)
library(topGO)
library(pathview)
library(ggplot2)
library(ggpubr)
library(Seurat)
XY_heatmap <- function (seurat_obj=seurat_obj, group=group,genes=genes,all_num=all_num,assay_sel=assay_sel,labels_rot=labels_rot,
  color=color,min_and_max_cut=num_cut,new_names=new_names,show_row_names=show_row_names,mark_gene=mark_gene,label_size=label_size){
  message("Processed data begain")
  ATAC <- GetAssayData(seurat_obj,slot="data",assay=assay_sel)
  ATAC_sel <- ATAC[genes,]
  ATAC_sel <- as.matrix(ATAC_sel)
  ATAC_sel_zscore <- t(apply(ATAC_sel, 1, function(x) (x-mean(x))/sd(x)))
  sel_cutoff <- min(abs(range(ATAC_sel_zscore)))
  if (is.null(min_and_max_cut)){
    ATAC_sel_zscore[ATAC_sel_zscore > sel_cutoff] <- sel_cutoff
    ATAC_sel_zscore[ATAC_sel_zscore < -sel_cutoff] <- -sel_cutoff
    } else {
      ATAC_sel_zscore[ATAC_sel_zscore > min_and_max_cut] <- min_and_max_cut
      ATAC_sel_zscore[ATAC_sel_zscore < -min_and_max_cut] <- -min_and_max_cut
    }
  meta_info <- seurat_obj@meta.data
  meta_info <- meta_info[order(meta_info[,group],decreasing=F),]
  annotation = data.frame(new_anno=meta_info[group],cell_id=rownames(meta_info))
  colnames(annotation) <- c("new_anno","cell_id")
  if (all_num == TRUE & !is.null(new_names)){
    annotation$new_anno <- paste(new_names,annotation$new_anno,sep="")
    aa <- as.data.frame(table(annotation$new_anno))
    aa <- aa[order(aa$Freq,decreasing=T),]
    aa$Var1 <- as.character(aa$Var1)
    annotation$new_anno <- factor(annotation$new_anno,levels=aa$Var1)
  }
  annotation <- annotation[order(annotation$new_anno),]
  annotation = data.frame(new_anno=annotation$new_anno,row.names=rownames(annotation))
  require(pheatmap)
  message("pheatmap printing start")
  require(ComplexHeatmap)
  require(BuenColors)
  require(scales)
  col1 <- jdb_palette("Darjeeling2")
  col2 <- jdb_palette("Darjeeling")
  col3 <- jdb_palette("Moonrise3")
  col_sel <- c(col1,col2,col3)
  col_sel <- hue_pal()(length(as.character(unique(annotation$new_anno))))
  col <- col_sel[1:length(as.character(unique(annotation$new_anno)))]
  names(col) <- as.character(unique(annotation$new_anno))
  top_anno <- HeatmapAnnotation(
  cluster = anno_block(gp = gpar(fill = col), # 设置填充色
  labels = as.character(unique(annotation$new_anno)),
  labels_gp = gpar(cex = label_size , col = "black"),labels_rot=labels_rot))
  if (is.null(mark_gene)){
  ph <- Heatmap(ATAC_sel_zscore[,rownames(annotation)],
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    show_row_names = show_row_names,
    top_annotation = top_anno,
   col = rev(color),
    column_split = annotation$new_anno,
    column_title_rot = 90)
    } else {
      both_gene <- intersect(rownames(ATAC_sel_zscore[,rownames(annotation)]),mark_gene)
      gene_pos <- which(rownames(ATAC_sel_zscore[,rownames(annotation)]) %in% both_gene)
      selected_gene <- rownames(ATAC_sel_zscore[,rownames(annotation)])[gene_pos]
      row_anno <-  rowAnnotation(mark_gene = anno_mark(at = gene_pos, labels = selected_gene))
      ph <- Heatmap(ATAC_sel_zscore[,rownames(annotation)],
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = show_row_names,
        top_annotation = top_anno,
        right_annotation = row_anno,
        col = rev(color),
        column_split = annotation$new_anno,
        column_title_rot = 90)
    }
    return(ph)
    print(ph)
}

# load data
RNAseq_PreB <- read.csv(file="1_RNAseq_PreB_shPhf23_vs_shRen_count_TPM_result.csv")
RNAseq_tumor <- read.csv(file="2_RNAseq_tumor_shPhf23_vs_shRen_count_TPM_result.csv")
RNAseq_Drug <- read.csv(file="3_RNAseq_HDAC_inhibitor_Drug_vs_DMSO_count_TPM_result.csv")
TCGA_DLBCL <- read.csv(file="4_TCGA_DLBCL_PHF23_RNA.csv")

*******************************************************************
*******************************************************************
*****************    GO plot figureS3 C-D    **********************
*****************    GO plot figureS3 C-D    **********************
*****************    GO plot figureS3 C-D    **********************
*******************************************************************
*******************************************************************
#figureS3C-D 
RNAseq_PreB_1 <- RNAseq_PreB[,c(9:14,16,20:22)]
RNAseq_PreB_1 <- subset(RNAseq_PreB_1,padj < 0.05)
RNAseq_PreB_up <- subset(RNAseq_PreB_1,log2FoldChange > 1)
RNAseq_PreB_down <- subset(RNAseq_PreB_1,log2FoldChange < -1)
RNAseq_PreB_up <- na.omit(RNAseq_PreB_up)
RNAseq_PreB_down <- na.omit(RNAseq_PreB_down)
dim(RNAseq_PreB_up)
dim(RNAseq_PreB_down)

ee  <-as.matrix(RNAseq_PreB_up$entrez)
dd <- as.vector(ee)
GO_RNAseq_PreB_up <- enrichGO(gene = dd, 
           OrgDb = org.Mm.eg.db,
      ont = "all", 
               pvalueCutoff = 0.01, 
                   pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05,
                   minGSSize = 10, 
                   maxGSSize = 500, 
                   readable = T, 
                   pool = FALSE)
GO_up_RNA <- as.data.frame(GO_RNAseq_PreB_up)
GO_up_RNA$log10_Pvalue <- -log(GO_up_RNA$pvalue,10)
GO_up_RNA_6 <- head(GO_up_RNA,6)
p1 <- ggbarplot(GO_up_RNA_6, 
  x = "Description", 
  y = "log10_Pvalue",
  color = "#FA6A63",            # Set bar border colors to white
  fill ="#FA6A63",
  sort.val = "asc",          # Sort the value in dscending order
  x.text.angle = 90,           # Rotate vertically x axis texts
  rotate = TRUE,
  title="PreB shPhf23_UP GO")
ggsave(p1,file="1_RNAseq_PreB_shPhf23_vs_shRen_up_GO.svg",width =10, height = 6,dpi=1080)

ee  <-as.matrix(RNAseq_PreB_down$entrez)
dd <- as.vector(ee)
GO_RNAseq_PreB_down <- enrichGO(gene = dd, 
           OrgDb = org.Mm.eg.db,
      ont = "all", 
               pvalueCutoff = 0.01, 
                   pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05,
                   minGSSize = 10, 
                   maxGSSize = 500, 
                   readable = T, 
                   pool = FALSE)
GO_down_RNA <- as.data.frame(GO_RNAseq_PreB_down)
GO_down_RNA$log10_Pvalue <- -log(GO_down_RNA$pvalue,10)
GO_down_RNA_6 <- head(GO_down_RNA,6)
p1 <- ggbarplot(GO_down_RNA_6, 
  x = "Description", 
  y = "log10_Pvalue",
  color = "#5B8FCF",            # Set bar border colors to white
  fill ="#5B8FCF",
  sort.val = "asc",          # Sort the value in dscending order
  x.text.angle = 90,           # Rotate vertically x axis texts
  rotate = TRUE,
  title="PreB shPhf23_DOWN GO")
ggsave(p1,file="1_RNAseq_PreB_shPhf23_vs_shRen_down_GO.svg",width =10, height = 6,dpi=1080)

*******************************************************************
*******************************************************************
*******    heatmap figure3D & figureS6 C-D     ********************
*******    heatmap figure3D & figureS6 C-D     ********************
*******    heatmap figure3D & figureS6 C-D     ********************
*******************************************************************
*******************************************************************

#figure3D  PreB_shPhf23_vs_shRen_heatmap
heatmap_gene <- rbind(RNAseq_PreB_up,RNAseq_PreB_down)
rownames(heatmap_gene) <- heatmap_gene$symbol
heatmap_gene <- heatmap_gene[order((heatmap_gene$log2FoldChange),decreasing = TRUE),]
heatmap_data <- heatmap_gene[,c(4:6,1:3)]
heatmap_resdeal_all_1_zscore <- t(apply(heatmap_data, 1, function(x) (x-mean(x))/sd(x)))
heatmap_resdeal_all_1_zscore <- na.omit(heatmap_resdeal_all_1_zscore)
preB_RNA <- CreateSeuratObject(counts = heatmap_resdeal_all_1_zscore, project = "preB")
preB_RNA$cluster <- c("PreB_Ren","PreB_Ren","PreB_Ren","PreB_shPhf23","PreB_shPhf23","PreB_shPhf23")
GO_Negative_apoptotic <- c("Bcl2a1a","Bcl2a1d","Cx3cl1")
diff <- c("Cebpb","Hk3","Ccr2","Itga1","Cd300a","Lyz2","Ccr5","Cebpa","Csf3r")
TSG <- c("Fpr2","Thbs1")
all_gene <- c(GO_Negative_apoptotic,diff,TSG)
mark_gene <- intersect(rownames(preB_RNA),all_gene)
pdf(file="1_RNAseq_PreB_shPhf23_vs_shRen_heatmap_gene.svg",width=5,height=7)
XY_heatmap(seurat_obj=preB_RNA,group="cluster",genes=rownames(preB_RNA),all_num=FALSE,new_names=NULL,labels_rot=90,
  assay_sel="RNA",color=colorRampPalette(c("firebrick3", "white", "navy"))(50),min_and_max_cut=2,show_row_names=FALSE,mark_gene=mark_gene,label_size=0)
dev.off()

#figureS6C  tumor_shPhf23_vs_shRen_heatmap
RNAseq_tumor_1 <- RNAseq_tumor[,c(9:14,16,20:21)]
RNAseq_tumor_1 <- subset(RNAseq_tumor_1,padj < 0.05)
RNAseq_tumor_1 <- subset(RNAseq_tumor_1,abs(log2FoldChange) >1)
RNAseq_tumor_1 <- na.omit(RNAseq_tumor_1)
RNAseq_tumor_up <- subset(RNAseq_tumor_1,log2FoldChange > 1)
RNAseq_tumor_down <- subset(RNAseq_tumor_1,log2FoldChange < -1)
dim(RNAseq_tumor_up)
dim(RNAseq_tumor_down)
RNAseq_Drug_1 <- RNAseq_Drug[,c(12:28)]
RNAseq_Drug_1 <- subset(RNAseq_Drug_1,padj < 0.05)
RNAseq_Drug_1 <- subset(RNAseq_Drug_1,abs(log2FoldChange) >1)
RNAseq_Drug_1 <- na.omit(RNAseq_Drug_1)
RNAseq_Drug_up <- subset(RNAseq_Drug_1,log2FoldChange > 1)
RNAseq_Drug_down <- subset(RNAseq_Drug_1,log2FoldChange < -1)
dim(RNAseq_Drug_up)
dim(RNAseq_Drug_down)
key_genes <- intersect(RNAseq_Drug_up$symbol,RNAseq_tumor_down$symbol)

rownames(RNAseq_tumor_1) <- RNAseq_tumor_1$symbol
RNAseq_tumor_heatmap <- RNAseq_tumor_1[order((RNAseq_tumor_1$log2FoldChange),decreasing = TRUE),]
RNAseq_tumor_heatmap <- RNAseq_tumor_heatmap[,c(1:6)]
heatmap_resdeal_all_1_zscore <- t(apply(heatmap_data, 1, function(x) (x-mean(x))/sd(x)))
heatmap_resdeal_all_1_zscore <- na.omit(heatmap_resdeal_all_1_zscore)
tumor_RNA_heatmap_se <- CreateSeuratObject(counts = heatmap_resdeal_all_1_zscore)
tumor_RNA_heatmap_se$group <- c("shRen","shRen","shRen","shPhf23","shPhf23","shPhf23")
tumor_RNA_heatmap_se$group <- factor(tumor_RNA_heatmap_se$group,levels=c("shRen","shPhf23"))
pdf(file="2_RNAseq_tumor_shPhf23_vs_shRen_heatmap_gene.svg",width=5,height=7) 
XY_heatmap(seurat_obj=tumor_RNA_heatmap_se,group="group",genes=rownames(tumor_RNA_heatmap_se),all_num=FALSE,new_names=NULL,labels_rot=90,
  assay_sel="RNA",color=colorRampPalette(brewer.pal(10, "RdBu"))(101),min_and_max_cut=2,show_row_names=FALSE,mark_gene=key_genes,label_size=0)
dev.off()

#figureS6D  HDAC_inhibitor_Drug_vs_DMSO_heatmap
rownames(RNAseq_Drug_1) <- RNAseq_Drug_1$symbol
RNAseq_Drug_heatmap <- RNAseq_Drug_1[order((RNAseq_Drug_1$log2FoldChange),decreasing = TRUE),]
RNAseq_Drug_heatmap <- RNAseq_Drug_heatmap[,c(1:9)]
heatmap_resdeal_all_1_zscore <- t(apply(RNAseq_Drug_heatmap, 1, function(x) (x-mean(x))/sd(x)))
heatmap_resdeal_all_1_zscore <- na.omit(heatmap_resdeal_all_1_zscore)
RNAseq_Drug_heatmap_se <- CreateSeuratObject(counts = heatmap_resdeal_all_1_zscore)
RNAseq_Drug_heatmap_se$group <- c("Chidamide","Chidamide","Chidamide","Entinostat","Entinostat","Entinostat","DMSO","DMSO","DMSO")
RNAseq_Drug_heatmap_se$group <- factor(RNAseq_Drug_heatmap_se$group,levels=c("DMSO","Entinostat","Chidamide"))
pdf(file="3_RNAseq_HDAC_inhibitor_Drug_vs_DMSO_heatmap_gene.svg",width=5,height=7) 
XY_heatmap(seurat_obj=RNAseq_Drug_heatmap_se,group="group",genes=rownames(RNAseq_Drug_heatmap_se),all_num=FALSE,new_names=NULL,labels_rot=90,
  assay_sel="RNA",color=colorRampPalette(brewer.pal(10, "RdBu"))(101),min_and_max_cut=2,show_row_names=FALSE,mark_gene=key_genes,label_size=0)
dev.off()

************************************************************************************
************************************************************************************
            boxplot TCGA DLBCL figure1B  
            boxplot TCGA DLBCL figure1B  
            boxplot TCGA DLBCL figure1B  
************************************************************************************
************************************************************************************
table(TCGA_DLBCL$CNV)
compare_means(Phf23_mRNA ~ CNV, data = TCGA_DLBCL)
range(TCGA_DLBCL$Phf23_mRNA)
TCGA_DLBCL$CNV <- factor(TCGA_DLBCL$CNV,levels=c("Diploid","Deletion"))
p1 <- ggboxplot(TCGA_DLBCL, x = "CNV", y = "Phf23_mRNA",
          color = "CNV", palette ="npg", add = "jitter",
         short.panel.labs = FALSE,size=1,ylim = c(-3, 3))+stat_compare_means(label = "p.format")
ggsave(p1, file="4_TCGA_DLBCL_PHF23_RNA_CNV_and_mRNA_RSEM.svg",width =7, height = 6,dpi=1080)

