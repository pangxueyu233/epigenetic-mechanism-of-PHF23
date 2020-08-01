# load packages
library(VennDiagram)
library(ggplot2)
library(scales)
library(dplyr)  
library(future)
library(future.apply)
library(ggpubr)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(clusterProfiler)

# load data
all_PHF23_all <- read.csv(file="1_Chipseq_PHF23_peaks_annotation.csv")
all_H3K4me3_all <- read.csv(file="1_Chipseq_H3K4me3_peaks_annotation.csv")
all_HDAC1_all <- read.csv(file="1_Chipseq_HDAC1_peaks_annotation.csv")
all_H3K27ac_all <- read.csv(file="1_Chipseq_H3K27ac_peaks_annotation.csv")
RNAseq_PreB <- read.csv(file="1_RNAseq_PreB_shPhf23_vs_shRen_count_TPM_result.csv")
PHF23_H3K4me3_PreB <- read.csv("1_PHF23_H3k4me3_PreB_basemean_binding_RNA.csv")
PHF23_H3K4me3_Baf3 <- read.csv("1_PHF23_H3k4me3_Baf3_basemean_binding_RNA.csv")
PHF23_HDAC1_TSS_H3K27ac <- read.csv(file="1_PHF23_HDAC1_TSS_H3K27ac_binding.csv")
PHF23_HDAC1_enhancer_H3K27ac <- read.csv(file="1_PHF23_HDAC1_enhancer_H3K27ac_binding.csv")
PHF23_HDAC1_TSS_PreB <- read.csv("1_PHF23_HDAC1_TSS_PreB_basemean_binding_RNA.csv")
PHF23_HDAC1_enhancer_PreB <- read.csv("1_PHF23_HDAC1_enhancer_PreB_basemean_binding_RNA.csv")
H3K4me3_H3K27ac_PreB <- read.csv("1_H3K4me3_H3K27ac_PreB_basemean_binding_RNA.csv")
H3K4me3_and_H3K27ac_PHF23_HDAC1_PreB <- read.csv("1_H3K4me3_and_H3K27ac_PHF23_HDAC1_PreB_basemean_binding_RNA.csv")
TP53_pathway <- read.csv("1_hallmark_tp53_pathway.csv")

*******************************************************************
*******************************************************************
******  Chipseq H3K4me3 PHF23 figure3B-C figureS3B   **************
******  Chipseq H3K4me3 PHF23 figure3B-C figureS3B   **************
******  Chipseq H3K4me3 PHF23 figure3B-C figureS3B   **************
*******************************************************************
*******************************************************************
#figure3B
my_plot <- venn.diagram(x=list(H3K4me3_all=na.omit(all_H3K4me3_all$SYMBOL),PHF23_all=na.omit(all_PHF23_all$SYMBOL)) , 
                          filename = NULL, 
                          fill=c('#FA6A63','#5B8FCF'), 
                          alpha=c(0.8,0.8),
                          euler.d = TRUE,
                          scaled = TRUE,
                          cat.fontfamily = "sans",
                          cat.cex = 1,
                          lwd = c(0.3, 0.3),
                          category.names=c("H3K4me3", "PHF23"), 
                          #cat.just=list(c(1,-40) , c(0.4,-40)),
                          print.mode = c("raw", "percent"),
                          fontfamily = "sans",
                          main="Overlap between H3K4me3 and PHF23",
                          main.fontfamily = "sans",
                          main.cex = 1,
                          cex = 0.9,
                          width = 10,
                          height = 10,
                          units = "in",
                          resolution = 400)
ggsave(my_plot, file="1_H3K4me3_PHF23_venn.svg", device = "svg")

#figure3C
co_binding_RNA <- subset(PHF23_H3K4me3_Baf3,cluster=="co_binding_RNA")[,c("BaFminus_01","BaFminus_02","BaFminus_03","BaFminus_04","BaFminus_05","cluster")]
only_PHF23_RNA <- subset(PHF23_H3K4me3_Baf3,cluster=="only_phf23_RNA")[,c("BaFminus_01","BaFminus_02","BaFminus_03","BaFminus_04","BaFminus_05","cluster")]
only_H3K4me3_RNA <- subset(PHF23_H3K4me3_Baf3,cluster=="only_H3K4me3_RNA")[,c("BaFminus_01","BaFminus_02","BaFminus_03","BaFminus_04","BaFminus_05","cluster")]
no_binding_RNA <- subset(PHF23_H3K4me3_Baf3,cluster=="no_binding_RNA")[,c("BaFminus_01","BaFminus_02","BaFminus_03","BaFminus_04","BaFminus_05","cluster")]
all_RNA <- rbind(co_binding_RNA,only_PHF23_RNA)
all_RNA <- rbind(all_RNA,only_H3K4me3_RNA)
all_RNA <- rbind(all_RNA,no_binding_RNA)
all_RNA$mean <- apply(all_RNA[,1:5],1,mean)
all_RNA$log_mean <- log(all_RNA$mean+0.1,2)
compare_means(log_mean ~ cluster, data = all_RNA)
all_RNA$cluster <- factor(all_RNA$cluster,levels=c("no_binding_RNA","only_H3K4me3_RNA","only_phf23_RNA","co_binding_RNA"))
p1 <- ggboxplot(all_RNA, x = "cluster", y = "log_mean",
          color = "cluster", palette = "jco",ylim = c(-4, 10),
         short.panel.labs = FALSE,size=1)
ggsave(p1, file="1_PHF23_H3K4me3_Baf3_basemean_binding_ggboxplot.svg",width =10, height = 10,dpi=1080)
median(subset(all_RNA,cluster=="co_binding_RNA")[,"log_mean"])
median(subset(all_RNA,cluster=="only_phf23_RNA")[,"log_mean"])
median(subset(all_RNA,cluster=="only_H3K4me3_RNA")[,"log_mean"])
median(subset(all_RNA,cluster=="no_binding_RNA")[,"log_mean"])

#figureS3B
co_binding_RNA <- subset(PHF23_H3K4me3_PreB,cluster=="co_binding_RNA")[,c("PreB_shRen_1_TPM","PreB_shRen_2_TPM","PreB_shRen_3_TPM","cluster")]
only_PHF23_RNA <- subset(PHF23_H3K4me3_PreB,cluster=="only_phf23_RNA")[,c("PreB_shRen_1_TPM","PreB_shRen_2_TPM","PreB_shRen_3_TPM","cluster")]
only_H3K4me3_RNA <- subset(PHF23_H3K4me3_PreB,cluster=="only_H3K4me3_RNA")[,c("PreB_shRen_1_TPM","PreB_shRen_2_TPM","PreB_shRen_3_TPM","cluster")]
no_binding_RNA <- subset(PHF23_H3K4me3_PreB,cluster=="no_binding_RNA")[,c("PreB_shRen_1_TPM","PreB_shRen_2_TPM","PreB_shRen_3_TPM","cluster")]
all_RNA <- rbind(co_binding_RNA,only_PHF23_RNA)
all_RNA <- rbind(all_RNA,only_H3K4me3_RNA)
all_RNA <- rbind(all_RNA,no_binding_RNA)
all_RNA$mean <- apply(all_RNA[,1:3],1,mean)
all_RNA$log_mean <- log(all_RNA$mean+0.1,2)
all_RNA$log_mean <- as.numeric(as.character(scale(all_RNA$log_mean)))
compare_means(log_mean ~ cluster, data = all_RNA)
all_RNA$cluster <- factor(all_RNA$cluster,levels=c("no_binding_RNA","only_H3K4me3_RNA","only_phf23_RNA","co_binding_RNA"))
p1 <- ggboxplot(all_RNA, x = "cluster", y = "log_mean",
          color = "cluster", palette = "jco",ylim = c(-1, 5),
         short.panel.labs = FALSE,size=1)
ggsave(p1, file="1_PHF23_H3k4me3_PreB_basemean_binding_ggboxplot.svg",width =10, height = 10,dpi=1080)
median(subset(all_RNA,cluster=="co_binding_RNA")[,"log_mean"])
median(subset(all_RNA,cluster=="only_phf23_RNA")[,"log_mean"])
median(subset(all_RNA,cluster=="only_H3K4me3_RNA")[,"log_mean"])
median(subset(all_RNA,cluster=="no_binding_RNA")[,"log_mean"])

*******************************************************************
*******************************************************************
******           Chipseq H3K27ac HDAC1 PHF23         **************
******           Chipseq H3K27ac HDAC1 PHF23         **************
******           Chipseq H3K27ac HDAC1 PHF23         **************
*******************************************************************
*******************************************************************
#figure4F
all_HDAC1_TSS <- subset(all_HDAC1_all,annotation_TSS=="TSS")
all_HDAC1_enhancer <- subset(all_HDAC1_all,annotation_TSS=="enhancer")
all_PHF23_TSS <- subset(all_PHF23_all,annotation_TSS=="TSS")
my_plot <- venn.diagram(x=list(PHF23_all=na.omit(all_PHF23_TSS$SYMBOL),HDAC1_TSS=na.omit(all_HDAC1_TSS$SYMBOL)) , 
                          filename = NULL, 
                          fill=c('#FA6A63','#5B8FCF'), 
                          alpha=c(0.8,0.8),
                          euler.d = TRUE,
                          scaled = TRUE,
                          cat.fontfamily = "sans",
                          cat.cex = 1,
                          lwd = c(0.3, 0.3),
                          category.names=c("PHF23", "HDAC1_TSS"), 
                          #cat.just=list(c(1,-40) , c(0.4,-40)),
                          print.mode = c("raw", "percent"),
                          fontfamily = "sans",
                          main="Overlap between PHF23 and HDAC1_TSS",
                          main.fontfamily = "sans",
                          main.cex = 1,
                          cex = 0.9,
                          width = 10,
                          height = 10,
                          units = "in",
                          resolution = 400)
ggsave(my_plot, file="2_PHF23_HDAC1_TSS_venn.svg", device = "svg")

my_plot <- venn.diagram(x=list(PHF23_all=na.omit(all_PHF23_TSS$SYMBOL),HDAC1_enhancer=na.omit(all_HDAC1_enhancer$SYMBOL)) , 
                          filename = NULL, 
                          fill=c('#FA6A63','#5B8FCF'), 
                          alpha=c(0.8,0.8),
                          euler.d = TRUE,
                          scaled = TRUE,
                          cat.fontfamily = "sans",
                          cat.cex = 1,
                          lwd = c(0.3, 0.3),
                          category.names=c("PHF23", "HDAC1_enhancer"), 
                          #cat.just=list(c(1,-40) , c(0.4,-40)),
                          print.mode = c("raw", "percent"),
                          fontfamily = "sans",
                          main="Overlap between PHF23 and HDAC1_enhancer",
                          main.fontfamily = "sans",
                          main.cex = 1,
                          cex = 0.9,
                          width = 10,
                          height = 10,
                          units = "in",
                          resolution = 400)
ggsave(my_plot, file="2_PHF23_HDAC1_enhancer_venn.svg", device = "svg")

#figure4G
compare_means(log_mean ~ cluster, data = PHF23_HDAC1_TSS_H3K27ac)
PHF23_HDAC1_TSS_H3K27ac$cluster <- factor(PHF23_HDAC1_TSS_H3K27ac$cluster,levels=c("no_binding","only_HDAC1","only_PHF23","co_binding"))
p1 <- ggboxplot(PHF23_HDAC1_TSS_H3K27ac, x = "cluster", y = "fold_change",
          color = "cluster", palette = "jco",
         short.panel.labs = FALSE,size=1,ylim = c(0, 10))
ggsave(p1, file="2_PHF23_HDAC1_TSS_with_H3k27ac_boxplot.svg",width =10, height = 10,dpi=1080)
median(subset(PHF23_HDAC1_TSS_H3K27ac,cluster=="co_binding")[,"fold_change"])
median(subset(PHF23_HDAC1_TSS_H3K27ac,cluster=="only_PHF23")[,"fold_change"])
median(subset(PHF23_HDAC1_TSS_H3K27ac,cluster=="only_HDAC1")[,"fold_change"])
median(subset(PHF23_HDAC1_TSS_H3K27ac,cluster=="no_binding")[,"fold_change"])

#figure4H
compare_means(log_mean ~ cluster, data = PHF23_HDAC1_enhancer_H3K27ac)
PHF23_HDAC1_enhancer_H3K27ac$cluster <- factor(PHF23_HDAC1_enhancer_H3K27ac$cluster,levels=c("no_binding","only_HDAC1","only_PHF23","co_binding"))
p1 <- ggboxplot(PHF23_HDAC1_enhancer_H3K27ac, x = "cluster", y = "fold_change",
          color = "cluster", palette = "jco",
         short.panel.labs = FALSE,size=1,ylim = c(0, 10))
ggsave(p1, file="2_PHF23_HDAC1_TSS_with_H3k27ac_boxplot.svg",width =10, height = 10,dpi=1080)
median(subset(PHF23_HDAC1_enhancer_H3K27ac,cluster=="co_binding")[,"fold_change"])
median(subset(PHF23_HDAC1_enhancer_H3K27ac,cluster=="only_PHF23")[,"fold_change"])
median(subset(PHF23_HDAC1_enhancer_H3K27ac,cluster=="only_HDAC1")[,"fold_change"])
median(subset(PHF23_HDAC1_enhancer_H3K27ac,cluster=="no_binding")[,"fold_change"])

#figure4L
PHF23_HDAC1_co_binding <- intersect(unique(na.omit(all_PHF23_TSS$SYMBOL)),unique(na.omit(all_HDAC1_TSS$SYMBOL)))
PHF23_HDAC1_co_binding <- as.data.frame(PHF23_HDAC1_co_binding)
PHF23_HDAC1_co_binding$cluster <- c("PHF23_HDAC1_co_binding")
PHF23_HDAC1_co_binding$entrez <- mapIds(x = org.Mm.eg.db,
                        keys = as.character(PHF23_HDAC1_co_binding$PHF23_HDAC1_co_binding),
            keytype ="SYMBOL",
            column ="ENTREZID",
            multiVals="first")
ee	<-as.matrix(PHF23_HDAC1_co_binding$entrez)
	dd <- as.vector(ee)
KEGGPHF23_HDAC1_co_binding <- enrichKEGG(gene =dd, 
					organism = "mouse", 
					keyType = "ncbi-geneid",
					 pvalueCutoff = 0.05,
				       pAdjustMethod = "BH", 
				       minGSSize = 10, 
				       maxGSSize = 500,
				       qvalueCutoff = 0.2, 
				       use_internal_data = TRUE)
p1 <- dotplot(KEGGPHF23_HDAC1_co_binding,showCategory=20)
ggsave(p1, file="2_PHF23_HDAC1_TSS_co_binding_KEGG_plot.svg",width =7, height = 8,dpi=1080)

#figure4M
my_plot <- venn.diagram(x=list(all_PHF23=all_PHF23_all$SYMBOL,TP53_pathway=TP53_pathway$MGI.symbol) , 
                          filename = NULL, 
                          fill=c("#F8766D", "#619CFF"), 
                          alpha=c(0.8,0.8),
                          euler.d = TRUE,
                          scaled = TRUE,
                          cat.fontfamily = "sans",
                          cat.cex = 1,
                          lwd = c(0.3, 0.3),
                          category.names=c("PHF23", "TP53_pathway"), 
                          #cat.just=list(c(1,-40) , c(0.4,-40)),
                          print.mode = c("raw", "percent"),
                          fontfamily = "sans",
                          main="Overlap between PHF23 and TP53_pathway",
                          main.fontfamily = "sans",
                          main.cex = 1,
                          cex = 0.9,
                          width = 10,
                          height = 10,
                          units = "in",
                          resolution = 400)
ggsave(my_plot, file="2_TP53_PHF23_plot_venn.svg", device = "svg")

#figure4N
TP53_without_PHF23_and_H3K4me3_ratio <- length(intersect(setdiff(TP53_pathway$MGI.symbol,all_PHF23_all$SYMBOL),all_H3K4me3_all$SYMBOL))/length(setdiff(TP53_pathway$MGI.symbol,all_PHF23_all$SYMBOL))
TP53_with_PHF23_and_H3K4me3_ratio <- length(intersect(intersect(TP53_pathway$MGI.symbol,all_PHF23_all$SYMBOL),all_H3K4me3_all$SYMBOL))/length(intersect(TP53_pathway$MGI.symbol,all_PHF23_all$SYMBOL))
H3K4me3_ratio <- data.frame(cluster=c("Without_PHF23","With_PHF23"),H3K4me3_ratio=c(TP53_without_PHF23_and_H3K4me3_ratio,TP53_with_PHF23_and_H3K4me3_ratio))
H3K4me3_ratio$H3K4me3_per <- round((H3K4me3_ratio$H3K4me3_ratio)*100,1)
H3K4me3_ratio$cluster <- factor(H3K4me3_ratio$cluster,levels=c("Without_PHF23","With_PHF23"))
my_plot <- ggplot(H3K4me3_ratio, aes(x=cluster, y=H3K4me3_per, fill=cluster)) +
  geom_bar(stat="identity")+ ylim(0,100) + scale_fill_manual(values=c("#00BFC4","#F8766D"))
ggsave(my_plot, file="3_TP53_with_without_PHF23_H3K4me3_levels_bar.svg", device = "svg")

#figure4O
TP53_without_PHF23_and_H3K27ac_ratio <- length(intersect(setdiff(TP53_pathway$MGI.symbol,all_PHF23_all$SYMBOL),all_H3K27ac_all$SYMBOL))/length(setdiff(TP53_pathway$MGI.symbol,all_PHF23_all$SYMBOL))
TP53_with_PHF23_and_H3K27ac_ratio <- length(intersect(intersect(TP53_pathway$MGI.symbol,all_PHF23_all$SYMBOL),all_H3K27ac_all$SYMBOL))/length(intersect(TP53_pathway$MGI.symbol,all_PHF23_all$SYMBOL))
H3K27ac_ratio <- data.frame(cluster=c("Without_PHF23","With_PHF23"),H3K27ac_ratio=c(TP53_without_PHF23_and_H3K27ac_ratio,TP53_with_PHF23_and_H3K27ac_ratio))
H3K27ac_ratio$H3K27ac_per <- round((H3K27ac_ratio$H3K27ac_ratio)*100,1)
H3K27ac_ratio$cluster <- factor(H3K27ac_ratio$cluster,levels=c("Without_PHF23","With_PHF23"))
my_plot <- ggplot(H3K27ac_ratio, aes(x=cluster, y=H3K27ac_per, fill=cluster)) +
  geom_bar(stat="identity")+ ylim(0,100) + scale_fill_manual(values=c("#00BFC4","#F8766D"))
ggsave(my_plot, file="3_TP53_with_without_PHF23_H3K27ac_levels_bar.svg", device = "svg")

#figureS5D
common_2 <- intersect(all_PHF23_TSS$SYMBOL,all_HDAC1_TSS$SYMBOL)
co_binding_with_H3K27ac <- round((length(intersect(all_H3K27ac_all$SYMBOL,common_2))/length(common_2))*100,1)
co_binding_without_H3K27ac <- 100-co_binding_with_H3K27ac
only_PHF23_with_H3K27ac <- round((length(intersect(all_H3K27ac_all$SYMBOL,setdiff(all_PHF23_TSS$SYMBOL,common_2)))/length(setdiff(all_PHF23_TSS$SYMBOL,common_2)))*100,1)
only_PHF23_without_H3K27ac <- 100-only_PHF23_with_H3K27ac
only_HDAC1_TSS_with_H3K27ac <- round((length(intersect(all_H3K27ac_all$SYMBOL,setdiff(all_HDAC1_TSS$SYMBOL,common_2)))/length(setdiff(all_HDAC1_TSS$SYMBOL,common_2)))*100,1)
only_HDAC1_TSS_without_H3K27ac <- 100-only_HDAC1_TSS_with_H3K27ac
no_binding_with_H3K27ac <- round((length(intersect(all_H3K27ac_all$SYMBOL,setdiff(na.omit(RNAseq_PreB$symbol),union(all_PHF23_TSS$SYMBOL,all_HDAC1_TSS$SYMBOL))))/length(setdiff(na.omit(RNAseq_PreB$symbol),union(all_PHF23_TSS$SYMBOL,all_HDAC1_TSS$SYMBOL))))*100,1)
no_binding_without_H3K27ac <- 100-no_binding_with_H3K27ac
HDAC1_TSS_H3K27ac <- data.frame(cluster=c("co_binding","co_binding","only_PHF23","only_PHF23","only_HDAC1_TSS","only_HDAC1_TSS","no_binding","no_binding"),ratio=c(co_binding_with_H3K27ac,co_binding_without_H3K27ac,only_PHF23_with_H3K27ac,only_PHF23_without_H3K27ac,only_HDAC1_TSS_with_H3K27ac,only_HDAC1_TSS_without_H3K27ac,no_binding_with_H3K27ac,no_binding_without_H3K27ac),H3K27ac=c("with","without","with","without","with","without","with","without"))
HDAC1_TSS_H3K27ac$ratio <- as.numeric(as.character(HDAC1_TSS_H3K27ac$ratio))
HDAC1_TSS_H3K27ac$cluster <- factor(HDAC1_TSS_H3K27ac$cluster,levels=c("no_binding","only_HDAC1_TSS","only_PHF23","co_binding"))
HDAC1_TSS_H3K27ac$H3K27ac <- factor(HDAC1_TSS_H3K27ac$H3K27ac,levels=c("without","with"))
p1 <- ggplot(HDAC1_TSS_H3K27ac, aes(fill=H3K27ac, y=ratio, x=cluster)) + 
    geom_bar(position="fill", stat="identity")+ scale_fill_manual(values=c("#00BFC4","#F8766D"))
ggsave(p1, file="3_PHF23_HDAC1_TSS_with_H3k27ac_ratio_Stacked.svg", device = "svg")

common_2 <- intersect(all_PHF23_TSS$SYMBOL,all_HDAC1_enhancer$SYMBOL)
co_binding_with_H3K27ac <- round((length(intersect(all_H3K27ac_all$SYMBOL,common_2))/length(common_2))*100,1)
co_binding_without_H3K27ac <- 100-co_binding_with_H3K27ac
only_PHF23_with_H3K27ac <- round((length(intersect(all_H3K27ac_all$SYMBOL,setdiff(all_PHF23_TSS$SYMBOL,common_2)))/length(setdiff(all_PHF23_TSS$SYMBOL,common_2)))*100,1)
only_PHF23_without_H3K27ac <- 100-only_PHF23_with_H3K27ac
only_HDAC1_enhancer_with_H3K27ac <- round((length(intersect(all_H3K27ac_all$SYMBOL,setdiff(all_HDAC1_enhancer$SYMBOL,common_2)))/length(setdiff(all_HDAC1_enhancer$SYMBOL,common_2)))*100,1)
only_HDAC1_enhancer_without_H3K27ac <- 100-only_HDAC1_enhancer_with_H3K27ac
no_binding_with_H3K27ac <- round((length(intersect(all_H3K27ac_all$SYMBOL,setdiff(na.omit(RNAseq_PreB$symbol),union(all_PHF23_TSS$SYMBOL,all_HDAC1_enhancer$SYMBOL))))/length(setdiff(na.omit(RNAseq_PreB$symbol),union(all_PHF23_TSS$SYMBOL,all_HDAC1_enhancer$SYMBOL))))*100,1)
no_binding_without_H3K27ac <- 100-no_binding_with_H3K27ac
HDAC1_enhancer_H3K27ac <- data.frame(cluster=c("co_binding","co_binding","only_PHF23","only_PHF23","only_HDAC1_enhancer","only_HDAC1_enhancer","no_binding","no_binding"),ratio=c(co_binding_with_H3K27ac,co_binding_without_H3K27ac,only_PHF23_with_H3K27ac,only_PHF23_without_H3K27ac,only_HDAC1_enhancer_with_H3K27ac,only_HDAC1_enhancer_without_H3K27ac,no_binding_with_H3K27ac,no_binding_without_H3K27ac),H3K27ac=c("with","without","with","without","with","without","with","without"))
HDAC1_enhancer_H3K27ac$ratio <- as.numeric(as.character(HDAC1_enhancer_H3K27ac$ratio))
HDAC1_enhancer_H3K27ac$cluster <- factor(HDAC1_enhancer_H3K27ac$cluster,levels=c("no_binding","only_HDAC1_enhancer","only_PHF23","co_binding"))
HDAC1_enhancer_H3K27ac$H3K27ac <- factor(HDAC1_enhancer_H3K27ac$H3K27ac,levels=c("without","with"))
p1 <- ggplot(HDAC1_enhancer_H3K27ac, aes(fill=H3K27ac, y=ratio, x=cluster)) + 
    geom_bar(position="fill", stat="identity")+ scale_fill_manual(values=c("#00BFC4","#F8766D"))
ggsave(p1, file="3_PHF23_HDAC1_enhancer_with_H3k27ac_ratio_Stacked.svg", device = "svg")

#figureS5E
PHF23_HDAC1_TSS_PreB$mean <- apply(PHF23_HDAC1_TSS_PreB[,c("PreB_shRen_1_TPM","PreB_shRen_2_TPM","PreB_shRen_3_TPM")],1,mean)
PHF23_HDAC1_TSS_PreB$log_mean <- log(PHF23_HDAC1_TSS_PreB$mean+0.1,2)
compare_means(log_mean ~ cluster, data = PHF23_HDAC1_TSS_PreB)
PHF23_HDAC1_TSS_PreB$cluster <- factor(PHF23_HDAC1_TSS_PreB$cluster,levels=c("no_binding_RNA","only_HDAC1_TSS_RNA","only_PHF23_RNA","co_binding_RNA"))
p1 <- ggboxplot(PHF23_HDAC1_TSS_PreB, x = "cluster", y = "log_mean",
          color = "cluster", palette = "jco",ylim = c(-6, 6),
         short.panel.labs = FALSE,size=1)
ggsave(p1, file="3_PHF23_HDAC1_TSS_PreB_basemean_binding_ggboxplot.svg",width =10, height = 10,dpi=1080)
median(subset(PHF23_HDAC1_TSS_PreB,cluster=="no_binding_RNA")[,"log_mean"])
median(subset(PHF23_HDAC1_TSS_PreB,cluster=="only_PHF23_RNA")[,"log_mean"])
median(subset(PHF23_HDAC1_TSS_PreB,cluster=="only_HDAC1_TSS_RNA")[,"log_mean"])
median(subset(PHF23_HDAC1_TSS_PreB,cluster=="co_binding_RNA")[,"log_mean"])

#figureS5F
PHF23_HDAC1_enhancer_PreB$mean <- apply(PHF23_HDAC1_enhancer_PreB[,c("PreB_shRen_1_TPM","PreB_shRen_2_TPM","PreB_shRen_3_TPM")],1,mean)
PHF23_HDAC1_enhancer_PreB$log_mean <- log(PHF23_HDAC1_enhancer_PreB$mean+0.1,2)
compare_means(log_mean ~ cluster, data = PHF23_HDAC1_enhancer_PreB)
PHF23_HDAC1_enhancer_PreB$cluster <- factor(PHF23_HDAC1_enhancer_PreB$cluster,levels=c("no_binding_RNA","only_HDAC1_enhancer_RNA","only_PHF23_RNA","co_binding_RNA"))
p1 <- ggboxplot(PHF23_HDAC1_enhancer_PreB, x = "cluster", y = "log_mean",
          color = "cluster", palette = "jco",ylim = c(-5, 6),
         short.panel.labs = FALSE,size=1)
ggsave(p1, file="3_PHF23_HDAC1_enhancer_PreB_basemean_binding_ggboxplot.svg",width =10, height = 10,dpi=1080)
median(subset(PHF23_HDAC1_enhancer_PreB,cluster=="no_binding_RNA")[,"log_mean"])
median(subset(PHF23_HDAC1_enhancer_PreB,cluster=="only_PHF23_RNA")[,"log_mean"])
median(subset(PHF23_HDAC1_enhancer_PreB,cluster=="only_HDAC1_enhancer_RNA")[,"log_mean"])
median(subset(PHF23_HDAC1_enhancer_PreB,cluster=="co_binding_RNA")[,"log_mean"])

#figureS5G
my_plot <- venn.diagram(x=list(all_H3K4me3=all_H3K4me3_all$SYMBOL,all_H3K27ac=all_H3K27ac_all$SYMBOL) , 
                          filename = NULL, 
                          fill=c("#F8766D", "#619CFF"), 
                          alpha=c(0.8,0.8),
                          euler.d = TRUE,
                          scaled = TRUE,
                          cat.fontfamily = "sans",
                          cat.cex = 1,
                          lwd = c(0.3, 0.3),
                          category.names=c("H3K4me3", "H3K27ac"), 
                          #cat.just=list(c(1,-40) , c(0.4,-40)),
                          print.mode = c("raw", "percent"),
                          fontfamily = "sans",
                          main="Overlap between H3K4me3 and H3K27ac",
                          main.fontfamily = "sans",
                          main.cex = 1,
                          cex = 0.9,
                          width = 10,
                          height = 10,
                          units = "in",
                          resolution = 400)
ggsave(my_plot, file="3_H3K4me3_H3K27ac_venn.svg", device = "svg")

#figureS5H
H3K4me3_H3K27ac_PreB$mean <- apply(H3K4me3_H3K27ac_PreB[,c("PreB_shRen_1_TPM","PreB_shRen_2_TPM","PreB_shRen_3_TPM")],1,mean)
H3K4me3_H3K27ac_PreB$log_mean <- log(H3K4me3_H3K27ac_PreB$mean+0.1,2)
H3K4me3_H3K27ac_PreB$log_mean <- as.numeric(as.character(scale(H3K4me3_H3K27ac_PreB$log_mean)))
compare_means(log_mean ~ cluster, data = H3K4me3_H3K27ac_PreB)
H3K4me3_H3K27ac_PreB$cluster <- factor(H3K4me3_H3K27ac_PreB$cluster,levels=c("no_binding_RNA","only_H3K4me3_RNA","only_H3K27ac_RNA","co_binding_RNA"))
p1 <- ggboxplot(H3K4me3_H3K27ac_PreB, x = "cluster", y = "log_mean",
          color = "cluster", palette = "jco",ylim = c(-1, 5),
         short.panel.labs = FALSE,size=1)
ggsave(p1, file="3_H3K4me3_H3K27ac_PreB_basemean_binding_ggboxplot.svg",width =10, height = 10,dpi=1080)
median(subset(H3K4me3_H3K27ac_PreB,cluster=="no_binding_RNA")[,"log_mean"])
median(subset(H3K4me3_H3K27ac_PreB,cluster=="only_H3K27ac_RNA")[,"log_mean"])
median(subset(H3K4me3_H3K27ac_PreB,cluster=="only_H3K4me3_RNA")[,"log_mean"])
median(subset(H3K4me3_H3K27ac_PreB,cluster=="co_binding_RNA")[,"log_mean"])

#figureS5I
common_2 <- intersect(all_H3K4me3_all$SYMBOL,all_H3K27ac_all$SYMBOL)
my_plot <- venn.diagram(x=list(all_PHF23=all_PHF23_all$SYMBOL,common_2=common_2,all_HDAC1=all_HDAC1_all$SYMBOL) , 
                          filename = NULL, 
                          fill=c("#F8766D", "#619CFF","#d986ed"), 
                          alpha=c(0.8,0.8,0.8),
                          euler.d = TRUE,
                          scaled = TRUE,
                          cat.fontfamily = "sans",
                          cat.cex = 1,
                          lwd = c(0.3, 0.3, 0.3),
                          category.names=c("PHF23", "common_2", "HDAC1"), 
                          #cat.just=list(c(1,-40) , c(0.4,-40)),
                          print.mode = c("raw", "percent"),
                          fontfamily = "sans",
                          main="Overlap between PHF23 and H3K4me3_H3K27ac and HDAC1",
                          main.fontfamily = "sans",
                          main.cex = 1,
                          cex = 0.9,
                          width = 10,
                          height = 10,
                          units = "in",
                          resolution = 400)
ggsave(my_plot, file="3_PHF23_H3K4me3_and_H3K27ac_HDAC1_plot_venn.svg", device = "svg")

#figureS5J
H3K4me3_and_H3K27ac_PHF23_HDAC1_PreB$mean <- apply(H3K4me3_and_H3K27ac_PHF23_HDAC1_PreB[,c("PreB_shRen_1_TPM","PreB_shRen_2_TPM","PreB_shRen_3_TPM")],1,mean)
H3K4me3_and_H3K27ac_PHF23_HDAC1_PreB$log_mean <- log(H3K4me3_and_H3K27ac_PHF23_HDAC1_PreB$mean+0.1,2)
H3K4me3_and_H3K27ac_PHF23_HDAC1_PreB$log_mean <- as.numeric(as.character(scale(H3K4me3_and_H3K27ac_PHF23_HDAC1_PreB$log_mean)))
compare_means(log_mean ~ cluster, data = H3K4me3_and_H3K27ac_PHF23_HDAC1_PreB)
H3K4me3_and_H3K27ac_PHF23_HDAC1_PreB$cluster <- factor(H3K4me3_and_H3K27ac_PHF23_HDAC1_PreB$cluster,levels=c("H3K4me3_H3K27ac_RNA","only_HDAC1_RNA","only_PHF23_RNA","co_binding_RNA"))
p1 <- ggboxplot(H3K4me3_and_H3K27ac_PHF23_HDAC1_PreB, x = "cluster", y = "log_mean",
          color = "cluster", palette = "jco",ylim = c(-2,4),
         short.panel.labs = FALSE,size=1)
ggsave(p1, file="3_H3K4me3_and_H3K27ac_PHF23_HDAC1_PreB_basemean_binding_ggboxplot.svg",width =10, height = 10,dpi=1080)
median(subset(H3K4me3_and_H3K27ac_PHF23_HDAC1_PreB,cluster=="H3K4me3_H3K27ac_RNA")[,"log_mean"])
median(subset(H3K4me3_and_H3K27ac_PHF23_HDAC1_PreB,cluster=="only_PHF23_RNA")[,"log_mean"])
median(subset(H3K4me3_and_H3K27ac_PHF23_HDAC1_PreB,cluster=="only_HDAC1_RNA")[,"log_mean"])
median(subset(H3K4me3_and_H3K27ac_PHF23_HDAC1_PreB,cluster=="co_binding_RNA")[,"log_mean"])

#figureS5K
my_plot <- venn.diagram(x=list(all_PHF23=all_PHF23_all$SYMBOL,TP53_pathway=TP53_pathway$MGI.symbol,all_H3K4me3=all_H3K4me3_all$SYMBOL) , 
                          filename = NULL, 
                          fill=c("#F8766D", "#619CFF","#d986ed"), 
                          alpha=c(0.8,0.8,0.8),
                          euler.d = TRUE,
                          scaled = TRUE,
                          cat.fontfamily = "sans",
                          cat.cex = 1,
                          lwd = c(0.3, 0.3, 0.3),
                          category.names=c("PHF23", "TP53_pathway", "H3K4me3"), 
                          #cat.just=list(c(1,-40) , c(0.4,-40)),
                          print.mode = c("raw", "percent"),
                          fontfamily = "sans",
                          main="Overlap between PHF23 and TP53_pathway and H3K4me3",
                          main.fontfamily = "sans",
                          main.cex = 1,
                          cex = 0.9,
                          width = 10,
                          height = 10,
                          units = "in",
                          resolution = 400)
ggsave(my_plot, file="3_TP53_PHF23_H3K4me3_venn.svg", device = "svg")

#figureS5L
my_plot <- venn.diagram(x=list(all_PHF23=all_PHF23_all$SYMBOL,TP53_pathway=TP53_pathway$MGI.symbol,all_H3K27ac=all_H3K27ac_all$SYMBOL) , 
                          filename = NULL, 
                          fill=c("#F8766D", "#619CFF","#d986ed"), 
                          alpha=c(0.8,0.8,0.8),
                          euler.d = TRUE,
                          scaled = TRUE,
                          cat.fontfamily = "sans",
                          cat.cex = 1,
                          lwd = c(0.3, 0.3, 0.3),
                          category.names=c("PHF23", "TP53_pathway", "H3K27ac"), 
                          #cat.just=list(c(1,-40) , c(0.4,-40)),
                          print.mode = c("raw", "percent"),
                          fontfamily = "sans",
                          main="Overlap between PHF23 and TP53_pathway and H3K27ac",
                          main.fontfamily = "sans",
                          main.cex = 1,
                          cex = 0.9,
                          width = 10,
                          height = 10,
                          units = "in",
                          resolution = 400)
ggsave(my_plot, file="3_TP53_PHF23_H3K27ac_venn.svg", device = "svg")

