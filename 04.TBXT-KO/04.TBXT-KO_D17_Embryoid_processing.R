# Code used to analysis the scRNA-seq data of TBXT-KO D17 embryoid

################### 1. load required packages
###################
###################
###################
###################
###################
library(Seurat)
library(dplyr)
library(ggplot2)
library(SCpubr)
library(clustree)

dir.create("/data/embryo_dev/myprocess/T-2-D17/integrated_d17-1")
setwd("/data/embryo_dev/myprocess/T-2-D17/integrated_d17-1")


load("/data/embryo_dev/myprocess/integrate_cr6ncbi_d17_d18_rmC1C6/only_newdata/d17-1_d18-1_integrate_analysis_231111.RData")
rm(list=setdiff(ls(),"sc_combined"))

d17 <- subset(sc_combined, orig.ident=="d17_1")
#d18 <- subset(sc_combined, orig.ident=="d18_1")

meta.data <- d17@meta.data
meta.data <- meta.data[,c(1,2,3,21,22)]
d17 <- CreateSeuratObject(counts = d17@assays$RNA@counts,meta.data = meta.data)

#meta.data <- d18@meta.data
#meta.data <- meta.data[,c(1,2,3,21,22)]
#d18 <- CreateSeuratObject(counts = d18@assays$RNA@counts,meta.data = meta.data)

td17 <- readRDS("/data/embryo_dev/myprocess/T-2-D17/cluster/T-2-D17_rna3_seuobj_250401.rds")

meta.data <- td17@meta.data
meta.data <- meta.data[,c(1,2,3,22,23)]
meta.data$celltype <- paste0("T2-D17--",meta.data$celltype)
meta.data$celltype2 <- paste0("T2-D17--",meta.data$celltype2)
td17 <- CreateSeuratObject(counts = td17@assays$RNA@counts,meta.data = meta.data)

sc_list <- list(d17=d17, td17=td17)

# normalize and identify variable features for each dataset independently
sc_list <- lapply(X = sc_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})
features <- SelectIntegrationFeatures(object.list = sc_list, nfeatures = 3000)
sc_anchors <- FindIntegrationAnchors(object.list = sc_list, anchor.features = features)
sc_combined <- IntegrateData(anchorset = sc_anchors)

DefaultAssay(sc_combined) <- "integrated"
# Run the standard workflow for visualization and clustering
sc_combined <- ScaleData(sc_combined, verbose = FALSE)
sc_combined <- RunPCA(sc_combined, npcs = 100, verbose = FALSE)
ElbowPlot(sc_combined,ndims = 50, reduction = "pca")
sc_combined <- RunUMAP(sc_combined, reduction = "pca", dims = 1:30, umap.method = "uwot", metric = "cosine")
sc_combined <- RunTSNE(sc_combined, reduction = "pca", dims = 1:30)

DimPlot(sc_combined, group.by = "orig.ident", raster = FALSE)

VlnPlot(sc_combined, features = c("nFeature_RNA"), group.by = "orig.ident", pt.size = 0)

sc_combined <- FindNeighbors(sc_combined, reduction = "pca", dims = 1:30)
for(i in seq(0.3,1.2,by=0.1)){
  sc_combined<- FindClusters(sc_combined,verbose = FALSE,resolution=i)
}

sc_combined<- FindClusters(sc_combined,verbose = FALSE,resolution=1.3)
sc_combined<- FindClusters(sc_combined,
                           verbose = FALSE,
                           resolution=1.4)

sc_combined<- FindClusters(sc_combined,
                           verbose = FALSE,
                           resolution=1.8)

library(clustree)
ptree <- clustree(sc_combined, prefix = "integrated_snn_res.")
ptree

p1 <- DimPlot(sc_combined, 
              group.by = "integrated_snn_res.0.9",
        raster = FALSE, pt.size = 0.1,
        label = T) +NoLegend()
p2 <- DimPlot(sc_combined, 
              reduction = "tsne",
              group.by = "integrated_snn_res.0.9",
              raster = FALSE, pt.size = 0.1,
              label = T)+NoLegend()
p1|p2




# library(SCpubr)
# do_DimPlot(sc_combined, group.by = "integrated_snn_res.1", pt.size = 0.1,
#            idents.keep = "27")
# do_DimPlot(sc_combined, group.by = "integrated_snn_res.1", pt.size = 0.1,
#            idents.keep = "28")
# DimPlot(sc_combined, group.by = "integrated_snn_res.1",
#         raster = FALSE, pt.size = 0.1,
#         split.by = "orig.ident",
#         label = T) +NoLegend()


DefaultAssay(sc_combined) <- "RNA"

p1 <- FeaturePlot(sc_combined, 
                  features = c("T"),
                  raster = FALSE, pt.size = 0.1) 
p2 <- FeaturePlot(sc_combined, 
                  reduction = "tsne",
                  features = c("T"),
                  raster = FALSE, pt.size = 0.1)
p1|p2

VlnPlot(sc_combined, features = "T", group.by = "orig.ident") +NoLegend()

Idents(sc_combined) <- "integrated_snn_res.0.9"
DefaultAssay(sc_combined) <- "RNA"
sc_combined_res09_markers <- FindAllMarkers(sc_combined, only.pos = TRUE)



metadf <- sc_combined@meta.data
ery_cells <- row.names(metadf)[metadf$integrated_snn_res.1.8=="35"]
neu_cells <- row.names(metadf)[metadf$integrated_snn_res.1.8=="36"]

library(dplyr)
anno <- sc_combined$integrated_snn_res.0.9
anno <- recode(anno,
               "0"="EPI","1"="AM","2"="Amn","3"="EXMC",
               "4"="EXMC","5"="Gut","6"="EXMC","7"="TE",
               "8"="EXMC","9"="EPI","10"="EXMC","11"="NNE",
               "12"="AM","13"="NM","14"="Gast","15"="EXMC",
               "16"="DE/VE","17"="EPI","18"="EXMC",
               "19"="EXMC", "20"="EC","21"="Gut",
               "22"="PGC", "23"="ysEndo","24"="EXMC",
               "25"="TE","26"="Ery_Neu","27"="Mac")
anno <- as.character(anno)
anno[match(ery_cells, row.names(metadf))] <- "Ery"
anno[match(neu_cells, row.names(metadf))] <- "Neu"
names(anno) <- row.names(metadf)
anno <-factor(anno, levels = unique(anno))

sc_combined$anno <- anno
DimPlot(sc_combined, label = TRUE, pt.size = 0.01, 
        raster = FALSE,
        split.by = "orig.ident",
        group.by = "anno") +NoLegend() 
  
  
VlnPlot(sc_combined, 
        features = c("nFeature_RNA"),
        group.by = "orig.ident", pt.size = 0)


FeaturePlot(sc_combined, features = c("T"),
            split.by = "orig.ident",raster = FALSE)


anno_freq <- data.frame(table(sc_combined$anno, sc_combined$orig.ident))
colnames(anno_freq) <- c("cluster","sample","number")

anno_freq$percent <- anno_freq$number
anno_freq$percent[anno_freq$sample=="d17_1"] <- 
  anno_freq$percent[anno_freq$sample=="d17_1"] /sum(anno_freq$percent[anno_freq$sample=="d17_1"]) *100
anno_freq$percent[anno_freq$sample=="T-2-D17"] <- 
  anno_freq$percent[anno_freq$sample=="T-2-D17"] /sum(anno_freq$percent[anno_freq$sample=="T-2-D17"]) *100

ggplot(anno_freq[anno_freq$cluster%in%c("Gast","NM"),],
       aes(x=cluster,y=percent,fill=sample))+
  geom_bar(stat="identity",position = position_dodge2())

ggplot(anno_freq,
       aes(x=cluster,y=percent,fill=sample))+
  geom_bar(stat="identity",position = position_dodge2())

uni_anno <- unique(anno_freq$cluster)
uni_anno_log2change <- log2(anno_freq$percent[anno_freq$sample=="T-2-D17"]/
  anno_freq$percent[anno_freq$sample=="d17_1"])


uni_anno_log2change_df <- data.frame(cluster=levels(anno_freq$cluster),
                                     fc = uni_anno_log2change)
uni_anno_log2change_df$group <- rep("up",nrow(uni_anno_log2change_df))
uni_anno_log2change_df$group[uni_anno_log2change_df$fc<0] <- "down"

ggplot(uni_anno_log2change_df,aes(x=cluster,y=fc,colour=group)) +
  geom_point(size=6) +coord_flip() +ylab(label = "log2(KO/WT)") +xlab(label = "") +
  theme_classic() +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black")) +
  geom_hline(yintercept = 0,color="black") +NoLegend()

# Idents(sc_combined) <- "integrated_snn_res.1"
# DefaultAssay(sc_combined) <- "RNA"
# sc_combined_res1_markers <- FindAllMarkers(sc_combined, only.pos = TRUE)

Idents(sc_combined) <- "anno"
DefaultAssay(sc_combined) <- "RNA"
sc_combined_anno_markers <- FindAllMarkers(sc_combined, only.pos = TRUE)

csc_markers <- c("POU5F1","SOX2","DNMT3B","T","MIXL1", "ISL1",
                 "TFAP2A","TFAP2C","NANOS3", "PRDM1",
                 "SOX17","GATA4","FOXA2","OTX2","CDX2","OSR1",
                 "DPPA3","POSTN", "COL6A1","KRT7","KRT8","NR2F2",
                 "MESP1","EOMES",
                 "CENPV","HGF",
                 "SLC7A3","EFHD1",
                 "KDR","PECAM1",
                 "C1QA","C1QB","C1QC",
                 "APOB","AFP","TTR",
                 "S100A9","S100A8","MPO","LYZ",
                 "HBM","GYPA","HBZ","HBE1")
DotPlot(sc_combined, features = csc_markers,
        cols=c("grey","#E7241A","#940001")) + 
  theme(axis.text.x = element_text(angle=45,vjust = 1,hjust = 1)) +labs(x="")+
  scale_size_continuous(range = c(0,6))

VlnPlot(sc_combined, features = c("MESP1","TBX6","RHOA","CDC42"),
        split.by = "orig.ident",
        stack = TRUE, flip = TRUE)+
  theme(axis.text.x = element_text(angle=45,vjust = 1,hjust = 1)) +labs(x="")+
  scale_size_continuous(range = c(0,6))



# calculate DEGs between KO and WT in Gasting cells
gast <- subset(sc_combined,anno %in% c("Gast","NM"))
Idents(gast) <- "orig.ident"
gast_degs <- FindMarkers(gast,ident.1="d17_1",ident.2="T-2-D17")
gast_degs$gene <- row.names(gast_degs)
gast_degs_sig <- gast_degs[(gast_degs$p_val_adj<0.05 & abs(gast_degs$avg_log2FC)>=1),]

endo <- subset(sc_combined,anno %in% c("DE/VE","Gut","ysEndo"))
Idents(endo) <- "orig.ident"
endo_degs <- FindMarkers(endo,ident.1="d17_1",ident.2="T-2-D17")
endo_degs$gene <- row.names(endo_degs)
endo_degs_sig <- endo_degs[(endo_degs$p_val_adj<0.05 & abs(endo_degs$avg_log2FC)>=1),]


# 
library(readxl)
tGlist1 <- read_excel("/data/embryo_dev/myprocess/T-2-D17/TBXT_downstream_deepseek.xlsx")
tGlist1 <- tGlist1$gene

tgf <- read.csv("/data/embryo_dev/myprocess/T-2-D17/KEGG_TGF_BETA_SIGNALING_PATHWAY.v2024.1.Hs.tsv",
                sep="\t")
tgfGlist <- unlist(strsplit(tgf[17,2],","))
wnt <- read.csv("/data/embryo_dev/myprocess/T-2-D17/KEGG_WNT_SIGNALING_PATHWAY.v2024.1.Hs.tsv",
                sep="\t")
wntGlist <- unlist(strsplit(wnt[17,2],","))

intersect(gast_degs_sig$gene, tGlist1)
intersect(gast_degs_sig$gene, tgfGlist)
intersect(gast_degs_sig$gene, wntGlist)

VlnPlot(gast, features = c("CTNNB1","ID3","ID1"), stack = T,flip = T) +NoLegend()

intersect(gast_degs$gene[gast_degs$p_val_adj<0.05], tGlist1)
intersect(gast_degs$gene[gast_degs$p_val_adj<0.05], tgfGlist)
intersect(gast_degs$gene[gast_degs$p_val_adj<0.05], wntGlist)


VlnPlot(gast, features = intersect(gast_degs$gene[gast_degs$p_val_adj<0.05], tGlist1), stack = T,flip = T) +NoLegend()
VlnPlot(gast, features = intersect(gast_degs$gene[gast_degs$p_val_adj<0.05], tgfGlist), stack = T,flip = T) +NoLegend()
VlnPlot(gast, features = intersect(gast_degs$gene[gast_degs$p_val_adj<0.05], wntGlist), stack = T,flip = T) +NoLegend()


gast <- ScaleData(gast,features = row.names(gast))
#DotPlot(gast,group.by = "orig.ident",features = gast_degs_sig$gene)

VlnPlot(endo, features = intersect(endo_degs$gene[endo_degs$p_val_adj<0.05], tGlist1), stack = T,flip = T) +NoLegend()
VlnPlot(endo, features = intersect(endo_degs$gene[endo_degs$p_val_adj<0.05], tgfGlist), stack = T,flip = T) +NoLegend()
VlnPlot(endo, features = intersect(endo_degs$gene[endo_degs$p_val_adj<0.05], wntGlist), stack = T,flip = T) +NoLegend()


write.table(gast_degs, file="T-2-D17_D17-1_DEGS_gastNM.tsv",sep="\t",row.names = FALSE, quote=FALSE)
write.table(endo_degs, file="T-2-D17_D17-1_DEGS_endodermlineage.tsv",sep="\t",row.names = FALSE, quote=FALSE)



# functional analysis
library(clusterProfiler)
library(org.Hs.eg.db)

degs_list <- list(gast_KO_UP = gast_degs_sig$gene[gast_degs_sig$avg_log2FC>0],
                  gast_KO_DOWN = gast_degs_sig$gene[gast_degs_sig$avg_log2FC<0],
                  endo_KO_UP = endo_degs_sig$gene[endo_degs_sig$avg_log2FC>0],
                  endo_KO_DOWN = endo_degs_sig$gene[endo_degs_sig$avg_log2FC<0])

# transform mfas gene to human gene
feature_df <- read.csv("/data/embryo_dev/mfas_blastoid_cr6.1.2_ncbi/features_hg19_update.tsv",
                       sep="\t",header = T)

degs_list_hEntrez <- lapply(degs_list,function(x){
  entrez <- unique(na.omit(feature_df$hg19_entrez_id[feature_df$macFas5_gene_symbol %in% x]))
  return(entrez)
})

degs_list_CompareCluster_BP <- compareCluster(geneClusters = degs_list_hEntrez,
                                                             fun = "enrichGO",
                                                             OrgDb = "org.Hs.eg.db", 
                                                             ont ="BP", 
                                                             pvalueCutoff = 0.05, 
                                                             pAdjustMethod = "BH", 
                                                             qvalueCutoff = 0.2,
                                                             readable=T)
write.table(degs_list_CompareCluster_BP@compareClusterResult, 
            file="T-2-D17_D17-1_DEGS_BP.tsv",sep="\t",row.names = FALSE, quote=FALSE)

save.image("/data/embryo_dev/myprocess/T-2-D17/integrated_d17-1/T-2-D17_D17-1_IntegratedAnalysis_250406.RData")



tExp <- FetchData(sc_combined,vars = c("orig.ident","celltype","celltype2","anno","T"))
#
table(tExp$orig.ident)
sum(tExp$orig.ident=="d17_1" & tExp$T > 0)
sum(tExp$orig.ident=="T-2-D17" & tExp$T > 0)

347/11441
190/9170

barplot(c(3.03,2.07),names.arg = c("T-2-D17","D17-1"), ylim = c(0,3.5))



eryneu <- subset(sc_combined, anno=="Ery_Neu")
eryneu <- NormalizeData(eryneu)
eryneu <- FindVariableFeatures(eryneu, selection.method = "vst", nfeatures = 3000)
eryneu <- ScaleData(eryneu, verbose = FALSE)
eryneu <- RunPCA(eryneu, npcs = 100, verbose = FALSE)
ElbowPlot(eryneu,ndims = 50, reduction = "pca")
eryneu <- RunUMAP(eryneu, reduction = "pca", dims = 1:30, umap.method = "uwot", metric = "cosine")
eryneu <- RunTSNE(eryneu, reduction = "pca", dims = 1:30)
eryneu <- FindNeighbors(eryneu)
eryneu <- FindClusters(eryneu, resolution = 0.5 )

DimPlot(eryneu, group.by = "orig.ident")
DimPlot(eryneu)



DotPlot(eryneu, features = c("S100A9","LYZ","HBZ","HBM"))

FeaturePlot(eryneu, features = c("S100A9","LYZ"))
FeaturePlot(eryneu, features = c("HBZ","HBM1"))

Idents(eryneu) <-   eryneu$RNA_snn_res.0.2
eryneu   <- RenameIdents(object = eryneu, 
                      '0' = 'Neu',
                      '1' = 'Neu',
                      '2' = 'Ery')
DimPlot(eryneu)

# sum(tExp[tExp$anno%in%c("Gast","NM"),"T"]>0)
# sum(tExp$anno%in%c("Gast","NM") & tExp$orig.ident=="d17_1")
# sum(tExp[tExp$anno%in%c("Gast","NM") & tExp$orig.ident=="d17_1","T"]>0)
# 
# sum(tExp$anno%in%c("Gast","NM") & tExp$orig.ident=="d18_1")
# sum(tExp[tExp$anno%in%c("Gast","NM") & tExp$orig.ident=="d18_1","T"]>0)
# sum(tExp[tExp$anno%in%c("Gast","NM") & tExp$orig.ident=="T-2-D17","T"]>0)
# sum(tExp$anno%in%c("Gast","NM") & tExp$orig.ident=="T-2-D17")
# 
# 
# dat <- data.frame(sample=unique(tExp$orig.ident),
#                   Tpos = c(237/832, 339/842, 115/473),
#                   Tpos2 = c(347/11441,452/10115, 190/9170)*100,
#                   Tpos3 = c(347,453,190))
# dat$sample <- factor(dat$sample,levels = c("T-2-D17","d17_1","d18_1"))
# p1 <- ggplot(dat, aes(x=sample,y=Tpos, fill=sample))+geom_bar(stat="identity") +ylab(label = "Percentage of T+ cells in Gast/NM") +xlab(label = "")
# p2 <- ggplot(dat, aes(x=sample,y=Tpos2, fill=sample))+geom_bar(stat="identity") +ylab(label = "Percentage of T+ cells")+xlab(label = "")
# p3 <- ggplot(dat, aes(x=sample,y=Tpos3, fill=sample))+geom_bar(stat="identity") + ylab(label = "Number of T+ cells")+xlab(label = "")
# library(cowplot)
# 
# plot_grid(p1,p2,p3, nrow=3)
# 
# 
# S100P

# VlnPlot(sc_combined, features = csc_markers,
#         split.by = "orig.ident",
#         stack = TRUE, flip = TRUE,
#         cols=c("grey","#E7241A","#940001"))
# 
# + 
#   theme(axis.text.x = element_text(angle=45,vjust = 1,hjust = 1)) +labs(x="")+
#   scale_size_continuous(range = c(0,6))
# 
library(RColorBrewer)


setwd("/data/embryo_dev/myprocess/T-2-D17/integrated_d17-1")
load("/data/embryo_dev/myprocess/T-2-D17/integrated_d17-1/T-2-D17_D17-1_IntegratedAnalysis_250406.RData")

dir.create("figures_integration_250408")

colors1 <- c("#23908A","#FDA505")
p <- VlnPlot(sc_combined, features = c("nFeature_RNA"), 
        group.by = "orig.ident", pt.size = 0,
        cols = colors1) +
  xlab(label = "") +
  labs(title = "")+
  ylab(label = "Number of Detected Genes") + NoLegend() +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color="black"))
ggsave(p, file="figures_integration_250408/NumberofGenesPerSample_vlnplot.pdf",
       height = 4.5, width = 3)

p1 <- DimPlot(sc_combined, group.by = "orig.ident", 
             pt.size = 0.1, cols = colors1,
             reduction = "tsne",
             raster = FALSE) +
  xlab(label = "TSNE1") +
  labs(title = "")+
  ylab(label = "TSNE2") +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color="black"))
p2 <- DimPlot(sc_combined, group.by = "integrated_snn_res.0.9", reduction = "tsne",
              pt.size = 0.1, label = TRUE,
              raster = FALSE) + NoLegend() +
  xlab(label = "TSNE1") + labs(title = "") + ylab(label = "TSNE2") +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color="black"))
pp <- p1|p2
ggsave(pp, file="figures_integration_250408/TSNE_cplorbySample_seuratclusters.pdf",
       height = 4.5, width = 11)


p3 <- DimPlot(sc_combined, group.by = "anno", reduction = "tsne",
              split.by = "orig.ident",
              pt.size = 0.1, label = FALSE,
              raster = FALSE)  +
  xlab(label = "TSNE1") + labs(title = "") + ylab(label = "TSNE2") +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color="black"))
ggsave(p3, file="figures_integration_250408/TSNE_cplorbyAnnoClusters.pdf",
       height = 4.5, width = 10)

p4 <- ggplot(uni_anno_log2change_df,aes(x=cluster,y=fc,colour=group)) +
  geom_point(size=4) +coord_flip() +ylab(label = "log2(KO/WT)") +xlab(label = "") +
  theme_classic() +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black")) +
  geom_hline(yintercept = 0,color="black") +NoLegend()
ggsave(p4, file="figures_integration_250408/CellRatio.pdf",
       height = 4.5, width = 7)

anno_freq$log10CellNumber <- log10(anno_freq$number)
anno_freq$log10CellNumber[anno_freq$log10CellNumber=="-Inf"] <- 0
p5 <- ggplot(anno_freq,aes(x=cluster,y=log10CellNumber, color=sample))+
  geom_point(size=3)+
  ylab(label = "log10(Count)") +xlab(label = "") +
  scale_color_manual(values = colors1)+
  theme_classic() +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black")) 
ggsave(p5, file="figures_integration_250408/CellCount.pdf",
       height = 4.5, width = 7)

color2 <- rev(c("Red",
                "Red",
                "Orange",
                "#ffffbf",
                "#e0f3f8",
                "#91bfdb",
                "#4575b4"))



DefaultAssay(sc_combined) <- "RNA"
p6 <- FeaturePlot(sc_combined, features = c("T","MIXL1","WNT5B","DLL3","MESP1","EOMES"),
                  cols = c("lightgrey","blue"),
            pt.size = 0.1, raster = FALSE, split.by = "orig.ident",
            reduction="tsne", order = TRUE )+
  xlab(label = "TSNE1") +  ylab(label = "TSNE2") +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color="black") )
ggsave(p6, file="figures_integration_250408/markers_tsne.pdf",
       height = 14, width = 5.5)

 
gast <- subset(sc_combined,anno %in% c("Gast"))
Idents(gast) <- "orig.ident"
gast_degs <- FindMarkers(gast,ident.1="d17_1",ident.2="T-2-D17")
gast_degs$gene <- row.names(gast_degs)


write.table(gast_degs, file = "T-2-D17_D17-1_Gast_DEGs.txt", sep="\t", 
            row.names = FALSE, quote=FALSE)


gast_degs_sig <- gast_degs[(gast_degs$p_val_adj<0.05 & abs(gast_degs$avg_log2FC)>=1),]

gast <- ScaleData(gast, features = row.names(gast))

DoHeatmap(gast,features = gast_degs_sig[order(gast_degs_sig$avg_log2FC,decreasing = T),"gene"],
          group.by = "orig.ident")
ggsave(filename = "figures_integration_250408/Gast_WTvsKO_logfc1_degs_heatmap.pdf",
       height = 14, width = 6)


FeaturePlot(sc_combined, features = c("MTPN","DPPA4","EIF5","SPP1","ID3","ID1"),
            cols = c("lightgrey","blue"),
            pt.size = 0.1, raster = FALSE, split.by = "orig.ident",
            reduction="tsne", order = TRUE )+
  xlab(label = "TSNE1") +  ylab(label = "TSNE2") +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color="black") )

VlnPlot(gast, features = c("EOMES"),pt.size = 0, cols = colors1)


degs_list <- list(gast_KO_UP = gast_degs_sig$gene[gast_degs_sig$avg_log2FC>0],
                  gast_KO_DOWN = gast_degs_sig$gene[gast_degs_sig$avg_log2FC<0])

# transform mfas gene to human gene
feature_df <- read.csv("/data/embryo_dev/mfas_blastoid_cr6.1.2_ncbi/features_hg19_update.tsv",
                       sep="\t",header = T)

degs_list_hEntrez <- lapply(degs_list,function(x){
  entrez <- unique(na.omit(feature_df$hg19_entrez_id[feature_df$macFas5_gene_symbol %in% x]))
  return(entrez)
})

degs_list_CompareCluster_BP <- compareCluster(geneClusters = degs_list_hEntrez,
                                              fun = "enrichGO",
                                              OrgDb = "org.Hs.eg.db", 
                                              ont ="BP", 
                                              pvalueCutoff = 0.05, 
                                              pAdjustMethod = "BH", 
                                              qvalueCutoff = 0.2,
                                              readable=T)
write.table(degs_list_CompareCluster_BP@compareClusterResult, 
            file="T-2-D17_D17-1_Gast_DEGS_BP.txt",sep="\t",row.names = FALSE, quote=FALSE)


save.image("/data/embryo_dev/myprocess/T-2-D17/integrated_d17-1/T-2-D17_D17-1_IntegratedAnalysis_250409.RData")

load("/data/embryo_dev/myprocess/T-2-D17/integrated_d17-1/T-2-D17_D17-1_IntegratedAnalysis_250409.RData")


anno2 <- sc_combined$anno
anno2 <- recode(anno2, 
                "Gast"="Gast", "AM"="AM","EXMC"="EXMC",
                "Amn"="Amn", "EPI"="EPI", "NNE"="NNE",
                "Gut"="Endoderm","NM"="NM","EC"="HEP",
                "DE/VE"="Endoderm","TE"="TE","PGC"="PGC",
                "ysEndo"="Endoderm","Neu"="Neu","Ery"="Ery","Mac"="Mac")
anno2 <- factor(anno2, levels = c("EPI","Gast","NNE","Amn","TE",
                                  "NM","AM","EXMC",
                                  "Endoderm","PGC","HEP","Neu","Mac","Ery"))
sc_combined$anno2 <- anno2

dir.create("figures_integration_250410")

anno2cols <- c("#FF9600","#FF76B9","#F68DFF","#FF847A","#00C2CD",
               "#72D2DB","#94B600","#DBA800",
               "#00C982","#00BFFF","#8EA9FF","#DC143C","#0000FF","#20B2AA")
p3 <- DimPlot(sc_combined, group.by = "anno2", reduction = "tsne",
              split.by = "orig.ident",
              cols = anno2cols,
              pt.size = 0.1, label = TRUE,
              raster = FALSE)  +
  xlab(label = "TSNE1") + labs(title = "") + ylab(label = "TSNE2") +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color="black"))
ggsave(p3, file="figures_integration_250410/TSNE_cplorbyAnno2Clusters.pdf",
       height = 4.5, width = 10)

p4 <- ggplot(uni_anno_log2change_df,aes(x=cluster,y=fc,colour=group)) +
  geom_point(size=4) +coord_flip() +ylab(label = "log2(KO/WT)") +xlab(label = "") +
  theme_classic() +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black")) +
  geom_hline(yintercept = 0,color="black") +NoLegend()
ggsave(p4, file="figures_integration_250408/CellRatio.pdf",
       height = 4.5, width = 7)


anno_freq <- data.frame(table(sc_combined$anno2, sc_combined$orig.ident))
colnames(anno_freq) <- c("cluster","sample","number")

anno_freq$percent <- anno_freq$number
anno_freq$percent[anno_freq$sample=="d17_1"] <- 
  anno_freq$percent[anno_freq$sample=="d17_1"] /sum(anno_freq$percent[anno_freq$sample=="d17_1"]) *100
anno_freq$percent[anno_freq$sample=="T-2-D17"] <- 
  anno_freq$percent[anno_freq$sample=="T-2-D17"] /sum(anno_freq$percent[anno_freq$sample=="T-2-D17"]) *100

uni_anno <- unique(anno_freq$cluster)
uni_anno_log2change <- log2(anno_freq$percent[anno_freq$sample=="T-2-D17"]/
                              anno_freq$percent[anno_freq$sample=="d17_1"])


uni_anno_log2change_df <- data.frame(cluster=levels(anno_freq$cluster),
                                     fc = uni_anno_log2change)
uni_anno_log2change_df$group <- rep("up",nrow(uni_anno_log2change_df))
uni_anno_log2change_df$group[uni_anno_log2change_df$fc<0] <- "down"

p4 <- ggplot(uni_anno_log2change_df,aes(x=cluster,y=fc,colour=group)) +
  geom_point(size=4) +coord_flip() +ylab(label = "log2(KO/WT)") +xlab(label = "") +
  theme_classic() +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black")) +
  geom_hline(yintercept = 0,color="black") +NoLegend()
ggsave(p4, file="figures_integration_250410/CellRatio.pdf",
       height = 4.5, width = 7)

anno_freq$log10CellNumber <- log10(anno_freq$number)
anno_freq$log10CellNumber[anno_freq$log10CellNumber=="-Inf"] <- 0
p5 <- ggplot(anno_freq,aes(x=cluster,y=log10CellNumber, color=sample))+
  geom_point(size=3)+
  ylab(label = "log10(Count)") +xlab(label = "") +
  scale_color_manual(values = colors1)+
  theme_classic() +
  theme(axis.text.x = element_text(colour="black",angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(colour="black")) 
ggsave(p5, file="figures_integration_250410/CellCount.pdf",
       height = 4.5, width = 7)

anno2Cluster <- setdiff(levels(sc_combined$anno2), c("Neu","Ery","Mac"))
anno2_kowt_degs_list <- lapply(anno2Cluster,function(x){
  xx <- subset(sc_combined,anno2 %in% x)
  Idents(xx) <- "orig.ident"
  xx_degs <- FindMarkers(xx,ident.1="d17_1",ident.2="T-2-D17",
                         logfc.threshold = 0, min.pct = 0)
  xx_degs$gene <- row.names(xx_degs)
  return(xx_degs)
})
names(anno2_kowt_degs_list) <- anno2Cluster

# endoderm
endo <- subset(sc_combined,anno2 %in% c("Endoderm"))
Idents(endo) <- "orig.ident"
endo_degs <- FindMarkers(endo,ident.1="d17_1",ident.2="T-2-D17")
endo_degs$gene <- row.names(endo_degs)

write.table(endo_degs, file = "figures_integration_250410/T-2-D17_D17-1_endoderm_DEGs.txt", sep="\t", 
            row.names = FALSE, quote=FALSE)

endo_degs_sig <- endo_degs[(endo_degs$p_val_adj<0.05 & abs(endo_degs$avg_log2FC)>=1),]

endo <- ScaleData(endo, features = row.names(endo))

DoHeatmap(endo,features = endo_degs_sig[order(endo_degs_sig$avg_log2FC,decreasing = T),"gene"],
          group.by = "orig.ident")
ggsave(filename = "figures_integration_250410/endoderm_WTvsKO_logfc1_degs_heatmap.pdf",
       height = 14, width = 6)

library(clusterProfiler)
library(org.Hs.eg.db)

degs_list <- list(endo_KO_UP = endo_degs_sig$gene[endo_degs_sig$avg_log2FC>0],
                  endo_KO_DOWN = endo_degs_sig$gene[endo_degs_sig$avg_log2FC<0])

degs_list_hEntrez <- lapply(degs_list,function(x){
  entrez <- unique(na.omit(feature_df$hg19_entrez_id[feature_df$macFas5_gene_symbol %in% x]))
  return(entrez)
})

degs_list_CompareCluster_BP <- compareCluster(geneClusters = degs_list_hEntrez,
                                              fun = "enrichGO",
                                              OrgDb = "org.Hs.eg.db", 
                                              ont ="BP", 
                                              pvalueCutoff = 0.05, 
                                              pAdjustMethod = "BH", 
                                              qvalueCutoff = 0.2,
                                              readable=T)
write.table(degs_list_CompareCluster_BP@compareClusterResult, 
            file="figures_integration_250410/T-2-D17_D17-1_endoderm_DEGS_BP.txt",sep="\t",row.names = FALSE, quote=FALSE)

# nascent mesoderm
nasM <- subset(sc_combined,anno2 %in% c("NM"))
Idents(nasM) <- "orig.ident"
table(nasM$orig.ident)
nasM_degs <- FindMarkers(nasM,ident.1="d17_1",ident.2="T-2-D17")
nasM_degs$gene <- row.names(nasM_degs)
write.table(nasM_degs, file = "figures_integration_250410/T-2-D17_D17-1_nasMes_DEGs.txt", sep="\t", 
            row.names = FALSE, quote=FALSE)

nasM_degs_sig <- nasM_degs[(nasM_degs$p_val_adj<0.05 & abs(nasM_degs$avg_log2FC)>=1),]

nasM <- ScaleData(nasM, features = row.names(nasM))
DoHeatmap(nasM,features = nasM_degs_sig[order(nasM_degs_sig$avg_log2FC,decreasing = T),"gene"],
          group.by = "orig.ident")
ggsave(filename = "figures_integration_250410/nasMes_WTvsKO_logfc1_degs_heatmap.pdf",
       height = 14, width = 6)


degs_list <- list(nasM_KO_UP = nasM_degs_sig$gene[nasM_degs_sig$avg_log2FC>0],
                  nasM_KO_DOWN = nasM_degs_sig$gene[nasM_degs_sig$avg_log2FC<0])
degs_list_hEntrez <- lapply(degs_list,function(x){
  entrez <- unique(na.omit(feature_df$hg19_entrez_id[feature_df$macFas5_gene_symbol %in% x]))
  return(entrez)
})

degs_list_CompareCluster_BP <- compareCluster(geneClusters = degs_list_hEntrez,
                                              fun = "enrichGO",
                                              OrgDb = "org.Hs.eg.db", 
                                              ont ="BP", 
                                              pvalueCutoff = 0.05, 
                                              pAdjustMethod = "BH", 
                                              qvalueCutoff = 0.2,
                                              readable=T)
write.table(degs_list_CompareCluster_BP@compareClusterResult, 
            file="figures_integration_250410/T-2-D17_D17-1_nasMes_DEGS_BP.txt",sep="\t",row.names = FALSE, quote=FALSE)

VlnPlot(sc_combined, group.by = "anno2", split.plot = T, pt.size = 0.1,
        features = "EOMES", split.by = "orig.ident",cols = colors1) + 
  xlab(label = "") +NoLegend()


save.image("/data/embryo_dev/myprocess/T-2-D17/integrated_d17-1/T-2-D17_D17-1_IntegratedAnalysis_250412.RData")




intersect(Reduce(intersect, list(gast_degs$gene[gast_degs$p_val_adj<0.05], nasM_degs$gene[gast_degs$p_val_adj<0.05], endo_degs$gene[endo_degs$p_val_adj<0.05])),
          tGlist1)
tmp1 <- intersect(Reduce(intersect, list(gast_degs$gene[gast_degs$p_val_adj<0.05], nasM_degs$gene[gast_degs$p_val_adj<0.05], endo_degs$gene[endo_degs$p_val_adj<0.05])),
          wntGlist)
tmp2 <- intersect(Reduce(intersect, list(gast_degs$gene[gast_degs$p_val_adj<0.05], nasM_degs$gene[gast_degs$p_val_adj<0.05], endo_degs$gene[endo_degs$p_val_adj<0.05])),
          tgfGlist)



VlnPlot(sc_combined, group.by = "anno2", split.plot = T, pt.size = 0.1, stack = TRUE,flip = T,
        features = unique(c("SOX2","T","MIXL1","EOMES","MESP1","MESP2","WNT5B",
                     "CTNNB1","TCF3","EZH2","TBL1XR1",
                     tmp1,tmp2)), split.by = "orig.ident",cols = colors1) + 
  xlab(label = "") +NoLegend()

ggsave(filename = "selected_KOWT_degs_vlnplot.pdf", height = 8 ,width=12)




HOXGenes <- row.names(sc_combined)[grep("^HOX", row.names(sc_combined))]

FeaturePlot(sc_combined, features = HOXGenes,raster = FALSE, ncol = 5,
            pt.size = 0.01, cols = c("lightgrey","red"))
ggsave(filename="HOXgenes_KOWT_featureplot.pdf",
       height = 15, width = 15)
FeaturePlot(sc_combined, features = HOXGenes,raster = FALSE, ncol = 5, reduction = "tsne",
            pt.size = 0.01, cols = c("lightgrey","red"))
ggsave(filename="HOXgenes_KOWT_featureplot_tsne.pdf",
       height = 15, width = 15)

VlnPlot(sc_combined, group.by = "anno2", split.plot = T, pt.size = 0.1, stack = TRUE,flip = T,
        features = HOXGenes, split.by = "orig.ident",cols = colors1) + 
  xlab(label = "") +NoLegend()
ggsave(filename="HOXgenes_KOWT_vlnplot.pdf",
       height = 15, width = 8)

save.image("/data/embryo_dev/myprocess/T-2-D17/integrated_d17-1/T-2-D17_D17-1_IntegratedAnalysis_250415.RData")


FeaturePlot(sc_combined, 
            features = c("HOXA1","HOXA9","HOXA10","HOXB2","HOXB3","HOXB4","HOXB7","HOXB8","HOXB9"),
            raster = FALSE, 
            split.by = "orig.ident",
            ncol = 5, 
            reduction = "tsne",
            pt.size = 0.01, cols = c("lightgrey","red"))

ggsave(filename="selected_HOXGenes_featurePlots_tsne.pdf",
       width = 5.68, height = 20.78)


FeaturePlot(sc_combined, 
            features = c("CTNNB1","TCF3","VANGL2","PLCB3","SMAD4","MAP3K7"),
            raster = FALSE, 
            split.by = "orig.ident",
            reduction = "tsne",
            pt.size = 0.01, cols = c("lightgrey","red"))
ggsave(filename="selected_WNT_Genes_featurePlots_tsne.pdf",
       width = 5.68, height = 13.78)

FeaturePlot(sc_combined, 
            features = c("ID3","ID1","SMAD7","TGFBR1","SMAD5","THBS2"),
            raster = FALSE, 
            split.by = "orig.ident",
            reduction = "tsne",
            pt.size = 0.01, cols = c("lightgrey","red"))
ggsave(filename="selected_TGFB_Genes_featurePlots_tsne.pdf",
       width = 5.68, height = 13.78)



VlnPlot(sc_combined, group.by = "anno2", split.plot = T, pt.size = 0.1, stack = TRUE,flip = T,
        features = c("HOXA1","HOXA9","HOXA10","HOXB2","HOXB3","HOXB4","HOXB7","HOXB8","HOXB9"), 
        split.by = "orig.ident",cols = colors1) + 
  xlab(label = "") 
ggsave(filename="selectedHOXgenes_KOWT_vlnplot.pdf",
       height = 13, width = 8)

VlnPlot(sc_combined, group.by = "anno2", split.plot = T, pt.size = 0.1, stack = TRUE,flip = T,
        features = c("MESP1","CDH1","EZH2","WNT4","WNT8A","RAC1","CACYBP","SFRP1",
                     "OTX2","OSR1","NODAL","FGF17","FGF8","WNT5A","WNT3","WNT3A","KDR","SNAI1",
                     "MIXL1","FOXC1"), 
        split.by = "orig.ident",cols = colors1) + 
  xlab(label = "") 
ggsave(filename="selectedSignalgenes_250509_KOWT_vlnplot.pdf",
       height = 15, width = 10)

FeaturePlot(sc_combined, 
            features = c("MESP1","CDH1","EZH2","WNT4","WNT8A","RAC1","CACYBP","SFRP1",
                         "OTX2","OSR1","NODAL","FGF17","FGF8","WNT5A","WNT3","WNT3A","KDR","SNAI1",
                         "MIXL1","FOXC1"),
            raster = FALSE,  ncol = 4,
            split.by = "orig.ident",
            reduction = "tsne",
            pt.size = 0.01, cols = c("lightgrey","red"))
ggsave(filename="SelectedFeatures_KOWT_featurePlots_250509_tsne.pdf",
       width = 5.68, height = 40.78)



VlnPlot(sc_combined, group.by = "anno2", split.plot = T, pt.size = 0.1, stack = TRUE,flip = T,
        features = c("TAL1","KDR","KLF1","CDH5","TIE1","ITGA2B","EPOR","LMO2","GATA1","FLT1","DLL4","TEK","PECAM1","GFI1B"), 
        split.by = "orig.ident",cols = colors1) + 
  xlab(label = "") 
ggsave(filename="SelectedFeatures_KOWT_vlnplot_250503.pdf",
       height = 10, width = 8)

FeaturePlot(sc_combined, 
            features = c("TAL1","KDR","KLF1","CDH5","TIE1","ITGA2B","EPOR","LMO2","GATA1","FLT1","DLL4","TEK","PECAM1","GFI1B"),
            raster = FALSE,  ncol = 4,
            split.by = "orig.ident",
            reduction = "tsne",
            pt.size = 0.01, cols = c("lightgrey","red"))
ggsave(filename="SelectedFeatures_KOWT_featurePlots_tsne.pdf",
       width = 5.68, height = 31.78)






dir.create("cellchat")

library(CellChat)
library(patchwork)
library(BiocParallel)

# T-KO-D17
alldata <- sc_combined@assays$RNA@data
allmeta <- sc_combined@meta.data

# transform monkey gene to human gene
alldata_humangenes <- feature_df$hg19_gene_symbol[match(row.names(alldata),feature_df$macFas5_gene_symbol)]

alldata <- alldata[(!is.na(alldata_humangenes)),]
row.names(alldata) <- alldata_humangenes[(!is.na(alldata_humangenes))]

d17_cells <- row.names(allmeta)[allmeta$orig.ident=="d17_1"]
td17_cells <- row.names(allmeta)[allmeta$orig.ident=="T-2-D17"]

# D17_1
data.input <- alldata[,d17_cells]
meta <- allmeta[d17_cells,]
meta <- meta[,c("orig.ident","nFeature_RNA","anno2")]
colnames(meta) <- c("sample","geneNumber","labels")
unique(meta$labels) 
meta$labels = droplevels(meta$labels, exclude = setdiff(levels(meta$labels),unique(meta$labels)))
unique(meta$labels)

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
groupSize

CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB

cellchat <- subsetData(cellchat) 

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

saveRDS(cellchat,"cellchat/D17_1_Cellchat.rds")


# T-KO-D17
data.input <- alldata[,td17_cells]
meta <- allmeta[td17_cells,]
meta <- meta[,c("orig.ident","nFeature_RNA","anno2")]
colnames(meta) <- c("sample","geneNumber","labels")
unique(meta$labels) 
meta$labels = droplevels(meta$labels, exclude = setdiff(levels(meta$labels),unique(meta$labels)))
unique(meta$labels)

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
groupSize

CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB

cellchat <- subsetData(cellchat) 

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

saveRDS(cellchat,"cellchat/T-KO-D17_Cellchat.rds")

rm(cellchat)



cellchat.NTC <- readRDS("cellchat/D17_1_Cellchat.rds")
cellchat.KO <- readRDS("cellchat/T-KO-D17_Cellchat.rds")

object.list <- list(WT = cellchat.NTC, KO = cellchat.KO)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat

### Part I: Predict general principles of cell-cell communication
## Compare the total number of interactions and interaction strength
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

###Compare the number of interactions and interaction strength among different cell populations
#To identify the interaction between which cell populations showing significant changes, CellChat compares the number of interactions and interaction strength among different cell populations.

##Differential number of interactions or interaction strength among different cell populations
#The differential number of interactions or interaction strength in the cell-cell communication network between two datasets can be visualized using circle plot, where red
#(or blue) colored edges represent increased (or decreased) signaling in the second dataset compared to the first one.

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")


#We can also show differential number of interactions or interaction strength in a greater details using a heatmap. The top colored bar plot represents the sum of column of values displayed in the heatmap (incoming signaling). The right colored bar plot represents the sum of row of values (outgoing signaling). In the colorbar, red
#(or blue) represents increased(or decreased) signaling in the second dataset compared to the first one.
gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2


#The differential network analysis only works for pairwise datasets. If there are more datasets for comparison, we can directly show the number of interactions or interaction strength between any two cell populations in each dataset.
#To better control the node size and edge weights of the inferred networks across different datasets, we compute the maximum number of cells per cell group and the maximum number of interactions (or interaction weights) across all datasets.

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

### Compare the major sources and targets in 2D space
## Comparing the outgoing and incoming interaction strength in 2D space allows ready identification of the cell populations with significant changes in sending or receiving signals between different datasets.

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)



#Identify and visualize the conserved and context-specific signaling pathways
#By comparing the information flow/interaction strengh of each signaling pathway, we can identify signaling pathways, (i) turn off, (ii) decrease, (iii) turn on or (iv) increase, by change their information flow at one condition as compared to another condition.
#Compare the overall information flow of each signaling pathway
#We can identify the conserved and context-specific signaling pathways by simply comparing the information flow for each signaling pathway, which is defined by the sum of communication probability among all pairs of cell groups in the inferred network (i.e., the total weights in the network).
#This bar graph can be plotted in a stacked mode or not. Significant signaling pathways were ranked based on differences in the overall information flow within the inferred networks between NL and LS skin. The top signaling pathways colored red are enriched in NL skin, and these colored green were enriched in the LS skin.
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2


library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 8 ,height = 15)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 8, height = 15)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))


ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 8, height = 15, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 8, height = 15, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))


ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 8, height = 15, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 8, height = 15, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))




pathways.show <- c("FGF") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], 
                      edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}


pathways.show <- c("FGF") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, 
                               color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
#> Do heatmap based on a single object 
#> 
#> Do heatmap based on a single object
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))




cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("WT", "KO")) # set factor level
plotGeneExpression(cellchat, signaling = "FGF", split.by = "datasets", colors.ggplot = T,split.plot = T)
plotGeneExpression(cellchat, signaling = "TGFb", split.by = "datasets", colors.ggplot = T,split.plot = T)
plotGeneExpression(cellchat, signaling = "WNT", split.by = "datasets", colors.ggplot = T,split.plot = T)
