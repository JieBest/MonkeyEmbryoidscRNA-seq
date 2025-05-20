################################################################################################
###This script integrate human embryo PCW3(CS10) dataset (Zeng et al., 2023) and Monkey embryoid
################################################################################################




library(Seurat)
library(ggplot2)
library(SCpubr)
library(harmony)
library(cowplot)


features_df  <- read.csv("../../public_data/features_hg19.tsv",sep="\t",header=T)

sc_list <- readRDS("../../cr6ncbi_d17tod28_nFatureRNA1000_Singlet_nCountRNA30000.rds")
sc_list <- subset(sc_list, orig.ident %in% 
                    c("d17_1_cr6","d18_1_cr6","d21_1_cr6","d22_1_cr6",
                      "d23_1","d25_1_cr6","d25_2_cr6"))
new <- sc_list
new_counts <- new@assays$RNA@counts
new_meta <- new@meta.data

# change mfas gene to human gene
features <- features_df$hg19_gene_symbol[match(row.names(new_counts), features_df$macFas5_gene_symbol)]
index <- which(is.na(features)==FALSE)
sub_counts <- new_counts[index,]
row.names(sub_counts) <- features_df$hg19_gene_symbol[match(row.names(sub_counts), features_df$macFas5_gene_symbol)]
# sub_counts 16653 X 93596
new_counts <- sub_counts # 16653 X 93596


pcw3 <- load("../../public_data/pcw3_filtered_12323cells_seuobj.RData")
pcw3_counts <-  pcw3_filtered@assays$RNA@counts # 21801 X 41359
pcw3_meta <- pcw3_filtered@meta.data

sharedGs <- intersect(row.names(new_counts), row.names(pcw3_counts))

new_counts <- new_counts[sharedGs,]
pcw3_counts <- pcw3_counts[sharedGs,]


new_meta_sub <- new_meta[,c(1,2,3,9)]
colnames(new_meta_sub) <- c("sample","nCount_RNA","nFeature_RNA","annotation")
new_meta_sub$study <- rep("blastoid",nrow(new_meta_sub))

pcw3_meta_sub <- pcw3_meta[,c(1,2,3,6)]
colnames(pcw3_meta_sub) <- c("sample","nCount_RNA","nFeature_RNA","annotation")
pcw3_meta_sub$study <- rep("human",nrow(pcw3_meta_sub))


blastoid <- CreateSeuratObject(counts = new_counts,meta.data = new_meta_sub)


pcw3 <- CreateSeuratObject(counts = pcw3_counts, meta.data = pcw3_meta_sub) 
pcw3$orig.ident <- pcw3$sample
Idents(pcw3) <- 'orig.ident'

sc_list =list(blastoid=blastoid, pcw3=pcw3)
saveRDS(sc_list, file = "raw_mfas_blastoidsD17-1D25-2_humanPCW3_sc_list.RDS")


rm(list=ls())
gc()

sc_list <- readRDS("raw_mfas_blastoidsD17-1D25-2_humanPCW3_sc_list.RDS")

# normalize and identify variable features for each dataset independently
sc_list <- lapply(X = sc_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})
# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = sc_list, nfeatures = 3000)
sc_anchors <- FindIntegrationAnchors(object.list = sc_list, anchor.features = features)
saveRDS(sc_anchors, file = "mfas_blastoidsD17-1D25-2_humanPCW3_sc_anchors.RDS")

rm(sc_list);gc()
sc_anchors <- readRDS("mfas_blastoidsD17-1D25-2_humanPCW3_sc_anchors.RDS")
#sc_anchors <- readRDS("mfas-human_sc_anchors.RDS")
# this command creates an 'integrated' data assay
sc_combined <- IntegrateData(anchorset = sc_anchors)
rm(sc_anchors)
gc()

saveRDS(sc_combined, file = "mfas_blastoidsD17-1D25-2_humanPCW3_sc_combined0.RDS")


setwd("/data/embryo_dev/myprocess/integrate_cr6ncbi_d17-1_d25-2_nod22-2/blastoidD17-1_D25_2_humanpcw3/")

load("mfas_blastoidsD17-1D25-2_humanPCW3_sc_combined0.RData")
alldata_meta <- alldata@meta.data

load("/data/embryo_dev/myprocess/integrate_cr6ncbi_d17-1_d25-2_nod22-2/vstCCA/pc20/pc20_d17-1_d25-2_nod22-2_seuobj_withClusterInfo_meta.RData")
#sc_combined_meta$integrated_snn_res.1 <- paste0("blastoid_", sc_combined_meta$integrated_snn_res.1)
table(sc_combined_meta$integrated_snn_res.1 )

alldata_meta$annotation[match(row.names(sc_combined_meta), row.names(alldata_meta))] <- as.character(sc_combined_meta$integrated_snn_res.1 )
alldata_meta$annotation <- paste0(alldata_meta$study,"_", alldata_meta$annotation)
table(alldata_meta$annotation)

alldata$annotation <- alldata_meta$annotation

rm(alldata_meta)


# Scale the integrated data, perform PCA, and generate an elbow plot for dimensionality estimation
alldata <- ScaleData(alldata)
alldata <- RunPCA(object = alldata, npcs = 100)
ElbowPlot(alldata, ndims = 100)
ggsave(filename = "blastoidD17-1_D25-2_humanpcw3_ElbowPlot.png", width = 7.68, height = 3.98, bg = "white")

# UMAP Plotting and Clustering
# -----------------------------
# Setting parameters for UMAP analysis
pcause = 100  # Number of dimensions to consider before analysis
alldata <- RunUMAP(object = alldata, dims = 1:pcause)
#alldata <- RunTSNE(object = alldata, dims = 1:pcause)
alldata <- FindNeighbors(object = alldata, dims = 1:pcause)

alldata <- FindClusters(object = alldata, resolution = 0.4)

toplist <- names(alldata$study[alldata$study == "blastoid"])

p1 <- DimPlot(object = alldata, reduction = "umap", order = toplist, 
              raster = F, group.by = "study", pt.size = 0.0001)
p2 <- DimPlot(object = alldata, reduction = "umap", label = T,
              raster = F, group.by = "integrated_snn_res.0.4", pt.size = 0.0001)
pp <- p1+p2
ggsave(pp, filename = "blastoidD17-1_D25-2_humanpcw3_noharmony_res04.umap.png",
       width = 12.69, height=4.89)


save(alldata, file="blastoidD17-1_D25-2_humanpcw3_noharmony_withCluster.RData")

# Harmony Integration
# -------------------
# Running Harmony algorithm for data integration across different conditions
pdf(file = "alldata_harmony_use.pdf", width = 8, height = 6)
alldata_harmony <- RunHarmony(alldata, group.by.vars = "study", plot_convergence = TRUE)
dev.off()

rm(alldata);gc()

# Plotting Harmony integrated data and saving the output

DimPlot(object = alldata_harmony, raster = FALSE,  reduction = "harmony", group.by = "study", pt.size = 0.01)
ggsave(filename = "alldata_harmony_cluster.png", width = 8, height = 6)

# Further UMAP analysis on Harmony integrated data
alldata_harmony <- alldata_harmony %>% RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution = 0.5) %>%
  identity()


# Saving UMAP plots of Harmony integrated data
DimPlot(object = alldata_harmony, reduction = "umap", raster = F, label = TRUE, pt.size = 0.01)
ggsave(filename = "alldata_harmony_cluster_colByRes05cluster.png", width = 8, height = 6)


DimPlot(object = alldata_harmony, reduction = "umap",order = toplist, 
        raster = F, group.by = "study", pt.size = 0.01)
ggsave(filename = "alldata_harmony_cluster_colByStudy.png", width = 8, height = 6)

alldata_harmony <- FindClusters(alldata_harmony, resolution = 0.6) %>% identity()
DimPlot(object = alldata_harmony, reduction = "umap", 
        raster = F, label = TRUE, pt.size = 0.01)
ggsave(filename = "alldata_harmony_cluster_colByRes06cluster.png", width = 8, height = 6)

alldata_harmony <- FindClusters(alldata_harmony, resolution = 0.7) %>% identity()
DimPlot(object = alldata_harmony, reduction = "umap", 
        raster = F, label = TRUE, pt.size = 0.01)
ggsave(filename = "alldata_harmony_cluster_colByRes07cluster.png", width = 8, height = 6)

alldata_harmony <- FindClusters(alldata_harmony, resolution = 0.8) %>% identity()
DimPlot(object = alldata_harmony, reduction = "umap", 
        raster = F, label = TRUE, pt.size = 0.01)
ggsave(filename = "alldata_harmony_cluster_colByRes08cluster.png", width = 8, height = 6)


alldata_harmony <- FindClusters(alldata_harmony,resolution = 0.9) %>% identity()
DimPlot(object = alldata_harmony, reduction = "umap", 
        raster = F, label = TRUE, pt.size = 0.01)
ggsave(filename = "alldata_harmony_cluster_colByRes09cluster.png", width = 8, height = 6)

alldata_harmony <- FindClusters(alldata_harmony,resolution = 1) %>% identity()
DimPlot(object = alldata_harmony, reduction = "umap", 
        raster = F, label = TRUE, pt.size = 0.01)
ggsave(filename = "alldata_harmony_cluster_colByRes1cluster.png", width = 8, height = 6)


options(max.print = 10000)
table(alldata_harmony$annotation,alldata_harmony$integrated_snn_res.0.5)

save(alldata_harmony, file="mfas_blastoidD17-1_D25-2_humanpcw3_alldata_harmony.RData")


alldata_harmony_meta <- alldata_harmony@meta.data
save(alldata_harmony_meta, file = "mfas_blastoidD17-1_D25-2_humanpcw3_alldata_harmony_metainfo.RData")

library(SCpubr)
studyClusters <- unique(alldata_harmony$annotation)
for(i in 1:length(studyClusters)){
  pp <- do_DimPlot(alldata_harmony, raster = F, group.by = "annotation",idents.keep = studyClusters[i],pt.size=0.001)
  ggsave(pp ,filename=paste0("mfas_blastoidD17-1_D25-2_humanpcw3_alldata_harmony_",studyClusters[i],"_umap.png"), width = 8, height = 6)
}

res05Clusters <- unique(alldata_harmony$integrated_snn_res.0.5)
for(i in 1:length(res05Clusters)){
  pp <- do_DimPlot(alldata_harmony, raster = F, group.by = "integrated_snn_res.0.5",idents.keep = res05Clusters[i],pt.size=0.001)
  ggsave(pp ,filename=paste0("mfas_blastoidD17-1_D25-2_humanpcw3_alldata_harmony_res05Cluster",res05Clusters[i],"_umap.png"), width = 8, height = 6)
}

save.image("blastoidD17-1_D25-2_humanpcw3_integrated.RData")

# Adjusting metadata for Harmony data

alldata_harmony@meta.data$study <- factor(alldata_harmony@meta.data$study, levels = c("blastoid", "human"))
toplist <- names(alldata_harmony$study[alldata_harmony$study == "blastoid"])

# Customized UMAP plots with specific color schemes
# Changing colors of the clusters and generating multiple UMAP plots
pp <- DimPlot(object = alldata_harmony, reduction = "umap", cols = c("#2979A9", "#C04844"), raster = FALSE, 
              order = toplist, group.by = "study", pt.size = 0.00001) +labs(title = "")
pp[[1]]$layers[[1]]$mapping$alpha= 0.7
ggsave(pp, filename = "alldata_harmony_cluster_changecolor_threecolor1.pdf", width = 8, height = 6)
ggsave(pp, filename = "alldata_harmony_cluster_changecolor_threecolor1.png", width = 8, height = 6)

toplist <- names(alldata_harmony$study[alldata_harmony$study == "human"])
pp <- DimPlot(object = alldata_harmony, reduction = "umap", cols = c("#2979A9", "#C04844"), raster = FALSE, 
              order = toplist, group.by = "study", pt.size = 0.00001) +labs(title = "")
pp[[1]]$layers[[1]]$mapping$alpha= 0.7
ggsave(pp, filename = "alldata_harmony_cluster_changecolor_threecolor2.pdf", width = 8, height = 6)
ggsave(pp, filename = "alldata_harmony_cluster_changecolor_threecolor2.png", width = 8, height = 6)
