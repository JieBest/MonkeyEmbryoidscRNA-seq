################# 1. load required packages
#################
#################
#################
#################
#################
library(Seurat)
library(dplyr)
library(ggplot2)
library(rlist)
library(cowplot)
library(Stringr)
library(DoubletFinder)
library(clustree)
library(reshape2)
library(SCpubr)
library(ggpubr)
library(RColorBrewer)


################# 2. Load expression matrix data and creat seurat object
#################
#################
#################
#################
#################

# read in data
d17 <- Read10X("/data/embryo_dev/mfas_embryoid_cr6.1.2_ncbi/d17_1/embryoid_D17_1/outs/filtered_feature_bc_matrix/")
d18 <- Read10X("/data/embryo_dev/mfas_embryoid_cr6.1.2_ncbi/d18_1/embryoid_D18_1/outs/filtered_feature_bc_matrix/")
d21 <- Read10X("/data/embryo_dev/mfas_embryoid_cr6.1.2_ncbi/d21_1/embryoid_D21_1/outs/filtered_feature_bc_matrix/")
d22 <- Read10X("/data/embryo_dev/mfas_embryoid_cr6.1.2_ncbi/d22_1/embryoid_D22_1/outs/filtered_feature_bc_matrix/")
d23 <- Read10X("/data/embryo_dev/mfas_embryoid_cr6.1.2_ncbi/d23_1/embryoid_D23_1/outs/filtered_feature_bc_matrix/")
d25_1 <- Read10X("/data/embryo_dev/mfas_embryoid_cr6.1.2_ncbi/d25_1/embryoid_D25_1/outs/filtered_feature_bc_matrix/")
d25_2 <- Read10X("/data/embryo_dev/mfas_embryoid_cr6.1.2_ncbi/d25_2/embryoid_D25_2/outs/filtered_feature_bc_matrix/")

## Create the Seurat objects
d17 <- CreateSeuratObject(counts = d17,  min.cells = 3, min.features = 200, project = "D17")
d18 <- CreateSeuratObject(counts = d18,  min.cells = 3, min.features = 200, project = "D18")
d21 <- CreateSeuratObject(counts = d21,  min.cells = 3, min.features = 200, project = "D21")
d22 <- CreateSeuratObject(counts = d22,  min.cells = 3, min.features = 200, project = "D22")
d23 <- CreateSeuratObject(counts = d23,  min.cells = 3, min.features = 200, project = "D23")
d25_1 <- CreateSeuratObject(counts = d25_1,  min.cells = 3, min.features = 200, project = "D25_1")
d25_2 <- CreateSeuratObject(counts = d25_2,  min.cells = 3, min.features = 200, project = "D25_2")


################# 2.  QC for each sample 
#################
#################
#################
#################
#################
# Identify apoptotic cells by % Mitochondrial gene expression (NOTE that I do not use this to directly fiter out cells)
Mitochondrial_genes<- c("ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB", "ND1","ND2", "ND3", "ND4", "ND4L", "ND5", "ND6")
{
  d17[["percent.mt"]] <- PercentageFeatureSet(d17,features = Mitochondrial_genes)
  d18[["percent.mt"]] <- PercentageFeatureSet(d18,features = Mitochondrial_genes)
  d21[["percent.mt"]] <- PercentageFeatureSet(d21,features = Mitochondrial_genes)
  d22[["percent.mt"]] <- PercentageFeatureSet(d22,features = Mitochondrial_genes)
  d23[["percent.mt"]] <- PercentageFeatureSet(d23,features = Mitochondrial_genes)
  d25_1[["percent.mt"]] <- PercentageFeatureSet(d25_1,features = Mitochondrial_genes)
  d25_2[["percent.mt"]] <- PercentageFeatureSet(d25_2,features = Mitochondrial_genes)
}

### QUALITY CONTROL PLOTS ###
# Now that we can visualise the distibutions we can remove cells that have a very high feature count (potential doublets)
# or remove apoptotic cells which have high mitochondrial gene expression
{
  VlnPlot(d17, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
  VlnPlot(d18, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
  VlnPlot(d21, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
  VlnPlot(d22, features = c("nFeature_RNA", "nCount_RNA", "percent.mt")) 
  VlnPlot(d23, features = c("nFeature_RNA", "nCount_RNA", "percent.mt")) 
  VlnPlot(d25_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt")) 
  VlnPlot(d25_2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt")) 
}

#Identify and remove doublets with doublet finder
pc.num <- 1:30 # using top 30 PCs

# d17 
obj <- d17
obj <- obj %>% SCTransform() %>% RunPCA() %>% RunUMAP(dim = pc.num) %>% FindNeighbors() %>% FindClusters()
sweep.res.list=paramSweep_v3(obj,PCs = pc.num,sct = T) ## sct:Logical representing whether SCTransform was used during original seurat object pre-processing
## ParamSweep_ V3 is a function that calculates the combination of pN and pK parameters
sweep.stats=summarizeSweep(sweep.res.list,GT=F)
opt_pK=find.pK(sweep.stats)
mpK <- as.numeric(as.vector(opt_pK$pK[which.max(opt_pK$BCmetric)])) ##select the best pK
homotypic.prop=modelHomotypic(annotations = obj@meta.data$seurat_clusters)
DoubletRate = 0.075  ## Assuming 7.5% doublet formation rate 
nExp_poi <- round(DoubletRate*nrow(obj@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
obj <- doubletFinder_v3(obj, PCs = pc.num, pN = 0.25, pK = mpK,
                        nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
# plot
p1 <- DimPlot(obj, reduction = "umap",group.by = "SCT_snn_res.0.8")
p2 <- DimPlot(obj, reduction = "umap", group.by = "DF.classifications_0.25_0.03_961")
p1|p2
d17$pANN <- obj$pANN_0.25_0.03_961
d17$DF <- obj$DF.classifications_0.25_0.03_961

# d18 
obj <- d18
obj <- obj %>% SCTransform() %>% RunPCA() %>% RunUMAP(dim = pc.num) %>% FindNeighbors() %>% FindClusters()
sweep.res.list=paramSweep_v3(obj,PCs = pc.num,sct = T) ## sct:Logical representing whether SCTransform was used during original seurat object pre-processing
## ParamSweep_ V3 is a function that calculates the combination of pN and pK parameters
sweep.stats=summarizeSweep(sweep.res.list,GT=F)
opt_pK=find.pK(sweep.stats)
mpK<-as.numeric(as.vector(opt_pK$pK[which.max(opt_pK$BCmetric)])) ##select the best pK
homotypic.prop=modelHomotypic(annotations = obj@meta.data$seurat_clusters)
DoubletRate = 0.075  ## Assuming 7.5% doublet formation rate 
nExp_poi <- round(DoubletRate*nrow(obj@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
obj <- doubletFinder_v3(obj, PCs = pc.num, pN = 0.25, pK = mpK,
                        nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
# plot
p1 <- DimPlot(obj, reduction = "umap")
p2 <- DimPlot(obj, reduction = "umap", group.by = "DF.classifications_0.25_0.28_1274")
p1|p2
d18$pANN<- obj$pANN_0.25_0.28_1274
d18$DF <- obj$DF.classifications_0.25_0.28_1274

# d21
obj <- d21
obj <- obj %>% SCTransform() %>% RunPCA() %>% RunUMAP(dim = pc.num) %>% FindNeighbors() %>% FindClusters()
sweep.res.list=paramSweep_v3(obj,PCs = pc.num,sct = T) ## sct:Logical representing whether SCTransform was used during original seurat object pre-processing
## ParamSweep_ V3 is a function that calculates the combination of pN and pK parameters
sweep.stats=summarizeSweep(sweep.res.list,GT=F)
opt_pK=find.pK(sweep.stats)
mpK<-as.numeric(as.vector(opt_pK$pK[which.max(opt_pK$BCmetric)])) ##select the best pK
homotypic.prop=modelHomotypic(annotations = obj@meta.data$seurat_clusters)
DoubletRate = 0.075  ## Assuming 7.5% doublet formation rate 
nExp_poi <- round(DoubletRate*nrow(obj@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
obj <- doubletFinder_v3(obj, PCs = pc.num, pN = 0.25, pK = mpK,
                        nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
# plot
p1 <- DimPlot(obj, reduction = "umap",group.by = "SCT_snn_res.0.8",label = TRUE)
p2 <- DimPlot(obj, reduction = "umap", group.by = "DF.classifications_0.25_0.07_1520")
p1|p2
d21$pANN <- obj$pANN_0.25_0.07_1520
d21$DF <- obj$DF.classifications_0.25_0.07_1520

# d22
obj <- d22
obj <- obj %>% SCTransform() %>% RunPCA() %>% RunUMAP(dim = pc.num) %>% FindNeighbors() %>% FindClusters()
sweep.res.list=paramSweep_v3(obj,PCs = pc.num,sct = T) ## sct:Logical representing whether SCTransform was used during original seurat object pre-processing
## ParamSweep_ V3 is a function that calculates the combination of pN and pK parameters
sweep.stats=summarizeSweep(sweep.res.list,GT=F)
opt_pK=find.pK(sweep.stats)
mpK<-as.numeric(as.vector(opt_pK$pK[which.max(opt_pK$BCmetric)])) ##select the best pK
homotypic.prop=modelHomotypic(annotations = obj@meta.data$seurat_clusters)
DoubletRate = 0.075  ## Assuming 7.5% doublet formation rate 
nExp_poi <- round(DoubletRate*nrow(obj@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
obj <- doubletFinder_v3(obj, PCs = pc.num, pN = 0.25, pK = mpK,
                        nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
# plot
p1 <- DimPlot(obj, reduction = "umap",group.by = "SCT_snn_res.0.8",label = TRUE)
p2 <- DimPlot(obj, reduction = "umap", group.by = "DF.classifications_0.25_0.3_1272")
p1|p2
d22$pANN <- obj$pANN_0.25_0.3_1272
d22$DF <- obj$DF.classifications_0.25_0.3_1272

# d23
obj <- d23
obj <- obj %>% SCTransform() %>% RunPCA() %>% RunUMAP(dim = pc.num) %>% FindNeighbors() %>% FindClusters()
sweep.res.list=paramSweep_v3(obj,PCs = pc.num,sct = T) ## sct:Logical representing whether SCTransform was used during original seurat object pre-processing
## ParamSweep_ V3 is a function that calculates the combination of pN and pK parameters
sweep.stats=summarizeSweep(sweep.res.list,GT=F)
opt_pK=find.pK(sweep.stats)
mpK<-as.numeric(as.vector(opt_pK$pK[which.max(opt_pK$BCmetric)])) ##select the best pK
homotypic.prop=modelHomotypic(annotations = obj@meta.data$seurat_clusters)
DoubletRate = 0.075  ## Assuming 7.5% doublet formation rate 
nExp_poi <- round(DoubletRate*nrow(obj@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
obj <- doubletFinder_v3(obj, PCs = pc.num, pN = 0.25, pK = mpK,
                        nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
# plot
p1 <- DimPlot(obj, reduction = "umap",group.by = "SCT_snn_res.0.8",label = TRUE)
p2 <- DimPlot(obj, reduction = "umap", group.by = "DF.classifications_0.25_0.22_1249")
p1|p2
d23$pANN <- obj$pANN_0.25_0.22_1249
d23$DF <- obj$DF.classifications_0.25_0.22_1249

# d25_1
obj <- d25_1
obj <- obj %>% SCTransform() %>% RunPCA() %>% RunUMAP(dim = pc.num) %>% FindNeighbors() %>% FindClusters()
sweep.res.list=paramSweep_v3(obj,PCs = pc.num,sct = T) ## sct:Logical representing whether SCTransform was used during original seurat object pre-processing
## ParamSweep_ V3 is a function that calculates the combination of pN and pK parameters
sweep.stats=summarizeSweep(sweep.res.list,GT=F)
opt_pK=find.pK(sweep.stats)
mpK<-as.numeric(as.vector(opt_pK$pK[which.max(opt_pK$BCmetric)])) ##select the best pK
homotypic.prop=modelHomotypic(annotations = obj@meta.data$seurat_clusters)
DoubletRate = 0.075  ## Assuming 7.5% doublet formation rate 
nExp_poi <- round(DoubletRate*nrow(obj@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
obj <- doubletFinder_v3(obj, PCs = pc.num, pN = 0.25, pK = mpK,
                        nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
# plot
p1 <- DimPlot(obj, reduction = "umap",group.by = "SCT_snn_res.0.8",label = TRUE)
p2 <- DimPlot(obj, reduction = "umap", group.by = "DF.classifications_0.25_0.17_1603")
p1|p2
d25_1$pANN <- obj$pANN_0.25_0.17_1603
d2_1$DF <- obj$DF.classifications_0.25_0.17_1603

# d25_2
obj <- d25_2
obj <- obj %>% SCTransform() %>% RunPCA() %>% RunUMAP(dim = pc.num) %>% FindNeighbors() %>% FindClusters()
sweep.res.list=paramSweep_v3(obj,PCs = pc.num,sct = T) ## sct:Logical representing whether SCTransform was used during original seurat object pre-processing
## ParamSweep_ V3 is a function that calculates the combination of pN and pK parameters
sweep.stats=summarizeSweep(sweep.res.list,GT=F)
opt_pK=find.pK(sweep.stats)
mpK<-as.numeric(as.vector(opt_pK$pK[which.max(opt_pK$BCmetric)])) ##select the best pK
homotypic.prop=modelHomotypic(annotations = obj@meta.data$seurat_clusters)
DoubletRate = 0.075  ## Assuming 7.5% doublet formation rate 
nExp_poi <- round(DoubletRate*nrow(obj@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
obj <- doubletFinder_v3(obj, PCs = pc.num, pN = 0.25, pK = mpK,
                        nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
# plot
p1 <- DimPlot(obj, reduction = "umap",group.by = "SCT_snn_res.0.8",label = TRUE)
p2 <- DimPlot(obj, reduction = "umap", group.by = "DF.classifications_0.25_0.08_1525")
p1|p2
d25_2$pANN <- obj$pANN_0.25_0.08_1525
d2_2$DF <- obj$DF.classifications_0.25_0.08_1525


### calculate cell cycle socre ###
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
s.genes[s.genes=="MLF1IP"] <- "CENPU"

# d17
obj <- d17
obj <- NormalizeData(obj) %>% FindVariableFeatures(selection.method = "vst",nfeatures=3000) %>% ScaleData(features = row.names(obj))
obj <- CellCycleScoring(obj,s.features = s.genes, g2m.features = g2m.genes,set.ident = TRUE)
RidgePlot(obj, features = c("PCNA","TOP2A","MCM6","MKI67"),ncol=2)
d17$S.Score <- obj$S.Score
d17$G2M.Score <- obj$G2M.Score
d17$Phase <- obj$Phase

# d18
obj <- d18
obj <- NormalizeData(obj) %>% FindVariableFeatures(selection.method = "vst",nfeatures=3000) %>% ScaleData(features = row.names(obj))
obj <- CellCycleScoring(obj,s.features = s.genes, g2m.features = g2m.genes,set.ident = TRUE)
RidgePlot(obj, features = c("PCNA","TOP2A","MCM6","MKI67"),ncol=2)
d18$S.Score <- obj$S.Score
d18$G2M.Score <- obj$G2M.Score
d18$Phase <- obj$Phase

# d21
obj <- d21
obj <- NormalizeData(obj) %>% FindVariableFeatures(selection.method = "vst",nfeatures=3000) %>% ScaleData(features = row.names(obj))
obj <- CellCycleScoring(obj,s.features = s.genes, g2m.features = g2m.genes,set.ident = TRUE)
RidgePlot(obj, features = c("PCNA","TOP2A","MCM6","MKI67"),ncol=2)
d21$S.Score <- obj$S.Score
d21$G2M.Score <- obj$G2M.Score
d21$Phase <- obj$Phase

# d22
obj <- d22
obj <- NormalizeData(obj) %>% FindVariableFeatures(selection.method = "vst",nfeatures=3000) %>% ScaleData(features = row.names(obj))
obj <- CellCycleScoring(obj,s.features = s.genes, g2m.features = g2m.genes,set.ident = TRUE)
RidgePlot(obj, features = c("PCNA","TOP2A","MCM6","MKI67"),ncol=2)
d22$S.Score <- obj$S.Score
d22$G2M.Score <- obj$G2M.Score
d22$Phase <- obj$Phase

# d23
obj <- d23
obj <- NormalizeData(obj) %>% FindVariableFeatures(selection.method = "vst",nfeatures=3000) %>% ScaleData(features = row.names(obj))
obj <- CellCycleScoring(obj,s.features = s.genes, g2m.features = g2m.genes,set.ident = TRUE)
RidgePlot(obj, features = c("PCNA","TOP2A","MCM6","MKI67"),ncol=2)
d23$S.Score <- obj$S.Score
d23$G2M.Score <- obj$G2M.Score
d23$Phase <- obj$Phase

# d25_1
obj <- d25_1
obj <- NormalizeData(obj) %>% FindVariableFeatures(selection.method = "vst",nfeatures=3000) %>% ScaleData(features = row.names(obj))
obj <- CellCycleScoring(obj,s.features = s.genes, g2m.features = g2m.genes,set.ident = TRUE)
RidgePlot(obj, features = c("PCNA","TOP2A","MCM6","MKI67"),ncol=2)
d25_1$S.Score <- obj$S.Score
d25_1$G2M.Score <- obj$G2M.Score
d25_1$Phase <- obj$Phase

# d25_2
obj <- d25_2
obj <- NormalizeData(obj) %>% FindVariableFeatures(selection.method = "vst",nfeatures=3000) %>% ScaleData(features = row.names(obj))
obj <- CellCycleScoring(obj,s.features = s.genes, g2m.features = g2m.genes,set.ident = TRUE)
RidgePlot(obj, features = c("PCNA","TOP2A","MCM6","MKI67"),ncol=2)
d25_2$S.Score <- obj$S.Score
d25_2$G2M.Score <- obj$G2M.Score
d25_2$Phase <- obj$Phase

# save in to file
saveRDS(d17, file="/data/embryo_dev/myprocess/rawrds/d17_withDF_Phase.rds")
saveRDS(d18, file="/data/embryo_dev/myprocess/rawrds/d18_withDF_Phase.rds")
saveRDS(d21, file="/data/embryo_dev/myprocess/rawrds/d21_withDF_Phase.rds")
saveRDS(d22, file="/data/embryo_dev/myprocess/rawrds/d22_withDF_Phase.rds")
saveRDS(d23, file="/data/embryo_dev/myprocess/rawrds/d23_withDF_Phase.rds")
saveRDS(d25_1, file="/data/embryo_dev/myprocess/rawrds/d25_1_withDF_Phase.rds")
saveRDS(d25_2, file="/data/embryo_dev/myprocess/rawrds/d25_2_withDF_Phase.rds")


### Filter cells for each sample ###
d17 <- subset(d17, nFeature_RNA > 1000 & DF=="Siglet" & nCount_RNA < 30000)
d18 <- subset(d18, nFeature_RNA > 1000 & DF=="Siglet" & nCount_RNA < 30000)
d21 <- subset(d21, nFeature_RNA > 1000 & DF=="Siglet" & nCount_RNA < 30000)
d22 <- subset(d22, nFeature_RNA > 1000 & DF=="Siglet" & nCount_RNA < 30000)
d23 <- subset(d23, nFeature_RNA > 1000 & DF=="Siglet" & nCount_RNA < 30000)
d25_1 <- subset(d25_1, nFeature_RNA > 1000 & DF=="Siglet" & nCount_RNA < 30000)
d25_2 <- subset(d25_2, nFeature_RNA > 1000 & DF=="Siglet" & nCount_RNA < 30000)

sc_list <- list(d17=d17, d18=d18, d21=d21, d22=d22, d23=d23, d25_1=d25_1, d25_2=d25_2)
saveRDS(sc_list, file="/data/embryo_dev/myprocess/filteredrds/sc_list0.rds")

#################3. Integrate data
#################
#################
#################
#################
#################

# normalize and identify variable features for each dataset independently
sc_list <- lapply(X = sc_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})

features <- SelectIntegrationFeatures(object.list = sc_list, nfeatures = 3000)
sc_anchors <- FindIntegrationAnchors(object.list = sc_list, anchor.features = features)
sc_combined <- IntegrateData(anchorset = sc_anchors)
saveRDS(sc_combined, file="/data/embryo_dev/myprocess/filteredrds/sc_combined0.RDS")

DefaultAssay(sc_combined) <- "integrated"
sc_combined <- ScaleData(sc_combined, verbose = FALSE)
sc_combined <- RunPCA(sc_combined, npcs = 100, verbose = FALSE)
ElbowPlot(sc_combined, ndims = 50, reduction = "pca")

sc_combined <- RunUMAP(sc_combined, reduction = "pca", dims = 1:20)
sc_combined <- RunTSNE(sc_combined, reduction = "pca", dims = 1:20)
sc_combined <- FindNeighbors(sc_combined, reduction = "pca",  dims = 1:20)

for(i in seq(0.1,2,by=0.1)){
  sc_combined <- FindClusters(sc_combined,verbose = FALSE, algorithm = 1,resolution=i)
}

pp_tree <- clustree(sc_combined@meta.data, prefix = "integrated_snn_res.")
ggsave(plot=pp_tree, filename = "/data//embryo_dev/myprocess/cluster/sc_combined_clustree_diff_resolution.pdf"), height=16, width=10)

Idents(sc_combined) <- "integrated_snn_res.1"
DefaultAssay(sc_combined) <- "RNA"
sc_res1_markers <- FindAllMarkers(sc_combined, only.pos=TRUE)

# filter undefined clusters with no specific markers
sc2 <- subset(sc_combined, integrated_snn_res.1 != "15")




