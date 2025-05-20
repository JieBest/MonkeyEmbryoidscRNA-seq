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





#################3. Integrate data
#################
#################
#################
#################
#################
