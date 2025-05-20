#This script integrate monkey embryo CS8-9 dataset (Zhai et al., 2022) and monkey embryoid
#Link to monkey embryo CS8-9 dataset is provided at National Center for Biotechnology Information (NCBI) Gene Expression Omnibus (GEO) under accession numbers GSE193007

library(Seurat)
library(harmony)
library(Seurat)
library(harmony)
library(ggplot2)
library(cowplot)
library(SCpubr)
library(dplyr)

sc_list <- readRDS("../../cr6ncbi_d17tod28_nFatureRNA1000_Singlet_nCountRNA30000.rds")
sc_list <- subset(sc_list, orig.ident %in% 
                    c("d17_1_cr6","d18_1_cr6","d21_1_cr6","d22_1_cr6",
                      "d23_1","d25_1_cr6","d25_2_cr6"))
new <- sc_list
new_counts <- new@assays$RNA@counts
new_meta <- new@meta.data

# load monkey in vivo data
zhai <- readRDS("../../public_data/mfas_embryo_CS8_CS11_gse193007_56636cells_update.RDS")
zhai <- subset(zhai, theiler_stage %in% c("CS8","CS9"))
zhai_counts <- zhai@assays$RNA@counts
zhai_meta <- zhai@meta.data

sharedGs <- intersect(row.names(new_counts), row.names(zhai_counts) ) # 17983 genes
new_counts <- new_counts[sharedGs,]
zhai_counts <- zhai_counts[sharedGs,]

new_meta_sub <- new_meta[,c(1,2,3,9)]
colnames(new_meta_sub) <- c("sample","nCount_RNA","nFeature_RNA","annotation")
new_meta_sub$study <- rep("blastoid",nrow(new_meta_sub))

zhai_meta_sub <- zhai_meta[,c(1,2,3,10)]
colnames(zhai_meta_sub) <- c("sample","nCount_RNA","nFeature_RNA","annotation")
zhai_meta_sub$study <- rep("invivo",nrow(zhai_meta_sub))


blastoid <- CreateSeuratObject(counts = new_counts,meta.data = new_meta_sub)
blastoid$orig.ident <- blastoid$sample
Idents(blastoid) <- 'orig.ident'

invivo <- CreateSeuratObject(counts = zhai_counts,meta.data = zhai_meta_sub)


sc_list =list(blastoid=blastoid, invivo=invivo)
save(sc_list, file = "raw_mfas_blastoid_d17-d25-2nod22-2_invivoCS89_sc_list.RData")


rm(list=ls())
load("raw_mfas_blastoid_d17-d25-2nod22-2_invivoCS89_sc_list.RData")


sc_list <- lapply(X = sc_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
  return(x)
})

features <- SelectIntegrationFeatures(object.list = sc_list, nfeatures = 3000)
sc_anchors <- FindIntegrationAnchors(object.list = sc_list, anchor.features = features)
save(sc_anchors, file = "mfas_blastoid_d17-d25-2nod22-2_invivoCS89_sc_anchors.RData")

# this command creates an 'integrated' data assay
sc_combined <- IntegrateData(anchorset = sc_anchors)
save(sc_combined,file="mfas_blastoid_d17-d25-2nod22-2_invivoCS89_sc_combined0.RData")


# library(Seurat)
# library(harmony)
# library(ggplot2)
# library(cowplot)
# library(SCpubr)
# library(dplyr)

setwd("/data/embryo_dev/myprocess/integrate_cr6ncbi_d17-1_d25-2_nod22-2/blastoid_invivo/")

load("mfas_blastoid_d17-d25-2nod22-2_invivoCS89_sc_combined0.RData")
alldata <- sc_combined
rm(sc_combined);gc()

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
ggsave(filename = "blastoid_invivoCS89_ElbowPlot.png", width = 7.68, height = 3.98, bg = "white")



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
ggsave(pp, filename = "blastoidD17-1_D25-2_invivoCS89_noharmony_res04.umap.png",
       width = 12.69, height=4.89)


save(alldata, file="blastoidD17-1_D25-2_invivoCS89_noharmony_withCluster.RData")

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

alldata_harmony <- FindClusters(alldata_harmony,resolution = 0.4) %>% identity()
DimPlot(object = alldata_harmony, reduction = "umap", 
        raster = F, label = TRUE, pt.size = 0.01)
ggsave(filename = "alldata_harmony_cluster_colByRes04cluster.png", width = 8, height = 6)


alldata_harmony <- FindClusters(alldata_harmony,resolution = 0.8) %>% identity()
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
table(alldata_harmony$annotation,alldata_harmony$integrated_snn_res.0.4)

save(alldata_harmony, file="mfas_blastoid_D17-1-D25-2_invivoCS89_alldata_harmony.RData")


alldata_harmony_meta <- alldata_harmony@meta.data
save(alldata_harmony_meta, file = "mfas_blastoid_D17-1-D25-2_invivoCS89_alldata_harmony_metainfo.RData")

library(SCpubr)
studyClusters <- unique(alldata_harmony$annotation)
for(i in 1:length(studyClusters)){
  pp <- do_DimPlot(alldata_harmony, raster = F, group.by = "annotation",idents.keep = studyClusters[i],pt.size=0.001)
  ggsave(pp ,filename=paste0("mfas_blastoid_D17-1-D25-2_invivoCS89_alldata_harmony_",studyClusters[i],"_umap.png"), width = 8, height = 6)
}

res04Clusters <- unique(alldata_harmony$integrated_snn_res.0.4)
for(i in 1:length(res04Clusters)){
  pp <- do_DimPlot(alldata_harmony, raster = F, group.by = "integrated_snn_res.0.4",idents.keep = res04Clusters[i],pt.size=0.001)
  ggsave(pp ,filename=paste0("mfas_blastoid_D17-1-D25-2_invivoCS89_alldata_harmony_res04Cluster",res04Clusters[i],"_umap.png"), width = 8, height = 6)
}

save.image("blastoid_D17-1-D25-2_invivoCS89_harmony_integrated.RData")

# Adjusting metadata for Harmony data
alldata_harmony@meta.data$study <- factor(alldata_harmony@meta.data$study, levels = c("blastoid", "invivo"))
toplist <- names(alldata_harmony$study[alldata_harmony$study == "blastoid"])
# Customized UMAP plots with specific color schemes
# Changing colors of the clusters and generating multiple UMAP plots
pp <- DimPlot(object = alldata_harmony, reduction = "umap", cols = c("#2979A9", "#C04844"), raster = FALSE, 
              order = toplist, group.by = "study", pt.size = 0.001) +labs(title = "")
#pp[[1]]$layers[[1]]$mapping$alpha= 0.7
ggsave(pp, filename = "alldata_harmony_cluster_changecolor_threecolor1.pdf", width = 8, height = 6)
ggsave(pp, filename = "alldata_harmony_cluster_changecolor_threecolor1.png", width = 8, height = 6)



toplist <- names(alldata_harmony$study[alldata_harmony$study == "invivo"])
pp <- DimPlot(object = alldata_harmony, reduction = "umap", cols = c("#C04844","#2979A9"), raster = FALSE, 
              order = toplist, group.by = "study", pt.size = 0.001) +labs(title = "")
#pp[[1]]$layers[[1]]$mapping$alpha= 0.7
ggsave(pp, filename = "alldata_harmony_cluster_changecolor_threecolor2.pdf", width = 8, height = 6)
ggsave(pp, filename = "alldata_harmony_cluster_changecolor_threecolor2.png", width = 8, height = 6)

pp <- do_DimPlot(alldata_harmony, group.by = "study", idents.keep = "invivo", 
           raster = FALSE,
           pt.size = 0.001, colors.use = c("invivo"="#C04844"),
           na.value = "lightgrey",
           legend.position = "none")
ggsave(pp, filename = "alldata_harmony_cluster_changecolor_threecolor3.pdf", width = 8, height = 6)

pp <- do_DimPlot(alldata_harmony, group.by = "study", idents.keep = "invivo", 
                 raster = FALSE,
                 pt.size = 0.001, colors.use = c("invivo"="#2979A9"),
                 na.value = "lightgrey",
                 legend.position = "none")
ggsave(pp, filename = "alldata_harmony_cluster_changecolor_threecolor4.pdf", width = 8, height = 6)


mergedAnno <- alldata_harmony$integrated_snn_res.0.4
mergedAnno <- recode(mergedAnno,
                     "0"="Mes","1"="Al","2"="ys.Meso",
                     "3"="Progenitors(ECT/EPI/PS)",
                     "4"="LP.Meso","5"="SE","6"="AM",
                     "7"="DE/VE",
                     "8"="Progenitors(ECT/EPI/PS)",
                     "9"="Gut","10"="Exe.Meso","11"="ys.Endo1",
                     "12"="Nas.Meso","13"="Al",
                     "14"="EC","15"="Blood","16"="ys.Endo2",
                     "17"="PGC","18"="Node","19"="Mes","20"="Mes")
mergedAnno <- factor(mergedAnno, levels = c("Progenitors(ECT/EPI/PS)", "PGC",
                                "DE/VE","Gut","ys.Endo1", "ys.Endo2",
                                "Nas.Meso","LP.Meso","Exe.Meso","Mes",
                                "Al","ys.Meso","Node",
                                "AM","SE",
                                "EC","Blood"))
# mergedAnno_cols <- c("#605546","#926358",
#                      "#214781","#DD90B4","#969696","#1B1B1B",
#                      "#CD6697","#CAE5F5","#BF7A33","#C1A897",
#                      "#806DA8","#EBB78D","#BC9C71",
#                      "#E9BEC9","#DF873E",
#                      "#374D2C","#D9987C")
mergedAnno_cols <- c("#8A554B","#C05028","#7B2063","#DCAB6A",
                     "#52B234","#7BD79E","#21B5C2","#FA8072",
                     "#754E9D","#E4CBF2","#0C6233","#80D7E1",
                     "#BC9C71","#5C9851","#F29694","#EF7C21","#CE75AB")
alldata_harmony$mergedAnno <- mergedAnno

umap_theme <- theme( 
  axis.line=element_blank(), 
  axis.text.x=element_blank(), 
  axis.text.y=element_blank(), 
  axis.ticks=element_blank(), 
  axis.title.x=element_blank(), 
  axis.title.y=element_blank(), 
  panel.background=element_blank(), 
  panel.border=element_blank(), 
  panel.grid.major=element_blank(), 
  panel.grid.minor=element_blank() 
) 
alpha.use <- 0.6
p <- DimPlot(alldata_harmony, raster = F, pt.size = 0.001, label = F, 
             cols = mergedAnno_cols,
             group.by = "mergedAnno")+labs(title = "") + umap_theme
p$layers[[1]]$mapping$alpha <- alpha.use
p <- p + scale_alpha_continuous(range = alpha.use, guide = F)

p <- do_DimPlot(alldata_harmony, raster = F, pt.size = 0.001, label = F, 
           colors.use =  c("Progenitors(ECT/EPI/PS)"="#605546","PGC"="#926358",
                    "DE/VE"="#214781","Gut"="#DD90B4",
                    "ys.Endo1"="#969696","ys.Endo2"="#1B1B1B",
                    "Nas.Meso"="#CD6697","LP.Meso"="#CAE5F5",
                    "Exe.Meso"="#BF7A33","Mes"="#C1A897",
                    "Al"="#806DA8","ys.Meso"="#EBB78D","Node"="#BC9C71",
                    "AM"="#E9BEC9","SE"="#DF873E",
                    "EC"="#374D2C","Blood"="#D9987C"),
           group.by = "mergedAnno")
p$layers[[1]]$mapping$alpha <- alpha.use
p <- p + scale_alpha_continuous(range = alpha.use, guide = F)

save_plot("/data/embryo_dev/myprocess/integrate_cr6ncbi_d17-1_d25-2_nod22-2/blastoid_invivo/blastoid_pc20_d17-1_d25-2_nod22-2_invivoCS89_celltype_UMAP_clustering.png", p, base_height = 8, base_aspect_ratio = 1.1, base_width = NULL, dpi=600)
save_plot("/data/embryo_dev/myprocess/integrate_cr6ncbi_d17-1_d25-2_nod22-2/blastoid_invivo/blastoid_pc20_d17-1_d25-2_nod22-2_invivoCS89_celltype_UMAP_clustering.pdf", p, base_height = 8, base_aspect_ratio = 1)

alpha.use=0.6
toplist <- names(alldata_harmony$study[alldata_harmony$study == "invivo"])
pp <- DimPlot(object = alldata_harmony, reduction = "umap", cols = c( "#E9BEC9","#CAE5F5"), raster = FALSE, 
              order = toplist, group.by = "study", pt.size = 0.001) +labs(title = "")+umap_theme
pp$layers[[1]]$mapping$alpha <- alpha.use
pp <- pp + scale_alpha_continuous(range = alpha.use, guide = F)

save_plot("/data/embryo_dev/myprocess/integrate_cr6ncbi_d17-1_d25-2_nod22-2/blastoid_invivo/alldata_harmony_cluster_changecolor_threecolor3.png", 
          pp, base_height = 8, base_aspect_ratio = 1.2, base_width = NULL, dpi=600)


mergedAnnoPerSample <- as.data.frame(table(alldata_harmony$mergedAnno, alldata_harmony$study))
colnames(mergedAnnoPerSample) <- c("cluster","study","number")

p1 <- ggplot(mergedAnnoPerSample,aes(x=study,y=number))+geom_bar(stat = "identity",aes(fill=cluster),position = "fill", width = 0.3) +
  scale_fill_manual(values = mergedAnno_cols) +
  xlab(label = "") +ylab(label = "Percentage")+labs(title = "")+theme_bw() +NoLegend()+
  theme(axis.text.x = element_text(angle=45, hjust = 1,vjust = 1, colour = "black"),
        axis.title.y = element_text(colour="black"),
        axis.line=element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank())
p1$layers[[1]]$mapping$alpha <- alpha.use
p1 <- p1 + scale_alpha_continuous(range = alpha.use, guide = F)
p1
save_plot("/data/embryo_dev/myprocess/integrate_cr6ncbi_d17-1_d25-2_nod22-2/blastoid_invivo/alldata_harmony_clusterPercentagePerStudy.png", 
          p1, base_height = 8, base_aspect_ratio = 0.9, base_width = NULL, dpi=600)


save(alldata_harmony, file = "/data/embryo_dev/myprocess/integrate_cr6ncbi_d17-1_d25-2_nod22-2/blastoid_invivo/blastoid_D17-1-D25-2_invivoCS89_harmony_withClusterInfo.RData")

save.image("blastoid_D17-1-D25-2_invivoCS89_harmony_integrated_240129.RData")


mergedAnnoPerSample$percentage <- mergedAnnoPerSample$number
for(study in unique(mergedAnnoPerSample$study)){
  mergedAnnoPerSample$percentage[mergedAnnoPerSample$study==study] <- 
    mergedAnnoPerSample$percentage[mergedAnnoPerSample$study==study]/sum(mergedAnnoPerSample$percentage[mergedAnnoPerSample$study==study])
}

cor(mergedAnnoPerSample$percentage[mergedAnnoPerSample$study=="blastoid"], 
    mergedAnnoPerSample$percentage[mergedAnnoPerSample$study=="invivo"], method = "spearman")

mergedAnnoPerSample$cluster <- factor(mergedAnnoPerSample$cluster,levels = rev(levels(mergedAnnoPerSample$cluster)))
p1 <- ggplot(mergedAnnoPerSample,aes(x=cluster,y=percentage, fill=study))+
  geom_bar(stat = "identity", position = "dodge" ,width = 0.7) +
  scale_fill_manual(values = c("#2979A9","#B74743"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=1,color="black"),
        axis.text.y = element_text(color = "black"),
        axis.line=element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank())+
  xlab(label = "")+ylab(label = "") + RotatedAxis() + coord_flip()
p1$layers[[1]]$mapping$alpha <- 0.8
p1 <- p1 + scale_alpha_continuous(range = 0.8, guide = F)
p1
save_plot("/data/embryo_dev/myprocess/integrate_cr6ncbi_d17-1_d25-2_nod22-2/blastoid_invivo/alldata_harmony_clusterPercentagePerStudy_dodge.pdf", 
          p1, base_height = 8, base_aspect_ratio = 0.8, base_width = NULL, dpi=600)

library(ggpubr)
mergedAnnoPerSample$percentage <- as.numeric(format(mergedAnnoPerSample$percentage,digits = 2))
mergedAnnoPerSample$labels <- paste0(mergedAnnoPerSample$percentage*100,"%")
p1 <- ggdonutchart(data = mergedAnnoPerSample[mergedAnnoPerSample$study=="invivo",],
                   x="number", fill="cluster",label = "labels",
                   lab.pos = "in", color = "black",lab.font = "black",
                   palette = mergedAnno_cols) +theme(legend.position = "right") +NoLegend()
p2 <- ggdonutchart(data = mergedAnnoPerSample[mergedAnnoPerSample$study=="blastoid",],
                   x="number", fill="cluster",label = "labels",
                   lab.pos = "in", color = "black",lab.font = "black",
                   palette = mergedAnno_cols) +theme(legend.position = "right") 
pp <- p1|p2
save_plot("/data/embryo_dev/myprocess/integrate_cr6ncbi_d17-1_d25-2_nod22-2/blastoid_invivo/alldata_harmony_clusterPercentagePerStudy_donutchart_240802.pdf", 
          pp, base_height = 6.98, base_aspect_ratio =1.6, base_width = NULL, dpi=600)

write.table(mergedAnnoPerSample,
            file = "/data/embryo_dev/myprocess/integrate_cr6ncbi_d17-1_d25-2_nod22-2/blastoid_invivo/alldata_harmony_clusterPercentagePerStudy.tsv",
            sep="\t",quote=FALSE, row.names = FALSE)

# remove Mes and Al
mergedCluster2 <- setdiff( unique(mergedAnnoPerSample$cluster), c("Mes","Al"))
mergedAnnoPerSample2 <- mergedAnnoPerSample[mergedAnnoPerSample$cluster %in% mergedCluster2,]
for(study in unique(mergedAnnoPerSample2$study)){
  mergedAnnoPerSample2$percentage[mergedAnnoPerSample2$study==study] <- 
    as.numeric(format(mergedAnnoPerSample2$number[mergedAnnoPerSample2$study==study]/sum(mergedAnnoPerSample2$number[mergedAnnoPerSample2$study==study]),
           digits = 2))
}
cor(mergedAnnoPerSample2$percentage[mergedAnnoPerSample2$study=="blastoid"], 
    mergedAnnoPerSample2$percentage[mergedAnnoPerSample2$study=="invivo"], method = "pearson")

#---- bar plot
mergedCluster2_cols <- c("#8A554B","#C05028","#7B2063","#DCAB6A",
                     "#52B234","#7BD79E","#21B5C2","#FA8072",
                     "#754E9D","#80D7E1",
                     "#BC9C71","#5C9851","#F29694","#EF7C21","#CE75AB")
mergedAnnoPerSample2$study <- factor(mergedAnnoPerSample2$study,levels = c("invivo","blastoid"))
p1 <- ggplot(mergedAnnoPerSample2,aes(x=study,y=number))+geom_bar(stat = "identity",aes(fill=cluster),position = "fill", width = 0.3) +
  scale_fill_manual(values = mergedCluster2_cols) +
  xlab(label = "") +ylab(label = "Cell proportion")+labs(title = "")+theme_bw() +NoLegend()+
  theme(axis.text.x = element_text(angle=45, hjust = 1,vjust = 1, colour = "black",size = 10),
        axis.title.y = element_text(colour="black", size = 10),
        axis.line=element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank()) 
save_plot("/data/embryo_dev/myprocess/integrate_cr6ncbi_d17-1_d25-2_nod22-2/blastoid_invivo/alldata_harmony_clusterPercentagePerStudy_rmMes_Al.pdf", 
          p1, base_height = 8, base_aspect_ratio = 0.8, base_width = NULL, dpi=600)


mergedAnnoPerSample2$cluster <- factor(mergedAnnoPerSample2$cluster,levels = rev(levels(mergedAnnoPerSample2$cluster)))
p1 <- ggplot(mergedAnnoPerSample2,aes(x=cluster,y=percentage, fill=study))+
  geom_bar(stat = "identity", position = "dodge" ,width = 0.7) +
  scale_fill_manual(values = c("#2979A9","#B74743"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=1,color="black"),
        axis.text.y = element_text(color = "black"),
        axis.line=element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank())+
  xlab(label = "")+ylab(label = "") + RotatedAxis() + coord_flip()
p1$layers[[1]]$mapping$alpha <- 0.8
p1 <- p1 + scale_alpha_continuous(range = 0.8, guide = F)
p1
save_plot("/data/embryo_dev/myprocess/integrate_cr6ncbi_d17-1_d25-2_nod22-2/blastoid_invivo/alldata_harmony_clusterPercentagePerStudy_dodge_rmMesAl.pdf", 
          p1, base_height = 8, base_aspect_ratio = 0.8, base_width = NULL, dpi=600)

library(ggpubr)
mergedAnnoPerSample2$labels <- paste0(mergedAnnoPerSample2$percentage*100,"%")
p1 <- ggdonutchart(data = mergedAnnoPerSample2[mergedAnnoPerSample2$study=="invivo",],
             x="number", fill="cluster",label = "labels",
             lab.pos = "in", color = "black",lab.font = "black",
             palette = mergedCluster2_cols) +theme(legend.position = "right") +NoLegend()
p2 <- ggdonutchart(data = mergedAnnoPerSample2[mergedAnnoPerSample2$study=="blastoid",],
                   x="number", fill="cluster",label = "labels",
                   lab.pos = "in", color = "black",lab.font = "black",
                   palette = mergedCluster2_cols) +theme(legend.position = "right") 
pp <- p1|p2
save_plot("/data/embryo_dev/myprocess/integrate_cr6ncbi_d17-1_d25-2_nod22-2/blastoid_invivo/alldata_harmony_clusterPercentagePerStudy_rmMesAl_donutchart_240802.pdf", 
          pp, base_height = 6.98, base_aspect_ratio =1.6, base_width = NULL, dpi=600)

write.table(mergedAnnoPerSample2,
            file = "/data/embryo_dev/myprocess/integrate_cr6ncbi_d17-1_d25-2_nod22-2/blastoid_invivo/alldata_harmony_clusterPercentagePerStudy_rmMesAl.tsv",
            sep="\t",quote=FALSE, row.names = FALSE)

#-------------------------------------------------------------------------------
# Plot for figure3
p1 <- do_DimPlot(alldata_harmony, raster = F, pt.size = 0.001, label = F, 
           colors.use =  c("Progenitors(ECT/EPI/PS)"="#8A554B","PGC"="#C05028",
                           "DE/VE"="#7B2063","Gut"="#DCAB6A",
                           "ys.Endo1"="#52B234","ys.Endo2"="#7BD79E",
                           "Nas.Meso"="#21B5C2","LP.Meso"="#FA8072",
                           "Exe.Meso"="#754E9D","Mes"="#E4CBF2",
                           "Al"="#0C6233","ys.Meso"="#80D7E1","Node"="#BC9C71",
                           "AM"="#5C9851","SE"="#F29694",
                           "EC"="#EF7C21","Blood"="#CE75AB"),
           group.by = "mergedAnno")

mergedAnno_cols <- c("#8A554B","#C05028","#7B2063","#DCAB6A",
                     "#52B234","#7BD79E","#21B5C2","#FA8072",
                     "#754E9D","#E4CBF2","#0C6233","#80D7E1",
                     "#BC9C71","#5C9851","#F29694","#EF7C21","#CE75AB")
p1 <- DimPlot(alldata_harmony, raster = F, pt.size = 0.001, label = F, 
        cols = mergedAnno_cols,
        group.by = "mergedAnno")+labs(title = "")+xlab(label = "UMAP1") +ylab(label = "UMAP2")

save_plot("alldata_harmony_cluster_fig3D_240308.pdf", 
          p1, base_height = 8, base_aspect_ratio = 1.5)


p2 <- p1+NoLegend()
save_plot("alldata_harmony_cluster_fig3D_240308_noLegend.pdf", 
          p2, base_height = 8, base_aspect_ratio = 1.1)

p1 <- DimPlot(alldata_harmony, raster = F, pt.size = 0.0001, label = F, 
              cols = mergedAnno_cols,
              group.by = "mergedAnno",split.by = "study", ncol = 2)+labs(title = "")+xlab(label = "UMAP1") +ylab(label = "UMAP2") +NoLegend()
p1$layers[[1]]$mapping$alpha <- 0.6
p1 <- p1 + scale_alpha_continuous(range = 0.6, guide = F)
save_plot("alldata_harmony_cluster_fig3C_240509_noLegend.pdf", 
          p1, base_height = 4, base_aspect_ratio = 1.8)

mergedAnnoPerSample <- as.data.frame(table(alldata_harmony$mergedAnno, alldata_harmony$study))
colnames(mergedAnnoPerSample) <- c("cluster","study","number")

p1 <- ggplot(mergedAnnoPerSample,aes(x=study,y=number))+geom_bar(stat = "identity",aes(fill=cluster),position = "fill", width = 0.9) +
  scale_fill_manual(values = mergedAnno_cols) +
  xlab(label = "") +ylab(label = "Percentage")+labs(title = "")+theme_bw() +NoLegend()+
  theme(axis.text.x = element_text(angle=45, hjust = 1,vjust = 1, colour = "black"),
        axis.title.y = element_text(colour="black"),
        axis.line=element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank()) +coord_flip()
# p1$layers[[1]]$mapping$alpha <- alpha.use
# p1 <- p1 + scale_alpha_continuous(range = alpha.use, guide = F)
p1
save_plot("alldata_harmony_clusterPercentagePerStudy.pdf", 
          p1, base_height =3, base_aspect_ratio =2, base_width = NULL, dpi=600)


toplist <- names(alldata_harmony$study[alldata_harmony$study == "invivo"])
pp <- DimPlot(object = alldata_harmony, reduction = "umap", cols = c( "#2979A9",  "#C04844"), raster = FALSE, 
              order = toplist, group.by = "study", pt.size = 0.001) +
  labs(title = "")+xlab(label = "UMAP1") +ylab(label = "UMAP2")
pp$layers[[1]]$mapping$alpha <- alpha.use
pp <- pp + scale_alpha_continuous(range = alpha.use, guide = F)

save_plot("alldata_harmony_cluster_changecolor_threecolor3.pdf", 
          pp, base_height = 8, base_aspect_ratio = 1.2)

pp <- pp+NoLegend()
save_plot("alldata_harmony_cluster_changecolor_threecolor3_240308_noLegend.pdf", 
          pp, base_height = 8, base_aspect_ratio = 1.1)

DefaultAssay(alldata_harmony) <- "RNA"

featuresPlot <- c("POU5F1","SOX2",
                  "TFAP2C","NANOG",
                  "SOX17","FOXA2",
                  "ISL1", # Gut "CDX2","CDH2","OSR1",
                  "TTR","AFP", # ys.Endo
                  "T","MESP1","TBX6",
                  "GATA6","HAND1", 
                  "ANXA1",
                  "CHRD", "NOTO",
                  "TFAP2A", 
                  "KDR","PECAM1",
                  "HBE1","HBZ")
VlnPlot(alldata_harmony, features = featuresPlot, split.by = "study", group.by = "mergedAnno",
        pt.size = 0, stack = TRUE, combine = FALSE,flip = T) + 
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=1,color="black"),
        axis.text.y = element_text(color = "black"))+
  xlab(label = "")+ylab(label = "") +scale_fill_manual(values = c("#2979A9","#B74743")) #
ggsave(filename = "alldata_harmony_conservedFeatures_vlnplot_240802.pdf",
       height = 9.89,width = 7.89)


# adding OTX2, OSR1, NODAL, FGF17, WNT3 and WNT3A
featuresPlot <- c("POU5F1","SOX2",
                  "TFAP2C","NANOG",
                  "SOX17","FOXA2",
                  "ISL1", # Gut "CDX2","CDH2","OSR1",
                  "TTR","AFP", # ys.Endo
                  "T","MESP1","TBX6",
                  "GATA6","HAND1", 
                  "ANXA1",
                  "CHRD", "NOTO",
                  "OTX2","OSR1","NODAL","FGF17","WNT3","WNT3A",
                  "TFAP2A", 
                  "KDR","PECAM1",
                  "HBE1","HBZ")
VlnPlot(alldata_harmony, features = featuresPlot, split.by = "study", group.by = "mergedAnno",
        pt.size = 0, stack = TRUE, combine = FALSE,flip = T) + 
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=1,color="black"),
        axis.text.y = element_text(color = "black"))+
  xlab(label = "")+ylab(label = "") +scale_fill_manual(values = c("#2979A9","#B74743")) #
ggsave(filename = "alldata_harmony_conservedFeatures_vlnplot_250512.pdf",
       height = 9.89,width = 7.89)
ggsave(filename = "alldata_harmony_conservedFeatures_vlnplot_250512.tiff",
       bg = "white",
       height = 9.89,width = 7.89)

# DotPlot(alldata_harmony,
#         features = c("POU5F1","SOX2","TFAP2A","TFAP2C"),
#         group.by = "mergedAnno", 
#         cols = c("blue", "red"), 
#         dot.scale = 8,
#         split.by = "study") +  RotatedAxis() + xlab(label = "")+ylab(label = "")


load("blastoid_D17-1-D25-2_invivoCS89_harmony_integrated_240129.RData")

dir.create("AMAseletcedFeatures_umap")
ama_features2.1 <- c("NODAL", "WNT3", "WNT3A", "WNT5A","WNT5B", "LEFTY1", "LOC102144954",
                     "EOMES", "T", "SNAI1", "SNAI2",
                     "BMP4",
                     "MESP1", "TBX5", "ISL1", "NKX2-5", "GATA4", "GATA6",
                     "KDR", "CDH5",
                     "TBX1",
                     "FOXF1", "HAND1",
                     "OSR1", "LHX1",
                     "MSGN1", "TBX6",
                     "SHH", "NOTO", "FOXA2","CDH1", "EPCAM",
                     "TFAP2C", "NANOG","SOX17",
                     "SOX9", "TWIST1",
                     "SOX1", "SOX2", "FOXG1", "EN1", "OLIG2", "NKX2-1", "PAX6",
                     "TCF15", "MEOX1","MEOX2",
                     "PAX1", "PAX9", "TBX4", "EVX1", "EGR2", "FGF5", 
                     "MESP2",  "SIX1", "SIX3", "TWIST2", "DLX5")

DefaultAssay(alldata_harmony) <- "RNA"
# Gastrulation (signals in PS): NODAL, WNT3, WNT5, LEFTY1, 2.
p <- FeaturePlot(alldata_harmony,pt.size = 0.0001, 
            raster = FALSE, 
            features = c("NODAL","WNT3", "WNT3A", "WNT5A","WNT5B", "LEFTY1"), 
            split.by = "study",
            cols = c("lightgrey","red"))
save_plot(filename = "./AMAseletcedFeatures_umap/Gastrulation_genes.png",
          plot = p,
          base_height = 18,
          base_width = 7)

# Nas Mes: EOMES, TBXT, SNAI1, SNAI2
p <- FeaturePlot(alldata_harmony,pt.size = 0.0001, 
                 raster = FALSE, 
                 features = c("EOMES", "T", "SNAI1", "SNAI2"), 
                 split.by = "study",
                 cols = c("lightgrey","red"))
save_plot(filename = "./AMAseletcedFeatures_umap/NasMeso_genes.png",
          plot = p,
          base_height = 12,
          base_width = 7)

# ExtraEmbr Mes: BMP4
p <- FeaturePlot(alldata_harmony,pt.size = 0.0001, 
                 raster = FALSE, 
                 features = c("BMP4"), 
                 split.by = "study",
                 cols = c("lightgrey","red"))
save_plot(filename = "./AMAseletcedFeatures_umap/ExtraEmbrMes_genes.png",
          plot = p,
          base_height = 3,
          base_width = 7)

#Hemogenic: KDR, CDH5
p <- FeaturePlot(alldata_harmony,pt.size = 0.0001, 
                 raster = FALSE, 
                 features = c("KDR", "CDH5"), 
                 split.by = "study",
                 cols = c("lightgrey","red"))
save_plot(filename = "./AMAseletcedFeatures_umap/Hemogenic_genes.png",
          plot = p,
          base_height = 6,
          base_width = 7)

#Cardiac: MESP1, TBX5, ISL1, NKX2.5, GATA4, 6
p <- FeaturePlot(alldata_harmony,pt.size = 0.0001, 
                 raster = FALSE, 
                 features = c("MESP1", "TBX5", "ISL1", "NKX2-5", "GATA4", "GATA6"), 
                 split.by = "study",
                 cols = c("lightgrey","red"))
save_plot(filename = "./AMAseletcedFeatures_umap/Cardiac_genes.png",
          plot = p,
          base_height = 18,
          base_width = 7)

#Lat Plate: FOXF1, HAND1
p <- FeaturePlot(alldata_harmony,pt.size = 0.0001, 
                 raster = FALSE, 
                 features = c("FOXF1", "HAND1"), 
                 split.by = "study",
                 cols = c("lightgrey","red"))
save_plot(filename = "./AMAseletcedFeatures_umap/LatPlateMeso_genes.png",
          plot = p,
          base_height = 6,
          base_width = 7)

#Inter Mes: OSR1, LHX1
p <- FeaturePlot(alldata_harmony,pt.size = 0.0001, 
                 raster = FALSE, 
                 features = c("OSR1", "LHX1"), 
                 split.by = "study",
                 cols = c("lightgrey","red"))
save_plot(filename = "./AMAseletcedFeatures_umap/InterMeso_genes.png",
          plot = p,
          base_height = 6,
          base_width = 7)

#Paraxial Mes: MSGN, TBX6
p <- FeaturePlot(alldata_harmony,pt.size = 0.0001, 
                 raster = FALSE, 
                 features = c("MSGN1", "TBX6"), 
                 split.by = "study",
                 cols = c("lightgrey","red"))
save_plot(filename = "./AMAseletcedFeatures_umap/ParaxialMes_genes.png",
          plot = p,
          base_height = 6,
          base_width = 7)

#Node: SHH, FOXA2. NOTO
p <- FeaturePlot(alldata_harmony,pt.size = 0.0001, 
                 raster = FALSE, 
                 features = c("SHH", "FOXA2","NOTO"), 
                 split.by = "study",
                 cols = c("lightgrey","red"))
save_plot(filename = "./AMAseletcedFeatures_umap/Node_genes.png",
          plot = p,
          base_height = 9,
          base_width = 7)

#Endoderm: FOXA2, CDH1, EPCAM
p <- FeaturePlot(alldata_harmony,pt.size = 0.0001, 
                 raster = FALSE, 
                 features = c("CDH1", "FOXA2","EPCAM"), 
                 split.by = "study",
                 cols = c("lightgrey","red"))
save_plot(filename = "./AMAseletcedFeatures_umap/Endoderm_genes.png",
          plot = p,
          base_height = 9,
          base_width = 7)

#PGCs: TFAP2C, NANOG, SOX17
p <- FeaturePlot(alldata_harmony,pt.size = 0.0001, 
                 raster = FALSE, 
                 features = c("TFAP2C", "NANOG","SOX17"), 
                 split.by = "study",
                 cols = c("lightgrey","red"))
save_plot(filename = "./AMAseletcedFeatures_umap/PGCs_genes.png",
          plot = p,
          base_height = 9,
          base_width = 7)

#Neural Crest: SOX9, TWI1
p <- FeaturePlot(alldata_harmony,pt.size = 0.0001, 
                 raster = FALSE, 
                 features = c("SOX9", "TWIST1"), 
                 split.by = "study",
                 cols = c("lightgrey","red"))
save_plot(filename = "./AMAseletcedFeatures_umap/NeuralCrest_genes.png",
          plot = p,
          base_height = 6,
          base_width = 7)

#Neural: SOX1, SOX2, FOXG1, EN1, OLIG2, NKX2.1, PAX6
p <- FeaturePlot(alldata_harmony,pt.size = 0.0001, 
                 raster = FALSE, 
                 features = c("SOX1", "SOX2", "FOXG1", "EN1", "OLIG2", "NKX2-1", "PAX6"), 
                 split.by = "study",
                 cols = c("lightgrey","red"))
save_plot(filename = "./AMAseletcedFeatures_umap/Neural_genes.png",
          plot = p,
          base_height = 21,
          base_width = 7)

#Somites: TCF15, MEOX
p <- FeaturePlot(alldata_harmony,pt.size = 0.0001, 
                 raster = FALSE, 
                 features = c("TCF15", "MEOX1","MEOX2"), 
                 split.by = "study",
                 cols = c("lightgrey","red"))
save_plot(filename = "./AMAseletcedFeatures_umap/Somites_genes.png",
          plot = p,
          base_height = 9,
          base_width = 7)

#New ones: PAX1, PAX9, TBX4, EVX1, EGR2, FGF5, MESP1, MESP2,  SIX1, SIX3, TWIST2, DLX5
p <- FeaturePlot(alldata_harmony,pt.size = 0.0001, 
                 raster = FALSE, 
                 features = c("PAX1", "PAX9", "TBX4", "EVX1", "EGR2", "FGF5"), 
                 split.by = "study",
                 cols = c("lightgrey","red"))
save_plot(filename = "./AMAseletcedFeatures_umap/New1_genes.png",
          plot = p,
          base_height = 18,
          base_width = 7)

p <- FeaturePlot(alldata_harmony,pt.size = 0.0001, 
                 raster = FALSE, 
                 features = c("MESP2",  "SIX1", "SIX3", "TWIST2", "DLX5"), 
                 split.by = "study",
                 cols = c("lightgrey","red"))
save_plot(filename = "./AMAseletcedFeatures_umap/New2_genes.png",
          plot = p,
          base_height = 15,
          base_width = 7)

p <- FeaturePlot(alldata_harmony,pt.size = 0.0001, 
                 raster = FALSE, 
                 features = c("OSR1",  "LHX1", "TCF15"), 
                 split.by = "study",
                 cols = c("lightgrey","red"))
save_plot(filename = "./AMAseletcedFeatures_umap/New3_genes.pdf",
          plot = p,
          base_height = 9,
          base_width = 7)



load("/data/embryo_dev/myprocess/integrate_cr6ncbi_d17-1_d25-2_nod22-2/blastoid_invivo/blastoidD17-1_D25-2_invivoCS89_noharmony_withCluster.RData")
load("/data/embryo_dev/myprocess/integrate_cr6ncbi_d17-1_d25-2_nod22-2/blastoid_invivo/blastoid_D17-1-D25-2_invivoCS89_harmony_withClusterInfo.RData")

table(row.names(alldata@meta.data)==row.names(alldata_harmony@meta.data))
table(alldata$study)

alldata$alldataHarmonyAnno <-alldata_harmony$mergedAnno 
DimPlot(alldata, raster = F, pt.size = 0.001, label = T, 
        cols = mergedAnno_cols,
        group.by = "alldataHarmonyAnno")+labs(title = "")+xlab(label = "UMAP1") +ylab(label = "UMAP2")

toplist <- names(alldata$study[alldata$study == "invivo"])
p1 <- do_DimPlot(alldata, reduction = "pca", group.by = "study", idents.keep = "invivo", 
           raster = FALSE, pt.size = 0.001, colors.use = c("invivo"="#C04844"),
           na.value = "lightgrey", legend.position = "none")

p2 <- do_DimPlot(alldata, reduction = "umap", group.by = "study", idents.keep = "invivo", 
                 raster = FALSE, pt.size = 0.001, colors.use = c("invivo"="#C04844"),
                 na.value = "lightgrey", legend.position = "none")

alldata$alldataHarmonyAnno <-alldata_harmony$mergedAnno 
p3 <- do_DimPlot(alldata, raster = F, pt.size = 0.001, label = F, 
                 colors.use =  c("Progenitors(ECT/EPI/PS)"="#8A554B","PGC"="#C05028",
                                 "DE/VE"="#7B2063","Gut"="#DCAB6A",
                                 "ys.Endo1"="#52B234","ys.Endo2"="#7BD79E",
                                 "Nas.Meso"="#21B5C2","LP.Meso"="#FA8072",
                                 "Exe.Meso"="#754E9D","Mes"="#E4CBF2",
                                 "Al"="#0C6233","ys.Meso"="#80D7E1","Node"="#BC9C71",
                                 "AM"="#5C9851","SE"="#F29694",
                                 "EC"="#EF7C21","Blood"="#CE75AB"),
                 group.by = "alldataHarmonyAnno")+labs(title = "")

p4 <- do_DimPlot(alldata, raster = F, pt.size = 0.001, label = F, 
                 colors.use =  c("Progenitors(ECT/EPI/PS)"="#8A554B","PGC"="#C05028",
                                 "DE/VE"="#7B2063","Gut"="#DCAB6A",
                                 "ys.Endo1"="#52B234","ys.Endo2"="#7BD79E",
                                 "Nas.Meso"="#21B5C2","LP.Meso"="#FA8072",
                                 "Exe.Meso"="#754E9D","Mes"="#E4CBF2",
                                 "Al"="#0C6233","ys.Meso"="#80D7E1","Node"="#BC9C71",
                                 "AM"="#5C9851","SE"="#F29694",
                                 "EC"="#EF7C21","Blood"="#CE75AB"),
                 group.by = "alldataHarmonyAnno",legend.position = "none")

ggsave(p3, file="/data/embryo_dev/Summary/Nature_revision/r3/blastoid_D17-1-D25-2_invivoCS89_alldataNoHarmony_umap_labelByHarmonyAnno.pdf",
       height = 6.18, width = 6.18)

pp <- plot_grid(p1,p2,p4, ncol = 3)
ggsave(pp, file="/data/embryo_dev/Summary/Nature_revision/r3/blastoid_D17-1-D25-2_invivoCS89_alldataNoHarmony_pca_umap_panel3NoLegned.pdf",
       height = 4.38, width = 15)

rm(p1,p2,p3,p4,pp); gc()



alldataHarmonyMergedAnnoMeanexp <- AverageExpression(alldata_harmony,group.by = )



load("/data/embryo_dev/myprocess/integrate_cr6ncbi_d17-1_d25-2_nod22-2/blastoid_invivo/blastoid_D17-1-D25-2_invivoCS89_harmony_withClusterInfo.RData")

# 
library(speckle)
library(SingleCellExperiment)
library(scater)

allcounts <- alldata_harmony@assays$RNA@counts
sce_all <- SingleCellExperiment(assays = list(counts = allcounts))
sce_all$sample <- alldata_harmony$orig.ident
sce_all$group <- alldata_harmony$study
sce_all$cluster <- alldata_harmony$mergedAnno

rm(allcounts, alldata, alldata_harmony)
gc()
# Visualise the data

sce_all <- scater::logNormCounts(sce_all)
sce_all <- scater::runPCA(sce_all)
sce_all <- scater::runUMAP(sce_all)

pca1 <- scater::plotReducedDim(sce_all, dimred = "PCA", colour_by = "cluster") +
  ggtitle("Cell Type")
pca2 <- scater::plotReducedDim(sce_all, dimred = "PCA", colour_by = "group") +
  ggtitle("Study")
pca1 + pca2


umap1 <- scater::plotReducedDim(sce_all, dimred = "UMAP", 
                                colour_by = "cluster") + 
  ggtitle("Cell Type")
umap2 <- scater::plotReducedDim(sce_all, dimred = "UMAP", colour_by = "group") +
  ggtitle("Study")
umap1 + umap2


# Test for differences 

# Perform logit transformation
propeller(sce_all)
# Perform arcsin square root transformation
propeller(sce_all, transform="asin")
# 
propeller(clusters=sce_all$cluster, sample=sce_all$sample, group=sce_all$group)


# Visualise the results
plotCellTypeProps(sce_all)

props <- getTransformedProps(sce_all$cluster, sce_all$sample, transform="logit")
barplot(props$Proportions, 
        col = mergedAnno_cols,
        legend=FALSE,
        ylab="Proportions")

par(mfrow=c(3,6))
for(i in seq(1,17,1)){
  stripchart(props$Proportions[i,]~rep(c("vitro","vivo"),c(7,4)),
             vertical=TRUE, pch=16, method="jitter",
             col = c("orange","purple"),cex=2, ylab="Proportions")
  title(rownames(props$Proportions)[i])
}

# scDC
library(scDC)
res_scDC_noClust <- scDC_noClustering(sce_all$cluster, sce_all$sample, calCI = TRUE,
                                      calCI_method = c("percentile", "BCa", "multinom"),
                                      nboot = 10000,
                                      ncores = 8)
barplotCI(res_scDC_noClust, condition = rep(c("vitro","vivo"),c(7,4)))
densityCI(res_scDC_noClust, condition = rep(c("vitro","vivo"),c(7,4)))

res_GLM <- fitGLM(res_scDC_noClust,
                  condition = rep(c("vitro","vivo"),c(7,4)),
                  pairwise = FALSE)
summary(res_GLM$pool_res_fixed)
summary(res_GLM$pool_res_random)


# 

# Step 3
cmdstanr::check_cmdstan_toolchain(fix = TRUE) # Just checking system setting
cmdstanr::install_cmdstan()



library(dplyr)
library(sccomp)
library(ggplot2)
library(forcats)
library(tidyr)

sccomp_result = 
  sce_all |>
  sccomp_estimate(
    .data=sce_all,
    formula_composition = ~cluster, 
    .sample = ~sample, 
    .cell_group = ~group, 
    cores = 1 
  )  |> 
  sccomp_remove_outliers(cores = 1)  |> 
  sccomp_test()



# 

library("sceasy")
library("reticulate")
use_condaenv('scrna')
sceasy::convertFormat(sce_all, from="sce", to="anndata",
                      outFile='/data/embryo_dev/myprocess/integrate_cr6ncbi_d17-1_d25-2_nod22-2/blastoid_invivo/blastoidD17-1_D25-2_invivoCS89_sce_withHormonyCluster.h5ad')


rm(list=ls())

library(Seurat)
library(SeuratDisk)

load("/data/embryo_dev/myprocess/integrate_cr6ncbi_d17-1_d25-2_nod22-2/blastoid_invivo/blastoid_D17-1-D25-2_invivoCS89_harmony_withClusterInfo.RData")
DefaultAssay(alldata_harmony) <- "RNA"

alldata_harmony_meta <- alldata_harmony@meta.data
alldata_harmony_meta <- alldata_harmony_meta[,c('orig.ident','mergedAnno','study')]
colnames(alldata_harmony_meta) <- c("sample","cluster",'group')
alldata_harmony_meta$sample <- stringr::str_split(alldata_harmony_meta$sample,pattern = "_cr6",simplify = T)[,1]

write.csv(alldata_harmony_meta,
          file="/data/embryo_dev/myprocess/integrate_cr6ncbi_d17-1_d25-2_nod22-2/blastoid_invivo/blastoid_D17-1-D25-2_invivoCS89_harmony_metadf_col3.csv",
          quote=FALSE)

SaveH5Seurat(alldata_harmony,filename="/data/embryo_dev/myprocess/integrate_cr6ncbi_d17-1_d25-2_nod22-2/blastoid_invivo/blastoidD17-1_D25-2_invivoCS89_sce_withHormonyCluster.h5seurat", overwrite = TRUE)

#数据转为最终h5ad格式
Convert("/data/embryo_dev/myprocess/integrate_cr6ncbi_d17-1_d25-2_nod22-2/blastoid_invivo/blastoidD17-1_D25-2_invivoCS89_sce_withHormonyCluster.h5seurat", dest = "/data/embryo_dev/myprocess/integrate_cr6ncbi_d17-1_d25-2_nod22-2/blastoid_invivo/blastoidD17-1_D25-2_invivoCS89_sce_withHormonyCluster.h5ad", overwrite = TRUE)





load("/data/embryo_dev/myprocess/integrate_cr6ncbi_d17-1_d25-2_nod22-2/blastoid_invivo/blastoidD17-1_D25-2_invivoCS89_noharmony_withCluster.RData")
DefaultAssay(alldata)
# UMAP on integrated Assay
p1 <- DimPlot(alldata, reduction = "pca",
              raster = FALSE, 
              cols = c("lightgrey","#B74643"),
              pt.size = 0.1,group.by = "study") + labs(title = "") +NoLegend()
p2 <- DimPlot(alldata, reduction = "umap",
              raster = FALSE,
              cols = c("lightgrey","#B74643"),
              pt.size = 0.1,group.by = "study") + labs(title = "") +NoLegend()
p3 <- DimPlot(alldata, reduction = "umap",
              raster = FALSE, label = TRUE,
              pt.size = 0.1,group.by = "integrated_snn_res.0.4") +NoLegend() + labs(title = "")
pp <- cowplot::plot_grid(p1,p2,p3, ncol=3, rel_widths = c(1,1,1))
ggsave(pp, filename="alldata_Integratedassay_pca_umap.pdf",
       height=3.58, width=12.68)

# PCA plot of RNA assay
DefaultAssay(alldata) <- "RNA"
alldata <- NormalizeData(alldata)
alldata <- FindVariableFeatures(alldata, selection.method = "vst", nfeatures = 3000)
alldata <- ScaleData(alldata)
alldata <- RunPCA(alldata, npcs = 100)
alldata <- RunUMAP(alldata,dims = 1:100, reduction = "pca")
p4 <- DimPlot(alldata, reduction = "pca",
              raster = FALSE, 
              cols = c("lightgrey","#B74643"),
              pt.size = 0.1,group.by = "study") + labs(title = "")
p9 <- DimPlot(alldata, reduction = "umap",
              raster = FALSE, 
              cols = c("lightgrey","#B74643"),
              pt.size = 0.1,group.by = "study") + labs(title = "")
pp <- p4|p9
ggsave(pp, filename="alldata_RNAassay_pca_umap.pdf",
       height=3.48, width=9.98)


load("/data/embryo_dev/myprocess/integrate_cr6ncbi_d17-1_d25-2_nod22-2/blastoid_invivo/blastoid_D17-1-D25-2_invivoCS89_harmony_withClusterInfo.RData")
DefaultAssay(alldata_harmony)
#
p5 <- DimPlot(alldata_harmony, reduction = "harmony",
              raster = FALSE, cols = c("lightgrey","#B74643"),
              pt.size = 0.001,group.by = "study") + labs(title = "") +NoLegend()
ggsave(p5, filename ="alldataharmony_harmony.pdf",height = 3.58, width = 3.98)
p6 <- DimPlot(alldata_harmony, reduction = "umap",
              raster = FALSE, cols = c("lightgrey","#B74643"),
              pt.size = 0.1,group.by = "study") + labs(title = "")
p7 <- DimPlot(alldata_harmony, reduction = "umap",
              raster = FALSE, label = TRUE,
              pt.size = 0.1,group.by = "integrated_snn_res.0.4") +NoLegend()
ggsave(p7, filename="alldataharmony_res04_umap.pdf",
       height=3.58, width=3.98)

load("/data/embryo_dev/myprocess/integrate_cr6ncbi_d17-1_d25-2_nod22-2/blastoid_invivo/blastoid_D17-1-D25-2_invivoCS89_harmony_withClusterInfo.RData")
DefaultAssay(alldata_harmony) <- "RNA"

library(dplyr)
study2 <- alldata_harmony$study
study2 <- recode(study2,
                 "blastoid"="Embryoid",
                 "invivo"="Embryo(CS8-9)")
alldata_harmony$study2 <- study2

FeaturePlot(alldata_harmony,
            features = c("TBX5","ISL1","GATA4","SNAI1","TCF21"),
            pt.size = 0.001, raster = FALSE,
            cols = c("lightgrey","red"),
            split.by = "study2")

FeaturePlot(alldata_harmony,
            features = c("CDX2","BMP4","IRX3"),
            pt.size = 0.001, raster = FALSE,
            cols = c("lightgrey","red"),
            split.by = "study2")


FeaturePlot(alldata_harmony,
            features = c("BMP4","HAND1","FOXF1"),
            pt.size = 0.001, raster = FALSE,
            cols = c("lightgrey","red"),
            split.by = "study2")


# Al markers 
known_Al_genes <- c("FOXF1", "HAND1", "COL3A1", "COL1A1", "ISL1", "TCF21",
                    "COL6A2","COL6A1","F5","PCOLCE","VCAN")
FeaturePlot(alldata_harmony,
            features = known_Al_genes,
            pt.size = 0.001, raster = FALSE,
            cols = c("lightgrey","red"),
            split.by = "study2")
ggsave(filename = "embryoid_embryoCS8-9_knownAlMarkers_featureplot_split.pdf",
       height = 28.89, width = 6.89)

# ysMeso markers 
known_ysMes_genes <- c("ANXA1", "ANXA8", "COL3A1", "COL1A1", "HAND1", "APOE",
                    "SLC9A3R1","MAF","GNG5","BNC1","EPHA7")
FeaturePlot(alldata_harmony,
            features = known_ysMes_genes,
            pt.size = 0.001, raster = FALSE,
            cols = c("lightgrey","red"),
            split.by = "study2")
ggsave(filename = "embryoid_embryoCS8-9_knownYsMesoMarkers_featureplot_split.pdf",
       height = 28.89, width = 6.89)

#Mes marker
known_Mes_genes <- c("ANXA8", "COL3A1", "COL1A1", 
                       "CSRP2","GPX4","COL21A1","TNNI1","DCN")
FeaturePlot(alldata_harmony,
            features = known_Mes_genes,
            pt.size = 0.001, raster = FALSE,
            cols = c("lightgrey","red"),
            split.by = "study2")
ggsave(filename = "embryoid_embryoCS8-9_knownMesMarkers_featureplot_split.pdf",
       height = 19.89, width = 6.89)

library(SCpubr)
p1 <- do_DimPlot(alldata_harmony,
           pt.size = 0.001,
           group.by = "annotation",
           idents.keep = c("invivo_AI","invivo_Mes","invivo_ys.Meso1","invivo_ys.Meso2"),
           plot.axes = TRUE,
           raster = FALSE)
p1

