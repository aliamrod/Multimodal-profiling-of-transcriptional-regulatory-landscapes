# Read in table
# unload gsl/2.7
# module load gsl/2.4
# R
# > install.packages("remotes")
# > remotes::install_github("RGLab/MAST")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("SingleCellExperiment")
library(SingleCellExperiment)
library(Seurat)
library(ggplot2)
library(dplyr)
library(config)
library(BiocManager)
library(dplyr)
library(Seurat)
library(ggplot2)
library(MAST)
library(config)

DE_Method='MAST'
top_DE <- 40
resolution <- 0.2

#### Initialize objects ####
merged.data <- readRDS('/work/Neuroinformatics_Core/s204365/outputDirectory/snRNA_analyses/snRNA_Ana_Oct31/merged_scRNA_filtered_umap.RDS')


#### Store Initial cluster Plot ####
pdf('/work/Neuroinformatics_Core/s204365/outputDirectory/snRNA_analyses/snRNA_Ana_Oct31/UMAPforDE.pdf', height=8,width=8)
DimPlot(merged.data, reduction = "umap",label = T,repel=T,pt.size = 0.2)
print(p)
dev.off()


DefaultAssay(merged.data) <- 'RNA'
data.markers <- FindAllMarkers(merged.data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,test.use = DE_Method)
saveRDS(data.markers, '/work/Neuroinformatics_Core/s204365/outputDirectory/snRNA_analyses/snRNA_DE_markers.RDS')

#data.markers %>% group_by(cluster) %>% top_n(2)
top50 <- data.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
top20 <- data.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.table(top50,file='/work/Neuroinformatics_Core/s204365/outputDirectory/snRNA_analyses/snRNA_Ana_Oct31/snRNA_top40_DEmarkers.tsv',quote=F,col.names=T,row.names=F,sep='\t')
DefaultAssay(merged.data) <- 'RNA'
pdf('/work/Neuroinformatics_Core/s204365/outputDirectory/snRNA_analyses/snRNA_Ana_Oct31/DE_top20.pdf', height=40,width=18)
DoHeatmap(merged.data, features = top20$gene,raster = T) + NoLegend()
dev.off()

pdf('/work/Neuroinformatics_Core/s204365/outputDirectory/snRNA_analyses/snRNA_Ana_Oct31/DE_top50.pdf', height=40,width=18)
DoHeatmap(merged.data, features = top50$gene,raster = T) + NoLegend()
dev.off()

pdf('plots/scRNA/CellCycle_VlnPLot.pdf',height=10,width=20)
print(VlnPlot(cortex, c("S.Score","G2M.Score"),pt.size = 0.001))
dev.off()        

ClusterID <- c("radial glia", "intermediate progenitors", "interneurons", "neurons", "UL neurons", "DL neurons", "L1")
cluster_annot <- read.table('/work/Neuroinformatics_Core/s204365/scRNA_0001/marker_genes_dev_ctx.txt',header=T,row.names="ClusterID", sep=":")
  'results/RNA_clusters_unfiltered.tsv',header=T,row.names = "ClusterID")
new.cluster.ids <- as.vector(cluster_annot$CellType)
names(new.cluster.ids) <- as.factor(levels(cortex))
cortex <- RenameIdents(cortex, new.cluster.ids)

cortex<- merged.data
levels(cortex) <- c('NSC','NSC_M','IPC','IPC_M','PN1','PN2','PN3','CR','IN','MG','Mural')

pdf('plots/scRNA/UMAPwithIDs.pdf',height=8,width=8)
DimPlot(cortex, reduction = "umap",label = T,repel=T,pt.size = 0.2) + NoLegend()
dev.off()
saveRDS(cortex,'data/merged_scRNA_unfiltered_IDs.RDS') 

DefaultAssay(cortex) <- 'RNA'

cortex.markers <- FindAllMarkers(cortex, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,test.use = DE_Method)
saveRDS(cortex.markers,'data/scRNA_DE_markers_IDs.RDS')

top40 <- cortex.markers %>% group_by(cluster) %>% top_n(n = 40, wt = avg_logFC)
top20 <- cortex.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
write.table(top40,file='results/scRNA_top40_DEmarkers_IDs.tsv',quote=F,col.names=T,row.names=F,sep='\t')
DefaultAssay(cortex) <- 'RNA'
pdf('plots/scRNA/DE_top20.pdf',height=40,width=20)
DoHeatmap(cortex, features = top20$gene,raster = T) + NoLegend()
dev.off()
pdf('plots/scRNA/DE_top40.pdf',height=40,width=20)
DoHeatmap(cortex, features = top40$gene,raster = T) + NoLegend()
dev.off()


saveRDS(merged.data.filtered,'data/merged_scRNA_filtered_IDs.RDS')  








### References ###
# [1] "MAST: a flexible statistical framework for assessing transcriptional changes and characterizing heterogeneity in single-cell RNA sequencing data" (2015)
# [2] "Divergent transriptional regulation of astrocyte reactivity across disorders" (2022)
# [3] "Differential integrated stress response and asparagine production drive symbiosis and therapy resistance of pancreatic adenocarcinoma cells" (2022)
# [4] "A novel insight into differential expression profiles of sporadic cerebral cavernouse malformation patients with different symptoms" (2021)
# [5] "Differential transcript usage unravels gene expression alterations in Alzheimer's disease human brains" (2021)
# [6] "A clustering-independent method for finding differentially expressed genes in single-cell transcriptome data" (2020)
# [7] "Gene regulation by gonadal hormone receptors underlies brain sex differences" (2022)
# [8] "Differential expression of an endogenous retroviral element [HERV-K(HML-6)] is associated with reduced survival in glioblastoma patients" (2022)
# [9] "Single-cell and bulk transcriptome seuqencing identifies two epithelial tumor cell states and refines the consensus molecular classification of colorectal cancer" (2022)
# [10] "Differential expression analysis of genes and long non-coding RNAs associated with KRAS mutation in colorectal cancer cells" (2022)
# [11] "The landscape of tumor cell states and spatial organization in H3-K27M mutant diffuse midline glioma across age and location" (2022)
# [12] "Molecular biomarkers for embryonic and adult neural stem cell and neurogenesis" (2015)
