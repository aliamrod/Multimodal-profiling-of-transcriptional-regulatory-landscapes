### Load libraries. ###
devtools::install_github(repo = "samuel-marsh/scCustomize", force = TRUE)
library(scCustomize)
library(dplyr)
library(Seurat)
library(Hmisc)
library(SummarizedExperiment)
require(edgeR)
library(tidyverse)
library(knitr)
library(devtools)
library(ggpubr)
library(qs)
library(patchwork)
require(LSD)
library(Matrix)
library(cowplot)
library(dplyr)
library(patchwork)
theme_set(theme_cowplot())

install.packages("config")
library(config)
source("/project/Neuroinformatics_Core/Konopka_lab/s204365/MACS2/plot_functions.R")
source("/project/Neuroinformatics_Core/Konopka_lab/s204365/MACS2/aux_functions.R")

############################################################
# Visualization
library(cowplot)
theme_set(theme_cowplot())
library(RColorBrewer)
library(viridis)
library(dplyr)
library(Seurat)
library(ggplot2)
library(harmony)
library(config)
set.seed(1)

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
  panel.grid.minor=element_blank(),
  plot.background=element_blank()
)
nDims <- c(20,25,30,35)
resolutions <- c(0.1, 0.3, 0.8, 1, 1.2)
umap_methods <- 1


#### Normal transformation + Log normalization ####
DefaultAssay(merged.data) <- 'RNA'
merged.data<- NormalizeData(merged.data, 
                            normalization.method = "LogNormalize",
                            scale.factor = 10000)
merged.data <- FindVariableFeatures(merged.data, selection.method = "vst",nfeatures=4000)
top50 <- head(VariableFeatures(merged.data), 50)
plot1 <- LabelPoints(
  VariableFeaturePlot(merged.data),
  points = head(VariableFeatures(merged.data),10),
  repel = TRUE
) + theme(legend.position="bottom")
png('/work/Neuroinformatics_Core/s204365/outputDirectory/snRNA_analyses/snRNA_Ana_Oct31//variable_features.png', width=5, height=5, res=200, units='in')
print(plot1)
dev.off()

merged.data <- ScaleData(merged.data,features =rownames(merged.data))
merged.data <- RunPCA(merged.data,npcs = 50)
pdf('/work/Neuroinformatics_Core/s204365/outputDirectory/snRNA_analyses/PCA.pdf',height=24,width=1
    
# PCA scatter plot colored by SampleID
pca2 <- DimPlot(merged.data, reduction = "pca", group.by='SampleID')
    
DimHeatmap(merged.data, dims = 1:20, cells = 500, balanced = TRUE)
DimHeatmap(merged.data, dims = 20:40, cells = 500, balanced = TRUE)
ElbowPlot(merged.data,ndims = 50)
dev.off()

write.table(merged.data[["pca"]]@feature.loadings,'/work/Neuroinformatics_Core/s204365/outputDirectory/snRNA_analyses/snRNA_Ana_Oct31/PCA.txt',quote=F,row.names=T,col.names=T,sep='\t')
merged.data <- FindNeighbors(merged.data, dims = 1:20,force.recalc = T)
merged.data <- FindClusters(merged.data, resolution = 0.6,n.start = 20,n.iter = 30,algorithm = 1)
merged.data <- RunUMAP(merged.data, dims = 1:20,min.dist = 0.5,spread = 1.5,n.components = 2L)
merged.data <- JackStraw(merged.data, num.replicate = 100)
merged.data <- ScoreJackStraw(merged.data, dims = 1:20)
JackStrawPlot(merged.data, dims = 1:15)

options(repr.plot.height = 2.5, repr.plot.width = 6)
merged.data <- merged.data %>% 
RunHarmony("library.batch", plot_convergence = TRUE)
merged.data <- RunUMAP(merged.data, reduction='harmony', dims = 1:30)
merged.data <- FindClusters(merged.data,resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0))

pdf('/work/Neuroinformatics_Core/s204365/outputDirectory/snRNA_analyses/snRNA_Ana_Oct31/', height=8,width=16)

DimPlot(merged.data, reduction = "umap",label = T,group.by='SampleID',cols=length(unique(merged.data@meta.data$SampleID)))
DimPlot(merged.data, reduction = "umap",label = T,group.by='genotype',cols=length(unique(merged.data@meta.data$genotype)))
dev.off()

# set up list of canonical cell type markers
canonical_markers <- list(
  'Radial glia' = c("Pax6", "Mki67", "Emx2"),
  'Intermediate progenitors' = c("Eomes" ),
  'Interneurons' = c("Gad2","Gad1","Sp9","Npy","Sst","Cck"),
  'Neurons' = c("Satb2", "Mef2c", "Rbfox3"),
  'UL neurons' = c("Cux1","Cux2", "Satb2", "Rorb"),
  'DL neurons' = c("Bcl11b", "Foxp2","Tle4","Tbr1"),
  'L1' = c('Reln')
)

png('/work/Neuroinformatics_Core/s204365/featureplot.png', width=10, height=10, units='in', res=200)
DoHeatmap(merged.data, group.by ="seurat_clusters", features=as.character(unlist(canonical_markers))) # + scale_fill_gradientn(colors=viridis(256))
dev.off()




######

######### DATA Exploration ####################
normalizations=c('harmony', 'VST')
nDims = c(20,25,30,25)
resolutions<- c(0.1,0.2,0.3,0.8, 1.0, 1.2)
umap_methods<- 1


  for (nDim in nDims){
    merged.data <- RunUMAP(merged.data, dims = 1:nDim,min.dist = 0.5,spread = 1.5,n.epochs=500)
    merged.data <- FindNeighbors(merged.data, dims = 1:nDim)
    for (resolution in resolutions){
      for (umap_method in umap_methods){
        merged.data <- FindClusters(merged.data, resolution = resolution,algorithm = umap_method)
        pdf(paste0('plots/scRNA/UMAP_',normalization,'_dims',nDim,'_res',resolution,'_umap',umap_method,'.pdf'),height=8,width=16)
        print(DimPlot(merged.data, reduction = "umap",label = T,group.by='orig.ident',cols=sample_colors(length(unique(merged.data$orig.ident)))))
        print(DimPlot(merged.data, reduction = "umap",label = T,repel=T,cols=cluster_colors(length(unique(Idents(merged.data))))))
        dev.off()
        for (set in c('diff_markers','layer_markers','temp_markers','other_markers','regional_markers')){
          pdf(paste0('plots/scRNA/UMAP_',normalization,'_dims',nDim,'_res',resolution,'_umap',umap_method,'_',set,'.pdf'),height=16,width=16)
          print(FeaturePlot(merged.data, features = get(set),cols = gene_colors(3),pt.size = 0.1,label = T))
          dev.off()
        }    
        pdf(paste0('/work/Neuroinformatics_Core/s204365/outputDirectory/snRNA_analyses/snRNA_Ana_Oct31/UMAP_','_dims',nDim,'_res',resolution,'_umap',umap_method,'_Features3.pdf'),height=16,width=16)
        print(VlnPlot(merged.data,'nFeature_RNA',pt.size = 0.001))
        print(VlnPlot(merged.data,'nCount_RNA',pt.size = 0.001))
        print(VlnPlot(cortex,'percent.mt',pt.size = 0.001))
        dev.off()
      }
    }
  }
}


########################################################
####### Choose values for subsequent analysis ##########
########################################################

normalization <- 'VST'
cluster_method <- 1
resolution <- 0.3
nDim <- 20
doHarmony=TRUE


DimPlot(merged.data, reduction = "umap",label = T,group.by='orig.ident',cols=c('darkgreen','darkred'))
DimPlot(merged.data, reduction = "umap",label = T,repel=T,cols=length(levels(merged.data)))
dev.off()

pdf(paste0('plots/scRNA/chosen_',normalization,'_dims',nDim,'_res',resolution,'_umap',umap_method,'_Features3.pdf'),height=16,width=16)
p1 <- VlnPlot(cortex,'nFeature_RNA',pt.size = 0.001)
p2 <- VlnPlot(cortex,'nCount_RNA',pt.size = 0.001)
p3 <- VlnPlot(cortex,'percent.mt',pt.size = 0.001)
print(CombinePlots(list(p1,p2,p3)))
dev.off()


# create feature plots, cutoff expression values for the 98th and 99th percentile
plot_list <- FeaturePlot(
  merged.data,
  features=unlist(canonical_markers),
  combine=FALSE, cols=viridis(256),
  max.cutoff='q98'
)

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
  panel.grid.minor=element_blank(),
  plot.background=element_blank()
)

# apply theme to each feature plot
for(i in 1:length(plot_list)){
  plot_list[[i]] <- plot_list[[i]] + umap_theme + NoLegend()
}

png('/work/Neuroinformatics_Core/s204365/basic_canonical_marker_featurePlot.png', width=10, height=10, units='in', res=200)
CombinePlots(plot_list)
dev.off()

library(BiocManager)
BiocManager::install('Nebulosa')
library(Nebulosa)

for(celltype in names(canonical_markers)){
  
  print(celltype)
  cur_features <- canonical_markers[[celltype]]
  
  # plot distributions for marker genes:
  p <- VlnPlot(
    merged.data,
    group.by='seurat_clusters',
    features=cur_features,
    pt.size = 0, ncol=1
  )
  png(paste0('/work/Neuroinformatics_Core/s204365/basic_canonical_marker_',celltype,'_vlnPlot.png'), width=10, height=3*length(cur_features), units='in', res=200)
  print(p)
  dev.off()
  
}


gene_list_plot <- c("Pax6", "Mki67", "Emx2", "Eomes",
                    "Gad2","Gad1","Sp9","Npy","Sst","Cck",
                    "Satb2", "Mef2c", "Rbfox3",
                    "Cux1","Cux2", "Satb2", "Rorb",
                    "Bcl11b", "Foxp2","Tle4","Tbr1","Reln")


human_colors_list <- c("dodgerblue", "navy", "forestgreen", "darkorange2", "darkorchid3", "orchid",
                       "orange", "gold", "gray", "red", "blue", "green", "pink")

# Create Plots
png('/work/Neuroinformatics_Core/s204365/stackedviolinplot.png')

Stacked_VlnPlot(seurat_object = merged.data, features = gene_list_plot, x_lab_rotate = TRUE,
                colors_use = human_colors_list)
dev.off()
sample_colors <- c("dodgerblue", "forestgreen", "firebrick1")

# Create Plots
Stacked_VlnPlot(seurat_object = merged.data, features = gene_list_plot, x_lab_rotate = TRUE,
                colors_use = sample_colors, split.by = "genotype")


VlnPlot_scCustom(seurat_object = merged.data, features = "Foxp2", raster = FALSE)


DotPlot_scCustom(seurat_object = merged.data, features = gene_list_plot[1:15], flip_axes = T,
                 remove_axis_titles = FALSE)


all_markers <- FindAllMarkers(object = merged.data)

top5_markers <- Extract_Top_Markers(marker_dataframe = all_markers, num_genes = 5, named_vector = FALSE,
                                    make_unique = TRUE)

Clustered_DotPlot(seurat_object = merged.data, features = top5_markers)

cluster_markers<- all_markers
##########
cluster_annotations <- list(
  '0' = 'UL',
  '1' = 'DL',
  '2' = 'UL',
  '3' = 'DL',
  '4' = 'Interneurons (Inhib)',
  '5' = 'UL',
  '6' = 'Radial',
  '7' = 'Interm. Progenitors',
  '8' = 'Radial',
  '9' = 'UL',
  '10' = 'L1',
  '11' = 'UL',
  '12' = 'UL'
)

# add CellType to seurat metadata
merged.data$CellType <- unlist(cluster_annotations[merged.data$seurat_clusters])
merged.data$CellType_cluster <- paste0(merged.data$CellType, '-', merged.data$seurat_clusters)

png('figures/basic_umap_celltypes.png', width=8, height=7, res=200, units='in')
DimPlot(merged.data, reduction = "umap", group.by='CellType') +
  umap_theme + ggtitle('UMAP colored by cell type annotations')
dev.off()

png('figures/basic_umap_celltype_clusters.png', width=8, height=8, res=200, units='in')
DimPlot(merged.data, reduction = "umap", group.by='CellType_cluster', label=TRUE) +
  umap_theme + ggtitle('UMAP colored by cell type + cluster') + NoLegend()
dev.off()


# plot the number of DEGs per cluster:
df <- as.data.frame(rev(table(cluster_markers$CellType_cluster)))
colnames(df) <- c('cluster', 'n_DEGs')

# bar plot of the number of cells in each sample
p <- ggplot(df, aes(y=n_DEGs, x=reorder(cluster, -n_DEGs), fill=cluster)) +
  geom_bar(stat='identity') +
  scale_y_continuous(expand = c(0,0)) +
  NoLegend() + RotatedAxis() +
  ylab(expression(italic(N)[DEGs])) + xlab('') +
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major.y=element_line(colour="lightgray", size=0.5),
  )

png('figures/basic_DEGs_barplot.png', width=9, height=4, res=300, units='in')
print(p)
dev.off()

# plot the top 3 DEGs per cluster as a heatmap:
top_DEGs <- cluster_markers %>%
  group_by(CellType_cluster) %>%
  top_n(3, wt=avg_logFC) %>%
  .$gene

png('figures/basic_DEGs_heatmap.png', width=10, height=10, res=300, units='in')
pdf('figures/basic_DEGs_heatmap.pdf', width=15, height=12, useDingbats=FALSE)
DoHeatmap(seurat_obj, features=top_DEGs, group.by='seurat_clusters', label=FALSE) + scale_fill_gradientn(colors=viridis(256)) + NoLegend()
dev.off()


#########
cluster_markers$CellType_cluster <- paste0(unlist(cluster_annotations[cluster_markers$cluster]), '-', cluster_markers$cluster)

df <- as.data.frame(rev(table(cluster_markers$CellType_cluster)))
colnames(df) <- c('cluster', 'n_DEGs')

# bar plot of the number of cells in each sample
p <- ggplot(df, aes(y=n_DEGs, x=reorder(cluster, -n_DEGs), fill=cluster)) +
  geom_bar(stat='identity') +
  scale_y_continuous(expand = c(0,0)) +
  NoLegend() + RotatedAxis() +
  ylab(expression(italic(N)[DEGs])) + xlab('') +
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major.y=element_line(colour="lightgray", size=0.5),
  )

png('figures/basic_DEGs_barplot.png', width=9, height=4, res=300, units='in')
print(p)
dev.off()

# plot the top 3 DEGs per cluster as a heatmap:
top_DEGs <- cluster_markers %>%
  group_by(CellType_cluster) %>%
  top_n(3, wt=avg_logFC) %>%
  .$gene

png('figures/basic_DEGs_heatmap.png', width=10, height=10, res=300, units='in')
pdf('figures/basic_DEGs_heatmap.pdf', width=15, height=12, useDingbats=FALSE)
DoHeatmap(seurat_obj, features=top_DEGs, group.by='seurat_clusters', label=FALSE) + scale_fill_gradientn(colors=viridis(256)) + NoLegend()
dev.off()



### Save final merged, filtered file. ###
saveRDS(merged.data,file='/work/Neuroinformatics_Core/s204365/outputDirectory/snRNA_analyses/merged_scRNA_filtered_umap.RDS')
