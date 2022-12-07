
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


files<- '/work/Neuroinformatics_Core/s204365/scRNA_0001/doubletfinder_seur_objects/'
setwd(files)


# Adjust CellRanger, CellBender, and DoubletFinder values, accordingly.
graph_plot_QC <- read.table(
  header=TRUE, text='Category        SampleID NucleiCount
1   a_CellRanger       AO11      12860
2   a_CellRanger       AO12      11840
3   a_CellRanger       AO19      3285
4   a_CellRanger       AO20      4240
5   a_CellRanger       AO21       4064
6   a_CellRanger       AO22       3551
7   a_CellRanger       AO23      4119
8   a_CellRanger        AO24      3740
9   b_CellBender      AO11       12421
10  b_CellBender      AO12      11598
11  b_CellBender      AO19      2516
12  b_CellBender      AO20      3915
13  b_CellBender      AO21      3887
14  b_CellBender      AO22      3365
15  b_CellBender      AO23      3889
16  b_CellBender      AO24      3346
17 DoubletFinder    AO11      9664
18 DoubletFinder    AO12      9251
19 DoubletFinder    AO19      2346
20 DoubletFinder    AO20      3644
21 DoubletFinder    AO21      3617
22 DoubletFinder    AO22      3133
23 DoubletFinder    AO23      3620
24 DoubletFinder    AO24      3112 
  ')
ggplot(graph_plot_QC, aes(SampleID, NucleiCount, fill = Category)) + 
  geom_bar(stat="identity", position = "dodge") + 
  scale_fill_brewer(palette = "Set1")


#install.packages("stringfish")
#remotes::install_cran("qs", type = "source", configure.args = "--with-simd=AVX2")

batch_correction=FALSE
add_control=FALSE
min.cells <- 3
min.features <- 200
minGenes <- 1000
minRNA <- 2500
maxGenes <- 7000
#maxRNA <- 50000
maxMT <- 10

input_dir<- "/work/Neuroinformatics_Core/s204365/scRNA_0001/"

# Load data into Seurat
ao11<- readRDS("sampleao11.rds")
ao12<- readRDS("sampleao12.rds")
ao19<- readRDS("sampleao19.rds")
ao20<- readRDS("sampleao20.rds")
ao21<- readRDS("sampleao21.rds")
ao22<- readRDS("sampleao22.rds")
ao23<- readRDS("sampleao23.rds")
ao24<- readRDS("sampleao24.rds")

ao11$genotype<- "CKO"
ao12$genotype<- "FLX"
ao19$genotype<- "CKO"
ao20$genotype<- "FLX"
ao21$genotype<- "CKO"
ao22$genotype<- "CKO"
ao23$genotype<- "FLX" #alphabetically organizes figures
ao24$genotype<- "FLX" 

ao11$SampleID<- "AO11"
ao12$SampleID<- "AO12"
ao19$SampleID<- "AO19"
ao20$SampleID<- "AO20"
ao21$SampleID<- "AO21"
ao22$SampleID<- "AO22"
ao23$SampleID<- "AO23"
ao24$SampleID<- "AO24"

ao11$library.batch<- "1"
ao12$library.batch<- "1"
ao19$library.batch<- "2"
ao20$library.batch<- "2"
ao21$library.batch<- "2"
ao22$library.batch<- "2"
ao23$library.batch<- "2"
ao24$library.batch<- "2"

### QC Plots ###
# Analyze median QC values per sample/library.
# Using gene name patterns

### Merge Seurat object ###
merged.data<- merge(ao11, y = c(ao12, ao19, ao20, ao21, ao22, ao23, ao24),
                    add.cell.ids=c("AO11", "AO12",
                                   "AO19", "AO20",
                                   "AO21", "AO22", 
                                   "AO23", "AO24"))


### Add Mitochondrial and Ribosomal Gene Percentages ###
merged.data<- Add_Mito_Ribo_Seurat(seurat_object = merged.data, species = "Mouse")
merged.data<- Add_Cell_Complexity_Seurat(seurat_object = merged.data, overwrite = TRUE)

merged.data@meta.data$group[merged.data@meta.data$SampleID == "AO11"|
                              merged.data@meta.data$SampleID == "AO19" |
                              merged.data@meta.data$SampleID == "AO21"|
                              merged.data@meta.data$SampleID == "AO22"]<- "CKO"
merged.data@meta.data$group[merged.data@meta.data$SampleID == "AO12"|
                              merged.data@meta.data$SampleID == "AO20" |
                              merged.data@meta.data$SampleID == "AO23" |
                              merged.data@meta.data$SampleID == "AO24"]<- "FLX"
merged.data<- Add_Mito_Ribo_Seurat(seurat_object = merged.data, species = "Mouse", overwrite=TRUE)


#Grouped by: CKO vs. FLX groups.
p1<- QC_Plots_Genes(seurat_object = merged.data, low_cutoff = 300,high_cutoff = 8000, group.by="group", pt.size=0)
p2<- QC_Plots_UMIs(seurat_object = merged.data, low_cutoff = 1000, high_cutoff = 55000,group.by="group",pt.size=0)
p3<- QC_Plots_Mito(seurat_object = merged.data, high_cutoff = 2.5, group.by="group",pt.size=0, plot_title="Mitochondrial Proportion (mtDNA%)")
p4<- QC_Plots_Complexity(seurat_object = merged.data, high_cutoff = 0.8, group.by="group",pt.size=0)
wrap_plots(p1,p2,p3,p4, ncol =4)


#Grouped by: Samples (8).
p5<- QC_Plots_Genes(seurat_object = merged.data, low_cutoff = 300,high_cutoff = 8000, group.by="SampleID", pt.size=0)
p6<- QC_Plots_UMIs(seurat_object = merged.data, low_cutoff = 1000, high_cutoff = 55000,group.by="SampleID",pt.size=0)
p7<- QC_Plots_Mito(seurat_object = merged.data, high_cutoff = 2.5, group.by="SampleID", pt.size=0,plot_title="Mitochondrial Proportion (mtDNA%)")
p8<- QC_Plots_Complexity(seurat_object = merged.data, high_cutoff = 0.8, group.by="SampleID",pt.size=0)
wrap_plots(p5,p6,p7,p8, ncol =4)


# Feature scatter plot.
QC_Plot_UMIvsGene(seurat_object = merged.data, low_cutoff_gene = 200, high_cutoff_gene = 8000, low_cutoff_UMI = 1000, group.by="SampleID")
QC_Plot_GenevsFeature(seurat_object = merged.data, feature1 = "percent_mito",low_cutoff_gene = 300,
                      high_cutoff_gene = 10000, high_cutoff_feature = 20,group.by="SampleID")

# Calculate median values and return data.frame
median_stats<- Median_Stats(seurat_object = merged.data, group_by_var="genotype",
                            median_var = "genotype")
Plot_Median_Genes(seurat_object = merged.data, group_by = "genotype")

Plot_Median_Genes(seurat_object = merged.data, group_by = "genotype")
Plot_Median_UMIs(seurat_object = merged.data, group_by = "SampleID")
Plot_Median_Mito(seurat_object = merged.data, group_by = "SampleID")
Plot_Median_Other(seurat_object = merged.data, median_var = "percent_ribo", group_by = "SampleID")



# apply filter
merged.data <- subset(merged.data, nCount_RNA >= 200 & nCount_RNA <= 50000 & mitoPercent <= 15)


# plot the number of cells in each sample post filtering
df <- as.data.frame(rev(table(merged.data$SampleID)))
colnames(df) <- c('SampleID', 'n_cells')
p <- ggplot(df, aes(y=n_cells, x=reorder(SampleID, -n_cells), fill=SampleID)) +
  geom_bar(stat='identity') +
  scale_y_continuous(expand = c(0,0)) +
  NoLegend() + RotatedAxis() +
  ylab(expression(italic(N)[cells])) + xlab('Sample ID') +
  ggtitle(paste('Total cells post-filtering:', sum(df$n_cells))) +
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major.y=element_line(colour="lightgray", size=0.5),
  )
