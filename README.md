# Multimodal-profiling-of-transcriptional-regulatory-landscapes

**gif**


Description of Scripts and Analysis.
snRNA-seq Analyses.
snRNA1-QC: Quality control and standard pre-processing of scRNA data using Seurat v3.1.5.
snRNA2-UMAP: Louvain clustering + UMAP visualization of scRNA data following transformation. 
snRNA3-DE: Assignment of cluster identities based DE genes (utilizing MAST- Finak et al., 2015). 
snRNA4-monocle3: Pseudotime analysis using Monocle3. Analyses include model fitting to identify gene expression changes as a function of pseudotime. 

snATAC-seq Analyses.
snATAC1-QC: Quality control of scATAC data based on (1) number of fragments and (2) TSS enrichment. 
snATAC2-aggregateBin: Final QC and code for generating count matrices based on genomic bins (for initial clustering), gene bodies and promoters. 
snATAC3-normalization: Code for initial clustering based on fragments from fixed-size genome wide bins.
snATAC4-peakNormalization: Final peak calling based on initial clusters to generate high-quality peak set, used for final clustering and visualization.
snATAC5-chromVAR_motifs: Computing motif accessibility deviations using chromVAR (Schep et al., 2017) implemented in Signac.
snATAC6-Compute_Gene_Scores: Computing gene activity scores using Cicero, used for subsequent integration analyses (observed below). 
snATAC7-cluster_unique_peaks: Identification of cluster specific peaks.

Integration Analyses of snRNA and snATAC.
snRNA_snATAC_Integration_01_Align_snATAC_snRNA: Integration of snRNA and snATAC data using Seurat CCA and identifcation of nearest neighbors (kNN). 
snRNA_snATAC_Integration_02_Create_Aggreagte_snATAC_snRNA: Aggreagte snRNA and snATAC data using nearest neighbor information.
snRNA_snATAC_Integration_04_P2G_analysis: Further characterization of identified enhancer-gene pairs.
snRNA_snATAC_Integration_05_P2G_monocle: Pseudotime analysis on integrated snRNA-snATAC object using Monocle3. Analyses also include model fitting to identify changes of accessibility and motif deviations as a function of pseudotime.
snRNA_snATAC_Integation_06_chromVAR: Motif analysis for integrated object using chromVAR.

List of Figures.


