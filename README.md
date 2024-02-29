**R scripts to generate panels from the main figures of Ginno et al., 2024**

R version 4.3.0

**Required packages:**
TxDb.Mmusculus.UCSC.mm10.knownGene, BSgenome.Mmusculus.UCSC.mm10, GenomicRanges, RColorBrewer, changepoint,
Repitools, randtests, Rsamtools, apcluster, regioneR, pheatmap, mixtools, ggplot2, viridis, scales, psych
QuasR, Gviz, grid

**Processed files:**
GEO submission GSE230579. 
PF1 mutations: GSE230579_PF1_SisterMutations.rds.gz, 
CAST/C3H liver tumour mutations: GSE230579_F1CastB6_mutations.rds.gz, 
ATAC peaks table in PF1: GSE230579_PF1_ATAC_Peaks.rds.gz, 
RNA genic counts in PF1: GSE230579_PF1_RNA_counts.rds.gz

**BAM file generation:**
Raw sequencing data for PF1 sisters and Mus castaneus x C3H tumours are in SRA submission PRJNA934746. Alignment parameters for generating and filtering bam files from raw sequencing data are described in materials and methods under the 
subheadings "Whole genome and ATAC alignment and filtering", "Dual-hybrid N-masked reference for WGS F1 tumour haplotype discrimination" and "PF1 RNA-seq data processing and analysis".

**Phenomex files:**
Tables of FUCCI fluorophore imaging on the Phenomex platform are included in this repository and necessary for generating Panel 1c.
