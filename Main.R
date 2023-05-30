
library(ggplot2)
library(ggrepel)
library(DESeq2)
library(ComplexHeatmap)
library(ggtext)
library(RColorBrewer)
library(org.Mm.eg.db)
library(enrichplot)
library(ggsci)
library(clusterProfiler)
library(ggupset)
source("Analysis_DGE.R")
source("Quality_Control.R")



dds <- Normalize_Counts(cts = cts, sampleInfo = sampleInfo,Variable_Of_Interest = "treatment", Groups_Selected = c("Control", "TumorBearing"),Folder_Name = "QC")

Distance_Clustering(dds = dds, Groups_Selected = c("Control", "TumorBearing"), Variable_Of_Interest = "treatment", Folder_Name = "QC" )

PCA_Plots(dds = dds,Variable_Of_Interest = "treatment", Samples_Column_name = "id", sampleInfo = sampleInfo, Groups_Selected = c("Control", "TumorBearing") , Variables_For_PCA = c("genotype","sex","type","RIN:numeric","mapping_percentage"), Folder_Name = "QC",Color_Choice = c("Blue","black"), Shape_Choice = c(12, 13, 14, 15))

Y_Reads(cts = cts, geneInfo = geneInfo, sampleInfo = sampleInfo, Folder_Name = "QC")

X_Reads(cts = cts, geneInfo = geneInfo, sampleInfo = sampleInfo, Folder_Name = "QC")

XIST_Counts(cts = cts, geneInfo = geneInfo, genes_column_name = "SYMBOL", sampleInfo = sampleInfo, Folder_Name = "QC")

twoGroupCompare(Feature_Counts = cts,Sample_Info = sampleInfo, Samples_Column_Name = "id", Gene_Info = geneInfo, Genes_Column_Name = "SYMBOL", Variable_Of_Interest = "treatment", Groups_Selected = c("Control","CalciumDeficient"), Covariates = c("sex","genotype"), Folder_name = "ControlVsCalciumDeficient", pvalue_Cutoff = 0.05,log2Fold_Cutoff = log2(1.5))

