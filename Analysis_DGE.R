#function formatDouble
#Arguments: x and digit
#Converts the value in x to have only two values after the decimal
formatDouble <- function(x, digit = 2){
  format(round(x, digit), nsmall = digit)
}


#function plotCt
#Arguments: dds, res, idxGene, variable and Samples column name
#Generates a plot of counts for the given gene in idxGene
plotCt <- function(dds, res, idxGene,variable,Samples_Column_Name){
  ct <- counts(dds, normalized = TRUE, replaced = FALSE)[idxGene, ]
  gn <- rownames(dds)[idxGene]
  dplot <- data.frame(
    sample = dds[[Samples_Column_Name]],
    ct = ct,
    group = dds[[variable]],
    stringsAsFactors = FALSE
  )
  
  pos <- position_jitter(width = 0.05, height = 0, seed = 2)
  
  gp <- ggplot(dplot, aes(x=group, y=ct)) + 
    geom_violin()+
    #geom_boxplot(outlier.shape = NA) + 
    geom_jitter(position = pos) + 
    ggrepel::geom_text_repel(aes(label = sample), position = pos) + 
    labs(x="", y="Normalized Counts", 
         title = paste0(
           gn, " (Padj = ", formatDouble(res$padj[idxGene]),  
           "; Log2Fold = ", formatDouble(res$log2FoldChange[idxGene]), ")")) + 
    scale_y_log10() +
    theme_light() + 
    theme(plot.title = element_text(hjust = 0.5, size = 10))
  
  gp
}


#function plotCtCon
#Arguments: dds, res, idxGene, variable and Samples column name
#Generates a correlation plot of normalized counts and given continuous variable(variable) for the given gene in idxGene
plotCtCon <- function(dds, res, idxGene,variable,Samples_Column_Name){
  ct <- counts(dds, normalized = TRUE, replaced = FALSE)[idxGene, ]
  sample = dds[[Samples_Column_Name]]
  gn <- rownames(dds)[idxGene]
  group <- dds[[variable]]
  df <- data.frame(x = group, y = ct )
  gp <- ggplot(df, aes(x = group, y = ct)) +
    ggrepel::geom_text_repel(aes(label = sample)) +
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "red") +
    theme_light() +
    labs(title = paste0(
      gn, " (Padj = ", formatDouble(res$padj[idxGene]),  
      "; Correlation Coefficient = ", cor(group,ct), ")"), x = "Age", y = "Normalized Counts") +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(gp)
}


#function plotHeat
#Arguments: vsd, res, geneInfoSel, pvCutoff, log2FoldCutoff, scale, column_km, cluster_columns and variable
#Extracts the significant genes based on pvalue and log2fold change (pvCutoff and log2FoldCutoff) and generates the heatmap of expression levels of these genes
plotHeat <- function(vsd, res, geneInfoSel, pvCutoff = 0.05, log2FoldCutoff = 1, scale = TRUE, column_km = 1, cluster_columns = TRUE,variable){
  if(sum(rownames(vsd) == rownames(res)) != nrow(vsd)){
    stop("Not Match")
  }
  idxSig <- which(res$pvalue < pvCutoff & abs(res$log2FoldChange) > log2FoldCutoff)
  ct <- assay(vsd)[idxSig, ]
  
  if(scale){
    ct <- t(scale(t(ct)))
  }
  
  idx <- match(rownames(res)[idxSig], geneInfoSel$ENSEMBL)
  gn <- geneInfoSel$SYMBOL[idx]
  
  rownames(ct) <- gn
  
  gCol <- pal_nejm("default")(length(unique(vsd[[variable]])))[seq(length(unique(vsd[[variable]])), 1, -1)]
  if(is.null(levels(vsd[[variable]]))){
    tmp <- factor(vsd[[variable]])
    names(gCol) <- levels(tmp)
  } else {
    names(gCol) <- levels(vsd[[variable]])
  }
  
  ha = HeatmapAnnotation(
    genotype = vsd[[variable]],
    col = list(genotype = gCol),
    annotation_legend_param = list(
      #genotype = list(nrow = 1, direction = "horizontal")
      title = ""
    ),
    show_annotation_name = FALSE
  )
  
  ht <- Heatmap(
    ct, 
    top_annotation = ha,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 4),
    
    column_names_gp = gpar(
      #col = gCol,
      col = "black",
      fontsize = 8),
    
    col = colorRampPalette(rev(brewer.pal(8, "PiYG")))(25),
    column_names_side = "bottom",
    row_names_side = "right",
    
    cluster_columns = cluster_columns,
    column_km = column_km, 
    
    heatmap_legend_param = list(
      title = "Scaled Gene Expression Level", at = c(-2, 0, 2), 
      direction = "horizontal",
      title_position = "topcenter",
      legend_width = unit(6, "cm")
    ),
  )
  
  ht
}


#function plotHeatCon
#Arguments: vsd, res, geneInfoSel, pvCutoff, scale, column_km, cluster_columns and variable
#Extracts the significant genes based on pvalue(pvCutoff) and generates the heatmap of expression levels of these genes
plotHeatCon <- function(vsd, res, geneInfoSel, pvCutoff = 0.05, scale = TRUE, column_km = 1, cluster_columns = TRUE,variable, Genes_Column_Name){
  
  if(sum(rownames(vsd) == rownames(res)) != nrow(vsd)){
    stop("Not Match")
  }
  idxSig <- which(res$pvalue < pvCutoff)
  ct <- assay(vsd)[idxSig, ]
  
  if(scale){
    ct <- t(scale(t(ct)))
  }
  
  idx <- match(rownames(res)[idxSig], geneInfoSel[[Genes_Column_Name]])
  
  gn <- geneInfoSel[[Genes_Column_Name]][idx]
  
  rownames(ct) <- gn
  
  ha = HeatmapAnnotation(
    genotype = vsd[[variable]],
    annotation_legend_param = list(
      #genotype = list(nrow = 1, direction = "horizontal")
      title = toupper(variable)
    ),
    show_annotation_name = FALSE
  )
  
  ht <- Heatmap(
    ct, 
    top_annotation = ha,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 4),
    
    column_names_gp = gpar(
      #col = gCol,
      col = "black",
      fontsize = 8),
    
    col = colorRampPalette(rev(brewer.pal(8, "PiYG")))(25),
    column_names_side = "bottom",
    row_names_side = "right",
    
    cluster_columns = cluster_columns,
    column_km = column_km, 
    
    heatmap_legend_param = list(
      title = "Scaled Gene Expression Level", at = c(-2, 0, 2), 
      direction = "horizontal",
      title_position = "topcenter",
      legend_width = unit(6, "cm")
    ),
  )
  
  ht
}




#function volcanoPlot
#Arguments: dplot, pvCutoff, log2FoldCutoff and  variable
#Generates a volcano plot
volcanoPlot <- function(dplot, pvCutoff = 0.05, log2FoldCutoff = 1,variable){
  tmp <- RColorBrewer::brewer.pal(9,"Set1")
  col <- c(tmp[9], tmp[2], tmp[3], tmp[1])
  names(col) <- c("G0", "G1", "G2", "G3")
  
  dplot[[variable]] <- "G0"
  dplot[[variable]][abs(dplot$fc) > log2FoldCutoff] <- "G1"
  dplot[[variable]][dplot$pv > -log10(pvCutoff)] <- "G2"
  dplot[[variable]][dplot$pv > -log10(pvCutoff) & abs(dplot$fc) > log2FoldCutoff] <- "G3"
  dplot[[variable]] <- factor(dplot[[variable]])
  
  dText <- dplot[dplot[[variable]]=="G3",]
  
  gp <- ggplot(dplot) +
    geom_point(aes(x=fc, y=pv, color = variable), size = 0.8) + 
    geom_text_repel(data = dText, aes(x=fc, y=pv, label = gene), color = col[4], size = 2.5) + 
    labs(x= "Log<sub>2</sub> fold change", y = "-Log<sub>10</sub> <i> P</i>") + 
    scale_color_manual(values = col, labels = c("NS", expression("Log"[2]*" FC"), 'p-value', 
                                                expression("p-value and log"[2]*" FC"))) + 
    theme_light() + 
    theme(axis.title.x = element_markdown(), axis.title.y = element_markdown(), 
          legend.title = element_blank(), legend.position = "top")
  gp
}


#function twoGroupCompare
#Arguments: Feature_Counts, Sample_Info, Samples_Column_Name, Gene_Info, Genes_Column_Name, Variale_Of_Interest, Groups_Selected, Covariates, Folder_name, pvCutoff, log2Fold_Cutoff, Extra_Filters
#Performs the DESeq2 analysis using the counts data(Feature_Counts), Sample Information (Sample_Info), Gene information (Gene_Info) using a categorical variable
#Generates the plots (Sample_Distance, PCA plot, MA plot, Volcano plot, Counts plot of top 10 genes, Heatmap of significant genes)
#Performs Go and KEGG pathway analysis and generates the result files
twoGroupCompare <- function(Feature_Counts,Sample_Info, Samples_Column_Name, Gene_Info, Genes_Column_Name, Variable_Of_Interest, Groups_Selected, Covariates = NULL , Folder_name, pvalue_Cutoff, log2Fold_Cutoff, Extra_Filters = NULL){
  
  

  #Check if Sample IDs are in same order in Feature counts and Sample Info files
  SampleIDs<-colnames(Feature_Counts)
  SampleIDs_SampleInfo<-Sample_Info[[Samples_Column_Name]]
  if(!identical(SampleIDs,SampleIDs_SampleInfo)){
    Feature_Counts<-Feature_Counts[, match(SampleIDs_SampleInfo,colnames(Feature_Counts))]
  }
  
  #Check if Gene IDs are  in same order in Feature counts and Gene Info files
  genes_geneinfo<-sapply(Gene_Info[[Genes_Column_Name]], function(x) sub("\\..*", "", x))
  genes_geneinfo<-data.frame(genes_geneinfo)
  genes_geneinfo<-genes_geneinfo[, 1]
  genes<-sapply(rownames(Feature_Counts), function(x) sub("\\..*", "", x))
  genes<-data.frame(genes)
  genes<-genes[, 1]
  if(!identical(genes_geneinfo,genes)){
    Feature_Counts<-Feature_Counts[match(genes_geneinfo,rownames(Feature_Counts)), ]
  }
  
  #Filters the samples based on given condition
  #Format: (column_name1:column1_values;column_name2:column2_values)
  if(!is.null(Extra_Filters)){
    filter_columns <- unlist(strsplit(Extra_Filters,";"))
    for (i in 1:length(filter_columns)){
      filter_column <- unlist(strsplit(filter_columns[i],":"))[1]
      selected_rows <- unlist(strsplit(filter_columns[i],":"))[2]
      selected_rows <- unlist(strsplit(selected_rows,","))
      for (j in 1:length(selected_rows)){
        Sample_Info <- Sample_Info[Sample_Info[[filter_column]] == selected_rows[j], ]
      }
    }
    
    #Extract the samples in feature counts based on filtered sample info
    sample_columns <- Sample_Info[[Samples_Column_Name]]
    Feature_Counts <- Feature_Counts[, intersect(sample_columns, colnames(Feature_Counts))]
    
  }
  
  #Extract the samples as given in Groups Selected variable.
  idx <- Sample_Info[[Variable_Of_Interest]] %in% Groups_Selected
  Sample_Info_Sel <- Sample_Info[idx,]
  Feature_Counts_Selected <- Feature_Counts[,idx]
  
  #Eliminate those genes that have sum of counts zero across all the samples
  idx <- rowSums(Feature_Counts_Selected == 0) < ncol(Feature_Counts_Selected)
  Gene_Info_Selected <- Gene_Info[idx,]
  Feature_Counts_Selected <- Feature_Counts_Selected[idx,]
  
  #Stop the analysis if the gene names in Gene_Info and Feature counts does not match
  if(sum(Gene_Info_Selected[[Genes_Column_Name]] == rownames(Feature_Counts_Selected), na.rm = TRUE) != nrow(Feature_Counts_Selected)){
    stop("Not Match")
  }
  
  
  #Conversion of categorical variables into factors
  Sample_Info_Sel[[Variable_Of_Interest]] <- factor(Sample_Info_Sel[[Variable_Of_Interest]], levels = Groups_Selected)
  
  
  
  #Design the formula for DESeq2 data object based on given inputs
  if(is.null(Covariates)){
    formula<-as.formula(paste("~ ", Variable_Of_Interest))
  }else{
    formula<-paste("~ ", Variable_Of_Interest)
    for (i in 1:length(Covariates)){
      formula<-paste(formula," + ", Covariates[i])
    }
    formula<-as.formula(formula)
  }
  
  #Construction of the DESeq2 data object
  dds <- DESeqDataSetFromMatrix(
    countData = Feature_Counts_Selected,
    colData = Sample_Info_Sel,
    design= formula)
  
  #Execute the DESeq2 function
  dds <- DESeq(dds)
  
  #resultsNames(dds) # lists the coefficients
  #Extract the results using the values in Groups Selected as arguments for name parameter of results function
  res <- results(dds, name= paste0(Variable_Of_Interest,"_", Groups_Selected[2], "_vs_", Groups_Selected[1]), alpha = 0.1)
  
  #Order the extracted results
  resOrdered <- res[order(res$pvalue),]
  
  #Eliminate those rows which have padj as NA
  resOrdered <- resOrdered[complete.cases(resOrdered$padj), ]
  
  
  # results
  ## dist clustering
  vsd <- vst(dds, blind=TRUE)
  sampleDists <- dist(t(assay(vsd)))
  
  sampleDistMatrix <- as.matrix(sampleDists)
  
  colGenotype <- c(ggsci::pal_nejm("default")(2)[2], ggsci::pal_nejm("default")(2)[1])
  names(colGenotype) <- Groups_Selected
  ha = HeatmapAnnotation(
    genotype = vsd[[Variable_Of_Interest]],
    col = list(genotype = colGenotype),
    annotation_legend_param = list(
      genotype = list(nrow = 1, direction = "horizontal", title = "")
    ),
    show_annotation_name = FALSE
  )
  
  ht <- Heatmap(
    sampleDistMatrix, 
    top_annotation = ha,
    cluster_rows = hclust(as.dist(sampleDistMatrix)),
    cluster_columns = hclust(as.dist(sampleDistMatrix)),
    name = "Dist", 
    show_row_names = FALSE,
    show_column_names = TRUE,
    column_names_gp = gpar(
      col = ifelse(vsd[[Variable_Of_Interest]] == Groups_Selected[1], ggsci::pal_nejm("default")(2)[2], ggsci::pal_nejm("default")(2)[1])),
    
    col = colorRampPalette(rev(brewer.pal(9, "Blues")) )(255),
    column_names_side = "bottom"
  )
  
  dir.create(Folder_name)
  pdf(paste0(Folder_name,"/sampleDist.pdf"), width = 5, height = 6)
  draw(ht,annotation_legend_side = "bottom")
  invisible(dev.off())
  
  
  # PCA
  pcaData <- plotPCA(vsd, intgroup=Variable_Of_Interest, returnData=TRUE)
  
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  gp <- ggplot(pcaData, aes(PC1, PC2, color=group, label = name)) +
    geom_point(size=2, show.legend = FALSE) +
    ggrepel::geom_text_repel(key_glyph = "rect", size = 3) + 
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    scale_color_manual(values = ggsci::pal_nejm("default")(2)[c(2,1)]) + 
    #ggsci::scale_color_nejm() + 
    #coord_fixed() +
    theme_light() +
    guides(color = guide_legend(byrow = T, keyheight = unit(.1, "cm"))) + 
    theme(legend.title = element_blank(), legend.position = "bottom")
  
  pdf(paste0(Folder_name, "/pcaVSD.pdf"), width = 4, height = 4.2)
  print(gp)
  invisible(dev.off())
  
  # MA plot
  pdf(paste0(Folder_name, "/MAPlot.pdf"), width = 5, height = 4)
  plotMA(res, colSig = "red")
  invisible(dev.off())
  
  
  # top 10 genes
  sigGN <- rownames(resOrdered)[1:10]
  gps <- list()
  i <- 1
  for(s in sigGN){
    idx <- which(Gene_Info_Selected[[Genes_Column_Name]] == s)
    gps[[i]] <- plotCt(dds, res, idx,Variable_Of_Interest, Samples_Column_Name)
    i <- i + 1
  }
  
  gp <- ggpubr::ggarrange(plotlist = gps, ncol = 2, nrow = 5)
  pdf(paste0(Folder_name, "/Top10.pdf"), width = 10, height = 25)
  print(gp)
  dev.off()
  
  
  # output results
  tmp <- as.data.frame(resOrdered)
  idx <- match(rownames(tmp), Gene_Info_Selected[[Genes_Column_Name]])
  if(sum( Gene_Info_Selected[[Genes_Column_Name]][idx] == rownames(tmp)) != nrow(tmp)){
    stop("Not Match")
  }
  
  tmp <- cbind(Gene_Info_Selected[idx,],tmp)
  
  write.csv(tmp, file=paste0(Folder_name, "/rlt.csv"), quote = TRUE, row.names = FALSE)
  
  
  # volcano plot
  dplot <- data.frame(
    fc = res$log2FoldChange,
    pv = -log10(res$pvalue),
    gene = Gene_Info_Selected[[Genes_Column_Name]],
    stringsAsFactors = FALSE
  )
  
  dplot <- dplot[!is.na(dplot$gene),]
  
  gp <- volcanoPlot(dplot, pvCutoff = pvalue_Cutoff, log2FoldCutoff = log2Fold_Cutoff,Variable_Of_Interest)
  
  pdf(paste0(Folder_name, "/volcano.pdf"), width = 5, height = 5)
  print(gp)
  invisible(dev.off())
  
  
  
  sigRlt <- as.data.frame(resOrdered[which(resOrdered$pvalue < pvalue_Cutoff & abs(resOrdered$log2FoldChange) > log2Fold_Cutoff),])
  idx <- match(rownames(sigRlt), Gene_Info_Selected[[Genes_Column_Name]])
  if(sum(Gene_Info_Selected[[Genes_Column_Name]][idx] == rownames(sigRlt)) != nrow(sigRlt)){
    stop("Not Match")
  }
  sigRlt <- cbind(Gene_Info_Selected[idx,], sigRlt)
  sigRlt <- sigRlt[!is.na(sigRlt[[Genes_Column_Name]]),]
  write.csv(sigRlt, file=paste0(Folder_name, "/sigRlt.csv"), quote = TRUE, row.names = FALSE)
  
  # heatmap
  flagNA <- is.na(Gene_Info_Selected$ENTREZID)
  ht <- plotHeat(vsd[!flagNA,], res[!flagNA,], Gene_Info_Selected[!flagNA,], pvCutoff = pvalue_Cutoff, log2FoldCutoff = log2Fold_Cutoff,variable=Variable_Of_Interest)
  pdf(paste0(Folder_name, "/HeatmapSigGenes.pdf"), width = 6, height = 10)
  ht <- draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "right")
  invisible(dev.off())
  
  
  # pathway analysis
  eKEGG <- clusterProfiler::enrichKEGG(
    gene         = sigRlt$ENTREZID,
    organism     = 'mmu',
    keyType = "ncbi-geneid")
  
  
  eGO <- clusterProfiler::enrichGO(
    gene         = sigRlt$ENTREZID,
    OrgDb         = "org.Mm.eg.db",
    keyType = "ENTREZID",
    ont           = "ALL",
    readable      = TRUE
  )
  
  fc <- res$log2FoldChange
  names(fc) <- Gene_Info_Selected$ENTREZID
  
  dim(eKEGG)[1]
  
  if(!(is.null(dim(eKEGG)))){
    if(dim(eKEGG)[1] > 0){ 
      pdf(paste0(Folder_name, "/KEGGDot.pdf"), width = 7, height = 7)
      print(enrichplot::dotplot(eKEGG))
      invisible(dev.off())
      
      pdf(paste0(Folder_name, "/KEGGUpset.pdf"), width = 9, height = 7)
      print(enrichplot::upsetplot(eKEGG))
      invisible(dev.off())
      
      
      eKEGGGS <- DOSE::setReadable(eKEGG, 'org.Mm.eg.db', 'ENTREZID')
      pdf(paste0(Folder_name, "/KEGGNet.pdf"), width = 9, height = 7)
      print(cnetplot(eKEGGGS, color.params = list(foldChange = fc)))
      invisible(dev.off())
      
      write.csv(eKEGGGS, file = paste0(Folder_name, "/eKEGG.csv"), row.names = FALSE, quote = TRUE)
    }
    
  }  
  
  if(dim(eGO)[1] > 0){
    pdf(paste0(Folder_name, "/GODot.pdf"), width = 7, height = 7)
    print(enrichplot::dotplot(eGO))
    invisible(dev.off())
    
    pdf(paste0(Folder_name, "/GOUpset.pdf"), width = 9, height = 7)
    print(enrichplot::upsetplot(eGO))
    invisible(dev.off())
    
    pdf(paste0(Folder_name, "/GONet.pdf"), width = 9, height = 7)
    print(cnetplot(eGO, max.overlaps = 100, color.params = list(foldChange = fc)))
    invisible(dev.off())
    
    write.csv(eGO, file = paste0(Folder_name, "/eGO.csv"), row.names = FALSE, quote = TRUE)
    
    
    
    
    
  }
  
  # Specify the file path and name
  readme_path <- paste0(Folder_name,"/readme.md")
  
  # Create the README file
  file.create(readme_path)
  
  # Write content to the README file
  readme_content <- c(
    paste0("Variable of Interest: ", Variable_Of_Interest),
    "",
    paste0("Groups Selected: ", Groups_Selected),
    "",
    paste0("Covariates: ", Covariates),
    "",
    paste0("pvalue cutoff: ", pvalue_Cutoff),
    "",
    paste0("log2fold cutoff: ", log2Fold_Cutoff),
    "",
    paste0("Extra Filters: ", Extra_Filters)
  )
  
  writeLines(readme_content, readme_path)
  
  save(res, dds, eKEGG, eGO, file = paste0(Folder_name, "/rlt.RData"))
}



#function continuousCompare
#Arguments: Feature_Counts, Sample_Info, Samples_Column_Name, Gene_Info, Genes_Column_Name, Variale_Of_Interest, Covariates, Folder_name, pvalue_Cutoff, Extra_Filters
#Performs the DESeq2 analysis using the counts data(Feature_Counts), Sample Information (Sample_Info), Gene information (Gene_Info) using a continuous variable
#Generates the plots (Sample_Distance, PCA plot, MA plot, Volcano plot, Normalized Counts plot of top 10 genes, Heatmap of significant genes)
#Performs Go and KEGG pathway analysis and generates the result files
continuousCompare <- function(Feature_Counts,Sample_Info, Samples_Column_Name, Gene_Info, Genes_Column_Name, Variable_Of_Interest, Covariates = NULL , Folder_name, pvalue_Cutoff, Extra_Filters = NULL){
  
  #Check if Sample IDs are in same order in Feature counts and Sample Info files
  SampleIDs<-colnames(Feature_Counts)
  SampleIDs_SampleInfo<-Sample_Info[[Samples_Column_Name]]
  if(!identical(SampleIDs,SampleIDs_SampleInfo)){
    Feature_Counts<-Feature_Counts[, match(SampleIDs_SampleInfo,colnames(Feature_Counts))]
  }
  
  #Check if Gene IDs are  in same order in Feature counts and Gene Info files
  genes_geneinfo<-sapply(Gene_Info[[Genes_Column_Name]] , function(x) sub("\\..*", "", x))
  genes_geneinfo<-data.frame(genes_geneinfo)
  genes_geneinfo<-genes_geneinfo[, 1]
  genes<-sapply(rownames(Feature_Counts), function(x) sub("\\..*", "", x))
  genes<-data.frame(genes)
  genes<-genes[, 1]
  if(!identical(genes_geneinfo,genes)){
    Feature_Counts<-Feature_Counts[match(genes_geneinfo,rownames(Feature_Counts)), ]
  }
  
  
  idx <- rowSums(Feature_Counts == 0) < ncol(Feature_Counts)
  Gene_Info_Selected <- Gene_Info[idx,]
  Feature_Counts_Selected <- Feature_Counts[idx,]
  
  if(sum(Gene_Info_Selected[[Genes_Column_Name]] == rownames(Feature_Counts_Selected), na.rm = TRUE) != nrow(Feature_Counts_Selected)){
    stop("Not Match")
  }
  
  
  
  formula<-as.formula(paste("~ ", Variable_Of_Interest))
  
  dds <- DESeqDataSetFromMatrix(
    countData = Feature_Counts_Selected,
    colData = Sample_Info,
    design= formula)
  
  
  dds <- DESeq(dds)
  
  resultsNames(dds) # lists the coefficients
  res <- results(dds, name= Variable_Of_Interest)
  
  resOrdered <- res[order(res$pvalue),]
  
  resOrdered <- resOrdered[complete.cases(resOrdered$padj), ]
  
  resOrdered
  
  class(Sample_Info[[Variable_Of_Interest]])
  
  # Apply variance-stabilizing transformation
  vsd <- vst(dds, blind = TRUE)
  
  # Calculate sample distances
  sampleDists <- dist(t(assay(vsd)))
  sampleDistMatrix <- as.matrix(sampleDists)
  
  
  # Create the heatmap
  heatmap = Heatmap(
    sampleDistMatrix,
    col = circlize::colorRamp2(c(0, max(sampleDistMatrix)), c("white", "red")),
    heatmap_legend_param = list(
      title = "Samp_Dis",
      at = seq(0, max(sampleDistMatrix), length.out = 4),
      labels = prettyNum(seq(0, max(sampleDistMatrix), length.out = 4), big.mark = ",")
    ),
    top_annotation = HeatmapAnnotation(
      variable = vsd[[Variable_Of_Interest]],
      annotation_legend_param = list(variable = list(title = Variable_Of_Interest, legend_direction = "vertical"))
    )
  )
  
  # Create a directory to save the output
  dir.create(Folder_name)
  
  # Save the heatmap as a PDF file
  pdf(paste0(Folder_name, "/sampleDist.pdf"), width = 7, height = 6)
  draw(heatmap, annotation_legend_side = "right")
  invisible(dev.off())
  


  # PCA
  pcaData <- plotPCA(vsd, intgroup= Variable_Of_Interest, returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  gp <- ggplot(pcaData, aes(PC1, PC2, color = get(Variable_Of_Interest),label = name)) +
    geom_point(size=2, show.legend = FALSE) +
    ggrepel::geom_text_repel(key_glyph = "rect", size = 3) + 
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    scale_color_gradient(low = ggsci::pal_nejm("default")(2)[1], high = ggsci::pal_nejm("default")(2)[2]) +
    #ggsci::scale_color_nejm() + 
    #coord_fixed() +
    theme_light() +
    guides(color = guide_legend(byrow = T, keyheight = unit(.1, "cm"))) + 
    guides(color = guide_legend(title = Variable_Of_Interest))
  
  pdf(paste0(Folder_name, "/pcaVSD.pdf"), width = 4, height = 4.2)
  print(gp)
  invisible(dev.off())
  
  
  # MA plot
  pdf(paste0(Folder_name, "/MAPlot.pdf"), width = 5, height = 4)
  plotMA(res, colSig = "red")
  invisible(dev.off())
  
  # top 10 genes
  sigGN <- rownames(resOrdered)[1:10]
  gps <- list()
  i <- 1
  for(s in sigGN){
    idx <- which(Gene_Info_Selected[[Genes_Column_Name]] == s)
    gps[[i]] <- plotCtCon(dds, res, idx,Variable_Of_Interest , Samples_Column_Name)
    i <- i + 1
  }
  
  gp <- ggpubr::ggarrange(plotlist = gps, ncol = 2, nrow = 5)
  pdf(paste0(Folder_name, "/Top10.pdf"), width = 20, height = 25)
  print(gp)
  invisible(dev.off())
  
  res
  
  # output results
  tmp <- as.data.frame(resOrdered)
  idx <- match(rownames(tmp), Gene_Info_Selected[[Genes_Column_Name]])
  if(sum( Gene_Info_Selected[[Genes_Column_Name]][idx] == rownames(tmp)) != nrow(tmp)){
    stop("Not Match")
  }
  
  tmp <- cbind(Gene_Info_Selected[idx,],tmp)
  
  write.csv(tmp, file=paste0(Folder_name, "/rlt.csv"), quote = TRUE, row.names = FALSE)
  
  
  
  sigRlt <- as.data.frame(resOrdered[which(resOrdered$pvalue < pvalue_Cutoff),])
  idx <- match(rownames(sigRlt), Gene_Info_Selected[[Genes_Column_Name]])
  if(sum(Gene_Info_Selected[[Genes_Column_Name]][idx] == rownames(sigRlt)) != nrow(sigRlt)){
    stop("Not Match")
  }
  sigRlt <- cbind(Gene_Info_Selected[idx,], sigRlt)
  sigRlt <- sigRlt[!is.na(sigRlt[[Genes_Column_Name]]),]
  write.csv(sigRlt, file=paste0(Folder_name, "/sigRlt.csv"), quote = TRUE, row.names = FALSE)
  
  
  # heatmap
  flagNA <- is.na(Gene_Info_Selected$ENTREZID)
  ht <- plotHeatCon(vsd[!flagNA,], res[!flagNA,], Gene_Info_Selected[!flagNA,], pvCutoff = pvalue_Cutoff, variable= Variable_Of_Interest, Genes_Column_Name = Genes_Column_Name)
  pdf(paste0(Folder_name, "/HeatmapSigGenes.pdf"), width = 6, height = 25)
  ht <- draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "right")
  invisible(dev.off())
  
  
  # pathway analysis
  eKEGG <- clusterProfiler::enrichKEGG(
    gene         = sigRlt$ENTREZID,
    organism     = 'mmu',
    keyType = "ncbi-geneid")
  
  
  eGO <- clusterProfiler::enrichGO(
    gene         = sigRlt$ENTREZID,
    OrgDb         = "org.Mm.eg.db",
    keyType = "ENTREZID",
    ont           = "ALL",
    readable      = TRUE
  )
  
  pv <- res$pvalue
  names(pv) <- Gene_Info_Selected$ENTREZID
  
  dim(eKEGG)[1]
  
  if(!(is.null(dim(eKEGG)))){
    if(dim(eKEGG)[1] > 0){ 
      pdf(paste0(Folder_name, "/KEGGDot.pdf"), width = 7, height = 7)
      print(enrichplot::dotplot(eKEGG))
      invisible(dev.off())
      
      pdf(paste0(Folder_name, "/KEGGUpset.pdf"), width = 9, height = 7)
      print(enrichplot::upsetplot(eKEGG))
      invisible(dev.off())
      
      
      eKEGGGS <- DOSE::setReadable(eKEGG, 'org.Mm.eg.db', 'ENTREZID')
      pdf(paste0(Folder_name, "/KEGGNet.pdf"), width = 9, height = 7)
      print(cnetplot(eKEGGGS, color.params = list(pval = pv)))
      invisible(dev.off())
      
      write.csv(eKEGGGS, file = paste0(Folder_name, "/eKEGG.csv"), row.names = FALSE, quote = TRUE)
    }
    
  }  
  
  if(dim(eGO)[1] > 0){
    pdf(paste0(Folder_name, "/GODot.pdf"), width = 7, height = 7)
    print(enrichplot::dotplot(eGO))
    invisible(dev.off())
    
    pdf(paste0(Folder_name, "/GOUpset.pdf"), width = 9, height = 7)
    print(enrichplot::upsetplot(eGO))
    invisible(dev.off())
    
    pdf(paste0(Folder_name, "/GONet.pdf"), width = 9, height = 7)
    print(cnetplot(eGO, max.overlaps = 100, color.params = list(pval = pv)))
    invisible(dev.off())
    
    write.csv(eGO, file = paste0(Folder_name, "/eGO.csv"), row.names = FALSE, quote = TRUE)
    
  }
  
  # Specify the file path and name
  readme_path <- paste0(Folder_name,"/readme.md")
  
  # Create the README file
  file.create(readme_path)
  
  # Write content to the README file
  readme_content <- c(
    paste0("Variable of Interest: ", Variable_Of_Interest),
    "",
    paste0("Covariates: ", Covariates),
    "",
    paste0("pvalue cutoff: ", pvalue_Cutoff),
    "",
    paste0("Extra Filters: ", Extra_Filters)
  )
  
  writeLines(readme_content, readme_path)
  
  save(res, dds, eKEGG, eGO, file = paste0(Folder_name, "/rlt.RData"))
  
  
}