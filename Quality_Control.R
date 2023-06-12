
Normalize_Counts <- function(cts,sampleInfo,Variable_Of_Interest,Groups_Selected,Folder_Name){
  
  idx <- sampleInfo[[Variable_Of_Interest]] %in% Groups_Selected
  sampleInfoSel <- sampleInfo[idx,]
  
  sampleInfoSel[[Variable_Of_Interest]] <- factor(sampleInfoSel[[Variable_Of_Interest]], levels = Groups_Selected)
  
  sampleInfoSel[[Variable_Of_Interest]]
  
  flag <- rowSums(cts == 0) < ncol(cts)
  cts <- cts[flag,]
  cts <- cts[, idx]
  
  
  formula<-as.formula(paste("~  ", Variable_Of_Interest))
  
  
  dds <- DESeqDataSetFromMatrix(
    countData = cts,
    colData = sampleInfoSel,
    design = formula)
  
  dds <- DESeq(dds)
  
  ctNorm <- counts(dds, normalized = TRUE, replaced = FALSE)
  
  dir.create(Folder_Name)
  write.csv(ctNorm, file = paste0(Folder_Name, "/normalizedReadCounts.csv"), quote = TRUE)
  
  return(dds)
}

Distance_Clustering <- function(dds,Groups_Selected,Variable_Of_Interest,Folder_Name){
  vsd <- vst(dds, blind=TRUE)
  sampleDist <- dist(t(assay(vsd)))
  
  sampleDistMatrix <- as.matrix(sampleDist)
  
  colTreatment <- c(ggsci::pal_nejm("default")(3)[3], ggsci::pal_nejm("default")(3)[2])
  names(colTreatment) <- Groups_Selected
  
  
  
  ha = HeatmapAnnotation(
    treatment = vsd[[Variable_Of_Interest]],
    col = list(treatment = colTreatment),
    annotation_legend_param = list(
      treatment = list(nrow = 1, direction = "horizontal", title = "")
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
      col =
        ifelse(vsd[[Variable_Of_Interest]] == "Control", ggsci::pal_nejm("default")(3)[3],
               ggsci::pal_nejm("default")(3)[2])),
    
    col = colorRampPalette(rev(brewer.pal(9, "Blues")) )(255),
    column_names_side = "bottom"
  )
  
  pdf(paste0(Folder_Name,"/sampleDist.pdf"), width = 5, height = 6)
  draw(ht,annotation_legend_side = "bottom")
  invisible(dev.off())
}

PCA_Plots <- function(dds, Variable_Of_Interest, Samples_Column_name, sampleInfo , Groups_Selected, Variables_For_PCA,Folder_Name,Color_Choice=NULL, Shape_Choice = NULL){
  
  idx <- sampleInfo[[Variable_Of_Interest]] %in% Groups_Selected
  sampleInfoSel <- sampleInfo[idx,]
  vsd <- vst(dds, blind=TRUE)
  pcaData <- plotPCA(vsd, intgroup=Variable_Of_Interest, returnData=TRUE)
  if(sum(pcaData$name == sampleInfoSel[[Samples_Column_name]]) != nrow(sampleInfoSel)){
    stop("Not Match")
  }
  
  for (i in 1:length(Variables_For_PCA)){
    if(grepl(":",Variables_For_PCA[i])){
      if(strsplit(Variables_For_PCA[i], ":")[[1]][2] == "numeric"){
        pcaData[[strsplit(Variables_For_PCA[i], ":")[[1]][1]]] <- as.numeric(sampleInfoSel[[strsplit(Variables_For_PCA[i], ":")[[1]][1]]])
      }
    }else{
      pcaData[[Variables_For_PCA[i]]] <- sampleInfoSel[[Variables_For_PCA[i]]]
    }
  }
  
  
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  gp <- ggplot(pcaData, aes(PC1, PC2, color=group, label = name, shape = type)) +
    geom_point(size=2, show.legend = TRUE) +
    ggrepel::geom_text_repel(key_glyph = "rect", size = 3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    (if (is.null(Color_Choice)){scale_color_manual(values = c("red","blue"))}
     else{scale_color_manual(values = Color_Choice)}) +
    (if (is.null(Shape_Choice)){}
     else{scale_shape_manual(values = Shape_Choice)}) +
    #ggsci::scale_color_nejm() +
    #coord_fixed() +
    theme_light() +
    guides(color = guide_legend(byrow = T, keyheight = unit(.1, "cm"))) +
    theme(legend.title = element_blank(), legend.box="vertical", legend.position = "bottom")
  
  pdf(paste0(Folder_Name,"/pcaVSD.pdf"), width = 6, height = 6.5)
  print(gp)
  invisible(dev.off())
  
  
  gp <- ggplot(pcaData, aes(PC1, PC2, color=RIN, label = name, shape = group)) +
    geom_point(size=2, show.legend = TRUE) +
    ggrepel::geom_text_repel(size = 3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    (if(is.null(Color_Choice)){scale_color_gradient(low = "red", high = "blue") }
     else{scale_color_gradient(low = Color_Choice[1], high = Color_Choice[2])})+
    (if(is.null(Shape_Choice)){}
     else{scale_shape_manual(values = Shape_Choice)})+
    #scale_color_manual(values = ggsci::pal_nejm("default")(3)[c(3,2)]) +
    #ggsci::scale_color_nejm() +
    #coord_fixed() +
    labs(shape = "") +
    theme_light() +
    theme(legend.box="vertical", legend.position = "bottom")
  
  pdf(paste0(Folder_Name,"/pcaVSDRIN.pdf"), width = 6, height = 6.5)
  print(gp)
  invisible(dev.off())
  
}

Y_Reads <- function(cts,geneInfo,sampleInfo,Folder_Name, Chromosome_CN, gender_column_name, Samples_column_name){
  
  numTotalReads <- colSums(cts)
  flagChrY <- stringr::str_detect(geneInfo[[Chromosome_CN]], "chrY")
  ctsY <- cts[flagChrY,]
  dplot <- data.frame(
    ct = colSums(ctsY),
    ctPercent =  colSums(ctsY)/numTotalReads,
    sex = sampleInfo[[gender_column_name]],
    id = sampleInfo[[Samples_column_name]],
    stringsAsFactors = FALSE
  )
  
  
  pos <- position_jitter(width = 0.2, height = 0, seed = 2)
  
  gp <- ggplot(dplot, aes(x=sex, y=ct)) +
    geom_boxplot(outlier.shape = NA)+
    #geom_boxplot(outlier.shape = NA) +
    geom_jitter(position = pos) +
    ggrepel::geom_text_repel(aes(label = id), position = pos) +
    labs(x="", y="Total Reads on ChrY") +
    theme_light()
  
  pdf(paste0(Folder_Name,"/chrYReads.pdf"), width = 5, height = 4)
  print(gp)
  invisible(dev.off())
  
  
  gp <- ggplot(dplot, aes(x=sex, y=ctPercent)) +
    geom_boxplot(outlier.shape = NA)+
    #geom_boxplot(outlier.shape = NA) +
    geom_jitter(position = pos) +
    ggrepel::geom_text_repel(aes(label = id), position = pos) +
    labs(x="", y="Chr Y Reads Proportion") +
    theme_light()
  
  pdf(paste0(Folder_Name,"/chrYReadsProportion.pdf"), width = 5, height = 4)
  print(gp)
  invisible(dev.off())
}

X_Reads <- function(cts,geneInfo,sampleInfo,Folder_Name, Chromosome_CN, gender_column_name, Samples_column_name){
  
  numTotalReads <- colSums(cts)
  flagChrX <- stringr::str_detect(geneInfo[[Chromosome_CN]], "chrX")
  ctsX <- cts[flagChrX,]
  dplot <- data.frame(
    ct = colSums(ctsX),
    ctPercent =  colSums(ctsX)/numTotalReads,
    sex = sampleInfo[[gender_column_name]],
    id = sampleInfo[[Samples_column_name]],
    stringsAsFactors = FALSE
  )
  
  pos <- position_jitter(width = 0.2, height = 0, seed = 2)
  
  gp <- ggplot(dplot, aes(x=sex, y=ct)) +
    geom_boxplot(outlier.shape = NA)+
    #geom_boxplot(outlier.shape = NA) +
    geom_jitter(position = pos) +
    ggrepel::geom_text_repel(aes(label = id), position = pos) +
    labs(x="", y="Total Reads on ChrX") +
    theme_light()
  
  pdf(paste0(Folder_Name,"/chrXReads.pdf"), width = 5, height = 4)
  print(gp)
  invisible(dev.off())
  
  
  gp <- ggplot(dplot, aes(x=sex, y=ctPercent)) +
    geom_boxplot(outlier.shape = NA)+
    #geom_boxplot(outlier.shape = NA) +
    geom_jitter(position = pos) +
    ggrepel::geom_text_repel(aes(label = id), position = pos) +
    labs(x="", y="Chr X Reads Proportion") +
    theme_light()
  
  pdf(paste0(Folder_Name,"/chrXReadsProportion.pdf"), width = 5, height = 4)
  print(gp)
  invisible(dev.off())
  
  
}

XIST_Counts <- function(cts,geneInfo,genes_column_name,sampleInfo,Folder_Name,gender_column_name,Samples_column_name){
  
  numTotalReads <- colSums(cts)
  dplot <- data.frame(
    Xist = cts[which(toupper(geneInfo[[genes_column_name]]) == toupper("Xist")),],
    XistNorm =  cts[which(toupper(geneInfo[[genes_column_name]]) == toupper("Xist")),]/numTotalReads,
    sex = sampleInfo[[gender_column_name]],
    id = sampleInfo[[Samples_column_name]]
  )
  
  pos <- position_jitter(width = 0.2, height = 0, seed = 2)
  
  gp <- ggplot(dplot, aes(x=sex, y=Xist)) +
    geom_boxplot(outlier.shape = NA)+
    geom_jitter(position = pos) +
    ggrepel::geom_text_repel(aes(label = id), position = pos) +
    #scale_y_log10() +
    labs(x="", y="Num Reads Xist") +
    theme_light()
  
  pdf(paste0(Folder_Name,"/Xist.pdf"), width = 5, height = 4)
  print(gp)
  invisible(dev.off())
  
  gp <- ggplot(dplot, aes(x=sex, y=XistNorm)) +
    geom_boxplot(outlier.shape = NA)+
    #geom_boxplot(outlier.shape = NA) +
    geom_jitter(position = pos) +
    ggrepel::geom_text_repel(aes(label = id), position = pos) +
    #scale_y_log10() +
    labs(x="", y="Normalized Xist") +
    theme_light()
  
  pdf(paste0(Folder_Name,"/XistNorm.pdf"), width = 5, height = 4)
  print(gp)
  invisible(dev.off())
  
  
}