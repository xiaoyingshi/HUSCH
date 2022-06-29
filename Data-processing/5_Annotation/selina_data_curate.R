clearFolder <- function(file, organism, idType){
  setwd(paste(strsplit(file,'/')[[1]][1:9], collapse='/'))
  rds <- readRDS(file)
  if(length(rds)==1){
    rds=rds
  }else{
    rds=rds$RNA
  }   
  # cell label
  allct <- sort(unique(rds$level3))
  allct <- allct[!grepl('Unknown',allct)]
  allct <- sapply(allct, function(x) gsub('_lv3','',x), USE.NAMES = F)
  
  # marker
  if(file.exists(gsub('.rds','_Marker.txt',file))){
    allsignatures <- read.table(gsub('.rds','_Marker.txt',file), header=T, sep = '\t')
    allcelltypes <- as.character(unique(allsignatures[, 1]))
    allcelltypes <- allcelltypes[allcelltypes %in% allct]
    allsignatures_list <- sapply(1:length(allcelltypes), function(x) {
      return(toupper(as.character(allsignatures[which(allsignatures[, 1] == allcelltypes[x]), 2])))
    })
    names(allsignatures_list) <- allcelltypes
    
    for(ct in unique(allct[allct %in% names(allsignatures_list)])){
      print(ct)
      res_marker=unique(allsignatures_list[ct][[1]])
      ggsave(paste0(rds@project.name,'_',ct,'_Feature.png'), FeaturePlot(rds, features = res_marker, slot = 'data'), height = 5, width = 7)    
    }
  }
  
  # clustering plot
  system('mkdir processed')
  if('MAESTRO'  %in% list.files('.')){
    system('mv MAESTRO processed')
  }
  ggsave(paste0('processed/',rds@project.name,'_cluster.png'), DimPlot(rds, group.by = "seurat_clusters"), height = 5, width = 7)
  ggsave(paste0('processed/',rds@project.name,'_celltype.png'), DimPlot(rds, group.by = "level3"), height = 5, width = 7)
  
  allctDF <- data.frame(allct, stringsAsFactors = F)
  colnames(allctDF) <- 'Celltype'
  write.table(allctDF, paste0('processed/',rds@project.name,'_Label.txt'), sep='\t', quote = F, row.names = F)
  
  # Feature Plot
  ct_genes <- FindMarkers(rds,rds$level3)
  write.table(ct_genes,file = paste0(rds@project.name,'_DiffGenes.tsv'),row.names = F,col.names = T,sep = '\t',quote = F)
  
  matCount <- rds@assays$RNA@counts
  write.table(matCount,file = paste0(rds@project.name,'_count.tsv'),col.names = T,sep = '\t',quote = F)
  
  matTPM <- RNACountToTPM(matCount, idType = "Symbol", organism = organism)
  write.table(matTPM,file = paste0(rds@project.name,'_TPM.tsv'),col.names = T,sep = '\t',quote = F)
  
  matTPMavg <- apply(matTPM, 1, function(t) tapply(t, rds$seurat_clusters, mean))
  write.table(matTPMavg,file = paste0(rds@project.name,'_avgTPM.tsv'),col.names = T,sep = '\t',quote = F)
}

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}

FindMarkers <- function(object, cluster, features = NULL, min.pct = 0.1, logfc.threshold = 0.25,
                        only.pos = FALSE, return.thresh = 1e-2,
                        slot = "data")
{
  matrix = GetAssayData(object, slot = slot)
  features = rownames(matrix)
  y <- cluster
  y = factor(y)
  test.res = wilcoxauc(matrix, y)
  
  # Calculate logFC
  if (slot != "scale.data"){
    if (slot == "data"){
      X = expm1(matrix)
    }
    group_sums = sumGroups(X, y, 1)
    group_means <- sweep(group_sums, 1, as.numeric(table(y)), "/") %>% t()
    cs <- colSums(group_sums)
    gs <- as.numeric(table(y))
    lfc <- Reduce(cbind, lapply(seq_len(length(levels(y))), function(g) {
      log(group_means[, g]+1) - log(((cs - group_sums[g, ]) / (length(y) - gs[g]))+1)
    }))
  }else{
    group_sums = sumGroups(X, y, 1)
    group_means <- sweep(group_sums, 1, as.numeric(table(y)), "/") %>% t()
    cs <- colSums(group_sums)
    gs <- as.numeric(table(y))
    lfc <- Reduce(cbind, lapply(seq_len(length(levels(y))), function(g) {
      group_means[, g] - ((cs - group_sums[g, ]) / (length(y) - gs[g]))
    }))
  }
  
  test.res$avg_logFC = as.vector(lfc)
  res = test.res[,c("pval", "avg_logFC", "pct_in", "pct_out", "padj", "group", "feature")]
  res[,c("pct_in", "pct_out")] = round(res[,c("pct_in", "pct_out")]/100, digits = 3)
  colnames(res) = c("p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene")
  res <- res %>% dplyr::filter(.data$p_val < return.thresh & 
                                 abs(.data$avg_logFC) > logfc.threshold &
                                 (.data$pct.1 > min.pct |
                                    .data$pct.2 > min.pct) &
                                 .data$gene %in% features)
  if(only.pos){
    res <- res %>% dplyr::filter(.data$avg_logFC > 0)
  }
  res <- res %>% dplyr::arrange(.data$cluster, .data$p_val, desc(.data$avg_logFC))
  return(res)
}

calMean <- function(data) {
  celltype <- unique(data$level3)
  expr <- data@assays$RNA@counts
  if (length(celltype)==1){
    avg = as.data.frame(Matrix::rowMeans(expr))
  } else {
    avg <- lapply(celltype, function(ct) {
      expr <- expr[, which(data$level3 == ct)]
      if (is.vector(expr)) {
        avg_expr <- expr
      } else {
        avg_expr <- Matrix::rowMeans(expr)
      }
      return(avg_expr)
    })
    avg <- avg %>% reduce(cbind)
  }
  colnames(avg) <- paste(celltype, data@project.name,sep=';')
  return(avg)
}


CCA <- function(object, batch) {
  cells <- ncol(object)
  cluster.res <- 1
  if (cells <= 100000) {
    npc <- 50
  } else {
    npc <- 100    
  }
  object <- NormalizeData(object, verbose = FALSE)
  object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  object <- ScaleData(object, verbose = FALSE)
  object <- RunPCA(object, npcs = npc, verbose = FALSE, features = VariableFeatures(object))
  pc.contribution <- object@reductions$pca@stdev / sum(object@reductions$pca@stdev) * 100
  pc.contribution.cum <- cumsum(pc.contribution)
  pc.first <- which(pc.contribution.cum > 75)[1]
  dims.use <- 1:pc.first
  data.list <- SplitObject(object, split.by = batch)
  for (i in 1:length(data.list)) {
    data.list[[i]] <- NormalizeData(data.list[[i]], verbose = FALSE)
    data.list[[i]] <- FindVariableFeatures(data.list[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  }
  anchors <- FindIntegrationAnchors(object.list = data.list, dims = dims.use, anchor.features = 2000)
  RNA.integrated <- IntegrateData(anchorset = anchors, dims = dims.use, k.weight=10)
  DefaultAssay(RNA.integrated) <- "integrated"
  
  RNA.integrated <- ScaleData(RNA.integrated, verbose = FALSE)
  RNA.integrated <- RunPCA(RNA.integrated, npcs = npc, verbose = FALSE, features = VariableFeatures(RNA.integrated))
  
  pc.contribution <- RNA.integrated@reductions$pca@stdev / sum(RNA.integrated@reductions$pca@stdev) * 100
  pc.contribution.cum <- cumsum(pc.contribution)
  pc.first <- which(pc.contribution.cum > 75)[1]
  dims.use <- 1:pc.first
  RNA.integrated <- RunUMAP(object = RNA.integrated, reduction = "pca", dims = dims.use)
  
  RNA.integrated <- FindNeighbors(object = RNA.integrated, reduction = "pca", dims = dims.use)
  RNA.integrated <- FindClusters(object = RNA.integrated, resolution = cluster.res)
  
  return(RNA.integrated)
}

count2seurat <- function(counts, npc) {
  data <- CreateSeuratObject(counts = counts)
  data <- NormalizeData(data, verbose = FALSE)
  data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  data <- ScaleData(data, verbose = FALSE)
  data <- RunPCA(data, npcs = npc, verbose = FALSE, features = VariableFeatures(data))
}


lineage_extract <- function(lineage_cells, level) {
  seurat_list = list()
  count = 0
  for (i in 1:length(data_list)) {
    if (sum(lineage_cells %in% data_list[[i]]@meta.data[, level]) != 0) {
      count = count + 1
      seurat_list[[count]] = data_list[[i]][, which(data_list[[i]]@meta.data[, level] %in% lineage_cells)]
      seurat_list[[count]] = SeuratObject::RenameCells(seurat_list[[count]], new.names = paste0(data_list[[i]]@project.name, ';', colnames(seurat_list[[count]]@assays$RNA@data)))
      seurat_list[[count]]$project = data_list[[i]]@project.name
    }
  }
  seurat_obj = seurat_list %>% reduce(merge)
  return(seurat_obj)
}

harmony <- function(object, batch) {
  cells <- ncol(object)
  if (cells <= 5000) {
    npc <- 50;
    cluster.res <- 0.6
  } else if (cells <= 100000) {
    npc <- 50;
    cluster.res <- 1
  } else {
    npc <- 100;
    cluster.res <- 1
  }
  object <- NormalizeData(object, verbose = FALSE)
  object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  object <- ScaleData(object, verbose = FALSE)
  object <- RunPCA(object, npcs = npc, verbose = FALSE, features = VariableFeatures(object))
  pc.contribution <- object@reductions$pca@stdev / sum(object@reductions$pca@stdev) * 100
  pc.contribution.cum <- cumsum(pc.contribution)
  pc.first <- which(pc.contribution.cum > 75)[1]
  dims.use <- 1:pc.first
  object <- RunHarmony(object, batch, verbose = FALSE) %>%
    RunUMAP(reduction = "harmony", dims = dims.use, verbose = FALSE) %>%
    FindNeighbors(reduction = "harmony", dims = dims.use, verbose = FALSE) %>%
    FindClusters(resolution = cluster.res, verbose = FALSE)
  return(object)
}

#level assign
subtype_assign <- function(subtype, level, Parent) {
  Children = unique(subtype)
  for (i in 1:length(data_list)) {
    for (j in 1:length(Children)) {
      cell_parent = colnames(data_list[[i]])[which(data_list[[i]]@meta.data[,level] == Parent)]
      Child = Children[j]
      cell_children = names(subtype[subtype == Child])
      names(cell_children) = unlist(lapply(cell_children, function(x) {
        x = unlist(str_split(x, ';'))[1]
      }))
      cell_children = unlist(lapply(cell_children, function(x) {
        x = unlist(str_split(x, ';'))[2]
      }))
      cell_children = cell_children[names(cell_children) == data_list[[i]]@project.name]
      cell_children = cell_parent[cell_parent %in% cell_children]
      data_list[[i]]@meta.data[cell_children, level] = Child
      data_list[[i]]@meta.data[cell_children, 'temp.assign'] = Child
    }
  }
  return(data_list)
}

maestro_subtype <- function(seurat_obj, signatures, Parent, level) {
  celltypes <- as.character(unique(signatures[, 1]))
  signature_list <- sapply(1:length(celltypes), function(x) {
    return(toupper(as.character(signatures[which(signatures[, 1] == celltypes[x]), 2])))
  })
  names(signature_list) <- celltypes
  if (Parent == 'CD4T') {
    signature_list = signature_list[c(2, 5, 7, 10:12)]
  } else if (Parent == 'CD8T') {
    signature_list = signature_list[c(1, 3, 4, 6, 8)]
  }
  if (length(unique(seurat_obj$seurat_clusters))>=2){
    genes <- FindAllMarkersMAESTRO(seurat_obj, test.use = 'presto', min.pct = 0.1, logfc.threshold = 0.25, only.pos = FALSE, verbose = FALSE, return.thresh = 1e-2, min.cells.feature = 3, min.cells.group = 3, slot = "data")
    cluster_celltype_score = sapply(as.integer(unique(seurat_obj$seurat_clusters)) - 1, function(x) {
      idx = genes$cluster == x
      avglogFC = genes$avg_logFC[idx]
      names(avglogFC) = toupper(genes$gene[idx])
      score_cluster = sapply(signature_list, function(y) {
        score = sum(avglogFC[y], na.rm = TRUE) / log2(length(y))
        return(score)
      })
    })
    colnames(cluster_celltype_score) = as.character(as.integer(unique(seurat_obj$seurat_clusters)) - 1)
    
    cellscore_max_celltype = apply(cluster_celltype_score, 2, function(x) {
      if (max(x) < 0) {
        return(paste0(Parent, '_Unknown'))
      } else {
        return(rownames(cluster_celltype_score)[which.max(x)])
      }
    })
    clusters = as.character(as.integer(seurat_obj$seurat_clusters) - 1)
    current.cluster.ids = unique(clusters)
    new.cluster.ids = cellscore_max_celltype
    subtype = plyr::mapvalues(x = clusters, from = current.cluster.ids, to = new.cluster.ids)
    if (Parent == 'T') {
      current.cluster.ids = sort(celltypes)
      new.cluster.ids = c("CD4T", "CD8T", "CD4T", "CD4T", "CD4T", "CD4T")
      subtype = plyr::mapvalues(x = subtype, from = current.cluster.ids, to = new.cluster.ids)
    }
    names(subtype) = colnames(seurat_obj)
    
    seurat_obj@meta.data[, level] = subtype
    print(DimPlot(seurat_obj, group.by = 'project') + theme(legend.text = element_text(size = 14)))
    options(repr.plot.height = 6, repr.plot.width = 10)
    
    print(DimPlot(seurat_obj, group.by = 'curated.assign', label = TRUE, repel = TRUE, label.size = 4.5) + theme(legend.text = element_text(size = 14)))
    options(repr.plot.height = 6, repr.plot.width = 8)
    
    print(DimPlot(seurat_obj, group.by = 'seurat_clusters', label = TRUE, repel = TRUE, label.size = 4.5) + theme(legend.text = element_text(size = 14)))
    options(repr.plot.height = 6, repr.plot.width = 8)
    
    print(DimPlot(seurat_obj, group.by = level, label = TRUE, repel = TRUE, label.size = 4.5) + theme(legend.text = element_text(size = 14)))
    options(repr.plot.height = 6, repr.plot.width = 9)
    return(subtype_assign(subtype, level, Parent))
  } else {
    seurat_obj@meta.data[, level] = Parent
    print(DimPlot(seurat_obj, group.by = 'project') + theme(legend.text = element_text(size = 14)))
    options(repr.plot.height = 6, repr.plot.width = 10)
    
    print(DimPlot(seurat_obj, group.by = 'curated.assign', label = TRUE, repel = TRUE, label.size = 4.5) + theme(legend.text = element_text(size = 14)))
    options(repr.plot.height = 6, repr.plot.width = 8)
    
    print(DimPlot(seurat_obj, group.by = 'seurat_clusters', label = TRUE, repel = TRUE, label.size = 4.5) + theme(legend.text = element_text(size = 14)))
    options(repr.plot.height = 6, repr.plot.width = 8)
    
    print(DimPlot(seurat_obj, group.by = level, label = TRUE, repel = TRUE, label.size = 4.5) + theme(legend.text = element_text(size = 14)))
    options(repr.plot.height = 6, repr.plot.width = 9)
    return(data_list)
  } 
}



#lv3
l3_assign <- function(seurat_obj, children_clusters, Parent, plot_size = 'small') {
  new.cluster.ids = c()
  for (i in 1:length(children_clusters)) {
    children_clusters[[i]] = children_clusters[[i]] + 1
    if (names(children_clusters)[i] == Parent) {
      new.cluster.ids[children_clusters[[i]]] = paste0(Parent, '_lv3')
    } else {
      new.cluster.ids[children_clusters[[i]]] = names(children_clusters)[i]
    }
  }
  clusters = as.character(as.integer(seurat_obj$seurat_clusters) - 1)
  current.cluster.ids = unique(as.character(sort(as.integer(seurat_obj$seurat_clusters) - 1)))
  subtype = plyr::mapvalues(x = clusters, from = current.cluster.ids, to = new.cluster.ids)
  names(subtype) = colnames(seurat_obj)
  seurat_obj$level3 = subtype
  
  print(DimPlot(seurat_obj, group.by = 'project') + theme(legend.text = element_text(size = 14)))
  options(repr.plot.height = 6, repr.plot.width = 10)
  
  print(DimPlot(seurat_obj, group.by = 'seurat_clusters', label = TRUE, repel = TRUE, label.size = 4.5) + theme(legend.text = element_text(size = 14)))
  options(repr.plot.height = 6, repr.plot.width = 8)
  
  if (plot_size == 'large') {
    print(DimPlot(seurat_obj, group.by = 'level3', label = TRUE, repel = TRUE, label.size = 4.5) + theme(legend.text = element_text(size = 14)))
    options(repr.plot.height = 6, repr.plot.width = 15)
  } else {
    print(DimPlot(seurat_obj, group.by = 'level3', label = TRUE, repel = TRUE, label.size = 4.5) + theme(legend.text = element_text(size = 14)))
    options(repr.plot.height = 6, repr.plot.width = 9)
  }
  
  names(subtype) = colnames(seurat_obj)
  return(subtype_assign(subtype, 'level3', Parent))
}

lymphocyte_subtype <- function(seurat_obj, Parent, level) {
  human.immune.CIBERSORT <- readRDS('/fs/home/shixiaoying/Project/Notebook/SELINA_Data_Curated/2_Annotation/Signature/human.immune.CIBERSORT.rds')
  lymphocyte_celltype=c("MemoryBcells","NaiveBcells","ActNK","RestNK","ActMemCD4Tcells","CD8Tcells","NaiveCD4Tcells","RestMemCD4Tcells","Tfh","Treg")
  lymphocyte_list <- sapply(1:length(lymphocyte_celltype),function(x){
    return(toupper(as.character(human.immune.CIBERSORT[which(human.immune.CIBERSORT[,1]==lymphocyte_celltype[x]),2])))
  })
  names(lymphocyte_list) <- lymphocyte_celltype
  
  if (length(unique(seurat_obj$seurat_clusters))>=2){
    genes <- FindAllMarkersMAESTRO(seurat_obj, test.use = 'presto', min.pct = 0.1, logfc.threshold = 0.25, only.pos = FALSE, verbose = FALSE, return.thresh = 1e-2, min.cells.feature = 3, min.cells.group = 3, slot = "data")
    cluster_celltype_score = sapply(as.integer(unique(seurat_obj$seurat_clusters)) - 1, function(x) {
      idx = genes$cluster == x
      avglogFC = genes$avg_logFC[idx]
      names(avglogFC) = toupper(genes$gene[idx])
      score_cluster = sapply(lymphocyte_list, function(y) {
        score = sum(avglogFC[y], na.rm = TRUE) / log2(length(y))
        return(score)
      })
    })
    colnames(cluster_celltype_score) = as.character(as.integer(unique(seurat_obj$seurat_clusters)) - 1)
    
    cellscore_max_celltype = apply(cluster_celltype_score, 2, function(x) {
      if (max(x) < 0) {
        return(paste0(Parent, '_Unknown'))
      } else {
        return(rownames(cluster_celltype_score)[which.max(x)])
      }
    })
    clusters = as.character(as.integer(seurat_obj$seurat_clusters) - 1)
    current.cluster.ids = unique(clusters)
    new.cluster.ids = cellscore_max_celltype
    subtype = plyr::mapvalues(x = clusters, from = current.cluster.ids, to = new.cluster.ids)
    
    current.cluster.ids = sort(lymphocyte_celltype)
    new.cluster.ids = c("CD4T","NK","CD8T","B","B","CD4T","CD4T","NK","Tfh","Treg")
    
    subtype = plyr::mapvalues(x = subtype, from = current.cluster.ids, to = new.cluster.ids)
    
    names(subtype) = colnames(seurat_obj)
    
    seurat_obj@meta.data[, level] = subtype
    print(DimPlot(seurat_obj, group.by = 'project') + theme(legend.text = element_text(size = 14)))
    options(repr.plot.height = 6, repr.plot.width = 10)
    
    print(DimPlot(seurat_obj, group.by = 'curated.assign', label = TRUE, repel = TRUE, label.size = 4.5) + theme(legend.text = element_text(size = 14)))
    options(repr.plot.height = 6, repr.plot.width = 8)
    
    print(DimPlot(seurat_obj, group.by = 'seurat_clusters', label = TRUE, repel = TRUE, label.size = 4.5) + theme(legend.text = element_text(size = 14)))
    options(repr.plot.height = 6, repr.plot.width = 8)
    
    print(DimPlot(seurat_obj, group.by = level, label = TRUE, repel = TRUE, label.size = 4.5) + theme(legend.text = element_text(size = 14)))
    options(repr.plot.height = 6, repr.plot.width = 9)
    return(subtype_assign(subtype, level, Parent))
  } else {
    seurat_obj@meta.data[, level] = Parent
    print(DimPlot(seurat_obj, group.by = 'project') + theme(legend.text = element_text(size = 14)))
    options(repr.plot.height = 6, repr.plot.width = 10)
    
    print(DimPlot(seurat_obj, group.by = 'curated.assign', label = TRUE, repel = TRUE, label.size = 4.5) + theme(legend.text = element_text(size = 14)))
    options(repr.plot.height = 6, repr.plot.width = 8)
    
    print(DimPlot(seurat_obj, group.by = 'seurat_clusters', label = TRUE, repel = TRUE, label.size = 4.5) + theme(legend.text = element_text(size = 14)))
    options(repr.plot.height = 6, repr.plot.width = 8)
    
    print(DimPlot(seurat_obj, group.by = level, label = TRUE, repel = TRUE, label.size = 4.5) + theme(legend.text = element_text(size = 14)))
    options(repr.plot.height = 6, repr.plot.width = 9)
    return(data_list)
  } 
}

stromal_list <- list()
stromal_list['Fibroblast'] = list(c('COL1A2','COL1A1'))
stromal_list['Endothelial'] = list(c('PECAM1', 'VWF'))

stromal_subtype <- function(seurat_obj, stromal_list, Parent, level) {
  if (length(unique(seurat_obj$seurat_clusters))>=2){
    genes <- FindAllMarkersMAESTRO(seurat_obj, test.use = 'presto', min.pct = 0.1, logfc.threshold = 0.25, only.pos = FALSE, verbose = FALSE, return.thresh = 1e-2, min.cells.feature = 3, min.cells.group = 3, slot = "data")
    cluster_celltype_score = sapply(as.integer(unique(seurat_obj$seurat_clusters)) - 1, function(x) {
      idx = genes$cluster == x
      avglogFC = genes$avg_logFC[idx]
      names(avglogFC) = toupper(genes$gene[idx])
      score_cluster = sapply(stromal_list, function(y) {
        score = sum(avglogFC[y], na.rm = TRUE) / log2(length(y))
        return(score)
      })
    })
    colnames(cluster_celltype_score) = as.character(as.integer(unique(seurat_obj$seurat_clusters)) - 1)
    
    cellscore_max_celltype = apply(cluster_celltype_score, 2, function(x) {
      if (max(x) < 0) {
        return(paste0(Parent, '_Unknown'))
      } else {
        return(rownames(cluster_celltype_score)[which.max(x)])
      }
    })
    clusters = as.character(as.integer(seurat_obj$seurat_clusters) - 1)
    current.cluster.ids = unique(clusters)
    new.cluster.ids = cellscore_max_celltype
    subtype = plyr::mapvalues(x = clusters, from = current.cluster.ids, to = new.cluster.ids)
    
    names(subtype) = colnames(seurat_obj)
    
    seurat_obj@meta.data[, level] = subtype
    print(DimPlot(seurat_obj, group.by = 'project') + theme(legend.text = element_text(size = 14)))
    options(repr.plot.height = 6, repr.plot.width = 10)
    
    print(DimPlot(seurat_obj, group.by = 'curated.assign', label = TRUE, repel = TRUE, label.size = 4.5) + theme(legend.text = element_text(size = 14)))
    options(repr.plot.height = 6, repr.plot.width = 8)
    
    print(DimPlot(seurat_obj, group.by = 'seurat_clusters', label = TRUE, repel = TRUE, label.size = 4.5) + theme(legend.text = element_text(size = 14)))
    options(repr.plot.height = 6, repr.plot.width = 8)
    
    print(DimPlot(seurat_obj, group.by = level, label = TRUE, repel = TRUE, label.size = 4.5) + theme(legend.text = element_text(size = 14)))
    options(repr.plot.height = 6, repr.plot.width = 9)
    return(subtype_assign(subtype, level, Parent))
  } else {
    seurat_obj@meta.data[, level] = Parent
    print(DimPlot(seurat_obj, group.by = 'project') + theme(legend.text = element_text(size = 14)))
    options(repr.plot.height = 6, repr.plot.width = 10)
    
    print(DimPlot(seurat_obj, group.by = 'curated.assign', label = TRUE, repel = TRUE, label.size = 4.5) + theme(legend.text = element_text(size = 14)))
    options(repr.plot.height = 6, repr.plot.width = 8)
    
    print(DimPlot(seurat_obj, group.by = 'seurat_clusters', label = TRUE, repel = TRUE, label.size = 4.5) + theme(legend.text = element_text(size = 14)))
    options(repr.plot.height = 6, repr.plot.width = 8)
    
    print(DimPlot(seurat_obj, group.by = level, label = TRUE, repel = TRUE, label.size = 4.5) + theme(legend.text = element_text(size = 14)))
    options(repr.plot.height = 6, repr.plot.width = 9)
    return(data_list)
  } 
}


## annotation

MarkerAnnot <- function(RNA, genes, clusterid, signatures = signature_list, min.score = 0){
  celltypes <- names(signatures)
  cluster_celltype_score =  sapply(clusterid, function(x){
    idx = genes$cluster==x
    avglogFC = genes$avg_logFC[idx]
    names(avglogFC) = toupper(genes$gene[idx])
    score_cluster = sapply(signature_list, function(y){
      score = sum(avglogFC[y], na.rm = TRUE) / log2(length(y))
      return(score)
    })
  })
  colnames(cluster_celltype_score) = clusterid
  cluster_celltype_score[is.na(cluster_celltype_score)] <- 0
  cellscore_max = apply(cluster_celltype_score, 2, max, na.rm = TRUE)
  cellscore_max_celltype = apply(cluster_celltype_score, 2, function(x){
    if (max(x) < min.score){
      return("Unassigned")
    }else{
      return(rownames(cluster_celltype_score)[which.max(x)])
    }
  })
  return(data.frame(cellscore_max_celltype, stringsAsFactors = F))  
}



##function of run batch effect removal since nn2 problem by adjust k.filter, k.weight
RNABatchCorrect_new <- function(RNA, batch, nfeatures = 2000, dims.use = 1:15, cluster.res = 0.6, only.pos = FALSE, genes.test.use = "presto",
                                genes.cutoff = 1E-5, genes.pct = 0.1, genes.logfc = 0.25,k.filter, k.weight,
                                runpca.agrs = list(), findneighbors.args = list(),
                                findclusters.args = list(), ...){
  RNA@meta.data$batch <- batch
  data.list <- SplitObject(RNA, split.by = "batch")
  for(i in 1:length(data.list)){
    data.list[[i]] <- NormalizeData(data.list[[i]], verbose = FALSE)
    data.list[[i]] <- FindVariableFeatures(data.list[[i]], selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
  }
  
  anchors <- FindIntegrationAnchors(object.list = data.list, dims = dims.use, anchor.features = nfeatures, k.filter = k.filter)
  RNA.integrated <- IntegrateData(anchorset = anchors, dims = dims.use,k.weight = k.weight)
  RNA.integrated@project.name <- paste0(RNA@project.name,'_CCA')
  DefaultAssay(RNA.integrated) <- "integrated"
  
  RNA.integrated <- ScaleData(RNA.integrated, verbose = FALSE)
  RNA.integrated <- fastDoCall("RunPCA", c(object = RNA.integrated, runpca.agrs))
  p = ElbowPlot(object = RNA.integrated, ndims = RNA.integrated@commands$RunPCA.integrated@params$npcs)
  ggsave(file.path(paste0(RNA.integrated@project.name, "_PCElbowPlot.png")), p,  width=10, height=4)
  
  RNA.integrated <- RunUMAP(object = RNA.integrated, reduction = "pca", dims = dims.use, ...)
  RNA.integrated <- fastDoCall("FindNeighbors", c(object = RNA.integrated, reduction = "pca", dims = dims.use, findneighbors.args))
  RNA.integrated <- fastDoCall("FindClusters", c(object = RNA.integrated, resolution = cluster.res, findclusters.args))
  
  p = DimPlot(object = RNA.integrated, label = TRUE, pt.size = 0.2)
  ggsave(file.path(paste0(RNA.integrated@project.name, "_cluster.png")), p, width=5, height=4)
  p = DimPlot(object = RNA.integrated, group = 'batch', label = TRUE, pt.size = 0.2)
  ggsave(file.path(paste0(RNA.integrated@project.name, "_batch.png")), p, width=5.5, height=4)
  
  cluster.genes <- FindAllMarkersMAESTRO(object = RNA.integrated, only.pos = only.pos, min.pct = genes.pct, test.use = genes.test.use, logfc.threshold = genes.logfc)
  cluster.genes <- cluster.genes[cluster.genes$p_val_adj < genes.cutoff, ]
  write.table(cluster.genes, paste0(RNA.integrated@project.name, "_DiffGenes.tsv"), quote = F, sep = "\t")
  
  return(list(RNA=RNA.integrated, genes=cluster.genes))
}
