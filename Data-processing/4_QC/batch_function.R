CalEntropy <- function(Seurat, BatchName, Sname) {
  message(paste0("[Info] Calculate Entropy for : ", Sname))
  knn <- kNN(Seurat@reductions[["umap"]]@cell.embeddings, k = 30)
  neighbor_matrix <- knn$id
  #get a vector of batch conditions: ["patient1","patient1","patient2"...]
  Vbatch <- Seurat@meta.data[, BatchName]
  #batch compositions: ["patient1","patient2","patient3"...]
  batch <- unique(Vbatch)
  Entropies <- c()
  for (i in 1:ncol(Seurat)) {
    sum = 0
    cell_neighbors <- as.character(unlist(lapply(neighbor_matrix[i,], function(x) {
      x <- Vbatch[x]
    })))
    for (j in 1:length(batch)) {
      Nbatch <- length(which(cell_neighbors == batch[j]))
      if (Nbatch == 0) {
        sum <- sum + 0
      } else {
        ratio <- Nbatch / 30
        sum <- sum + ratio * log2(ratio)
      }
    }
    Entropies[i] <- -sum
  }
  # return entropies(a vector of number): [1,0.5,0.3 ...]
  return(Entropies)
}

batch_correct <- function(file_path) {
  RNA <- readRDS(file_path)$RNA
  Sname <- gsub("_res.rds", "", basename(file_path))
  message(paste0("[Info] Run batch correction for : ", Sname))
  # choose npc
  cells <- ncol(RNA)
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
  pc.contribution <- RNA@reductions$pca@stdev / sum(RNA@reductions$pca@stdev) * 100
  pc.contribution.cum <- cumsum(pc.contribution)
  pc.first <- which(pc.contribution.cum > 75)[1]
  dims.use = 1:pc.first

  data.list <- SplitObject(RNA, split.by = BatchName)
  for (i in 1:length(data.list)) {
    data.list[[i]] <- NormalizeData(data.list[[i]], verbose = FALSE)
    data.list[[i]] <- FindVariableFeatures(data.list[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  }
  anchors <- FindIntegrationAnchors(object.list = data.list, dims = dims.use, anchor.features = 2000)
  RNA.integrated <- IntegrateData(anchorset = anchors, dims = dims.use)
  RNA.integrated@project.name <- paste0(Sname, '_CCA')
  DefaultAssay(RNA.integrated) <- "integrated"

  RNA.integrated <- ScaleData(RNA.integrated, verbose = FALSE)
  RNA.integrated <- RunPCA(RNA.integrated, npcs = npc, verbose = FALSE, features = VariableFeatures(RNA.integrated))
  p = ElbowPlot(object = RNA.integrated, ndims = RNA.integrated@commands$RunPCA.integrated@params$npcs)
  ggsave(paste0(dirname(file_path), "/", Sname, "_CCA_PCElbowPlot.png"), p, width = 10, height = 4)
  RNA.integrated <- RunUMAP(object = RNA.integrated, reduction = "pca", dims = dims.use)

  RNA.integrated <- FindNeighbors(object = RNA.integrated, reduction = "pca", dims = dims.use)
  RNA.integrated <- FindClusters(object = RNA.integrated, resolution = cluster.res)

  p = DimPlot(object = RNA.integrated, group = BatchName, label = FALSE, pt.size = 0.2)
  ggsave(paste0(dirname(file_path), "/", Sname, "_", BatchName, "_CCA_.png"), p, width = 5.5, height = 4)

  cluster.genes <- FindAllMarkersMAESTRO(object = RNA.integrated, only.pos = FALSE, min.pct = 0.1, test.use = "presto", logfc.threshold = 0.25)
  cluster.genes <- cluster.genes[cluster.genes$p_val_adj < 1E-5,]
  write.table(cluster.genes, paste0(dirname(file_path), "/", Sname, "_CCA_DiffGenes.tsv"), quote = F, sep = "\t")
  saveRDS(list(RNA = RNA.integrated, genes = cluster.genes), file = paste0(dirname(file_path), "/", Sname, "_CCA_res.rds"))
}

batch_main <- function(file_path) {

  #################################################################################
  # Sname=sample name;Nbatch=number of batches;flag=the result of batch detecting #
  #################################################################################

  Seurat <- readRDS(file_path)$RNA
  Sname <- Seurat@project.names
  BatchName <- batchInfo$batchName[BatchName$rds_path==file_path]
  if (BatchName %in% colnames(Seurat@meta.data)) {
    Nbatch <- length(unique(Seurat@meta.data[, BatchName]))
    MaxEntropy <- log2(Nbatch)
    if (Nbatch >= 2) {
      Entropy <- CalEntropy(Seurat, BatchName, Sname)
      if (MaxEntropy / median(Entropy) >= 4) {
        flag <- "Yes"
        batch_correct(file_path)
      } else {
        flag <- "No"
      }
    }
    else {
      flag <- "Only one batch"
    }
  }
  else {
    flag <- "No batch information"
  }
  ### final result
  if (flag == "Yes") {
    return(c(Sname, flag, paste0(dirname(file_path), "/", Sname, "_CCA_res.rds")))
  } else {
    return(c(Sname, flag, file_path))
  }
}

