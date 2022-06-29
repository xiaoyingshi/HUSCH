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




##function of run batch effect removal by harmony
harmony <- function(rds, batch, nfeatures = 2000, dims.use = 1:15, cluster.res = 0.6, only.pos = FALSE, genes.test.use = "presto",
                    genes.cutoff = 1E-5, genes.pct = 0.1, genes.logfc = 0.25
                    ){
  cells <- ncol(rds$RNA)
  if (cells <= 5000) {
    npc <- 50; cluster.res <- 0.6
  } else if(cells <= 100000) {
    npc <- 50; cluster.res <- 1
  } else {
    npc <- 100; cluster.res <- 1
  }
  pc.contribution <- rds$RNA@reductions$pca@stdev / sum(rds$RNA@reductions$pca@stdev) * 100
  pc.contribution.cum <- cumsum(pc.contribution)
  pc.first <- which(pc.contribution.cum > 75)[1]
  dims.use = 1:pc.first 
  RNA.harmony <- RunHarmony(rds$RNA,batch) %>% 
    RunUMAP(reduction = "harmony", dims = dims.use) %>% 
    FindNeighbors(reduction = "harmony", dims = dims.use) %>% 
    FindClusters(resolution = cluster.res)
  RNA.harmony@project.name <- paste0(rds$RNA@project.name,'_harmony')
  
  p = DimPlot(object = RNA.harmony, label = TRUE, pt.size = 0.2)
  ggsave(file.path(paste0(RNA.harmony@project.name, "_cluster.png")), p, width=5, height=4)
  p = DimPlot(object = RNA.harmony, group = batch, label = TRUE, pt.size = 0.2)
  ggsave(file.path(paste0(RNA.harmony@project.name,'_',batch,"_batch.png")), p, width=5.5, height=4)
  
  cluster.genes <- FindAllMarkersMAESTRO(object = RNA.harmony, only.pos = FALSE, min.pct = genes.pct,  logfc.threshold = genes.logfc)
  cluster.genes <- cluster.genes[cluster.genes$p_val_adj < genes.cutoff, ]
  write.table(cluster.genes, paste0(RNA.harmony@project.name, "_DiffGenes.tsv"), quote = F, sep = "\t")
  
  return(list(RNA=RNA.harmony, genes=cluster.genes))
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