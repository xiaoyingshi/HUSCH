### basic QC module
#nfeature density plot
nfeature_plot <- function(Seurat){
  nFeature <- as.numeric(Seurat$nFeature_RNA)
  median <- round(median(nFeature),0)
  nFeature <- as.data.frame(nFeature)
  colnames(nFeature) <- "nFeature"
  plot <- ggplot(data = nFeature) +
    geom_density(aes(x = nFeature),color = "#6ecdb7") +
    geom_vline(xintercept = median,linetype="dashed", size=0.5)+
    labs(title = "nFeature", subtitle = paste0("median = ",median))+ 
    theme_bw()+
    theme(plot.title = element_text(hjust=0.5,vjust = 0.5, margin = margin(l=50,r=50,t=10,b=10), face = "bold", colour = "black"),
        plot.subtitle =  element_text(hjust=0.5,vjust = 0.5, margin = margin(l=50,r=50,t=0,b=5), colour = "black"),
        axis.text.x=element_text(size=10,face = "bold", hjust = 0.5,vjust = 0.5),
        axis.text.y=element_text(size=10,face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"))
    return(list(plot,median))
}

#nCount density plot
ncount_plot <- function(Seurat){
  nCount <- as.numeric(Seurat$nCount_RNA)
  median <- round(median(nCount),0)
  nCount <- as.data.frame(nCount)
  colnames(nCount) <- "nCount"
  plot <- ggplot(data = nCount) +
    geom_density(aes(x = nCount),color = "#6ecdb7") +
    geom_vline(xintercept = median,linetype="dashed", size=0.5)+
    labs(title = "nCount", subtitle = paste0("median = ",median))+ 
    theme_bw()+
    theme(plot.title = element_text(hjust=0.5,vjust = 0.5, margin = margin(l=50,r=50,t=10,b=10), face = "bold", colour = "black"),
        plot.subtitle =  element_text(hjust=0.5,vjust = 0.5, margin = margin(l=50,r=50,t=0,b=5), colour = "black"),
        axis.text.x=element_text(size=10,face = "bold", hjust = 0.5,vjust = 0.5),
        axis.text.y=element_text(size=10,face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"))
    return(list(plot,median))
}

#DEgene number bar plot
DE_gene_plot <- function(DEG){
  DEG <- DEG[which(abs(DEG$avg_logFC)>=0.25),]
  DEG <- DEG %>% mutate(Dif=ifelse(avg_logFC>=0.25,"up","down"))
  cluster_sorted <- unique(as.numeric(DEG$cluster))
  cluster_sorted <- cluster_sorted[order(cluster_sorted)]
  Palette <- c("blue","red")
  up_mean <- round(length(which(DEG$dif=="up"))/length(cluster_sorted),0)
  down_mean <- round(length(which(DEG$dif=="down"))/length(cluster_sorted),0)
  plot <- ggplot(DEG)+
    geom_bar(aes(x=factor(cluster,levels = cluster_sorted),fill=Dif),position = "dodge")+
    scale_fill_manual(values=Palette)+
    labs(x="Cluster",title = "DEgenes_number",subtitle = paste0("up mean= ",up_mean," | ","down mean = ",down_mean))+
    theme_bw()+
    theme(plot.title = element_text(hjust=0.5,vjust = 0.5, margin = margin(l=50,r=50,t=10,b=10),face = "bold", colour = "black"),
          plot.subtitle =  element_text(hjust=0.5,vjust = 0.5, margin = margin(l=50,r=50,t=0,b=5), colour = "black"),
          axis.text.x=element_text(size=10,face = "bold", hjust = 0.5,vjust = 0.5),
          axis.text.y=element_text(size=10,face = "bold"),
          axis.title.x = element_text(size = 12, face = "bold"),
          axis.title.y = element_text(size = 12, face = "bold"))
    return(list(plot,up_mean))
}

basic_qc_main <- function(file_path){
  # file_path="/Users/morsouron/Desktop/rds/database/pbmc_healthy_9519_res.rds"
  Sname <- gsub("_res.rds","",basename(file_path))
  message(paste0("[Info] Running basic_qc module for : ",Sname))
  rds <- readRDS(file_path)
  Seurat <- rds$RNA
  DEG <- rds$genes

  #basic information
  ncell <- dim(Seurat)[2]
  ncluster <- length(unique(Seurat$seurat_clusters))

  #UMAP plot of batch
  if (BatchName %in% colnames(Seurat@meta.data)){
    Batch_plot <- DimPlot(Seurat,group.by = BatchName,label = TRUE,pt.size = 0.1,label.size = 4,repel = TRUE)+labs(title = "Batch UMAP")+theme(plot.title = element_text(hjust=0.5,vjust = 0.5, margin = margin(l=100,r=50,t=10,b=10),face = "bold", colour = "black"))
  }

  #basic information plot
  nFeature_P <- nfeature_plot(Seurat)[[1]]
  nCount_P <- ncount_plot(Seurat)[[1]]
  DE_gene_P <- DE_gene_plot(DEG)[[1]]

  #statistic summary
  nfeature_median <- nfeature_plot(Seurat)[[2]]
  nCount_median <- ncount_plot(Seurat)[[2]]
  DEG_up_mean <- DE_gene_plot(DEG)[[2]]

  #estimate data quality
  Ncell <- ifelse(ncell>=1000,1,0)
  Nfeature <- ifelse(nfeature_median>=600,1,0)
  Ncount <- ifelse(nCount_median>=1500,1,0)
  Data_quality <- Ncell+Nfeature+Ncount

  #export plot
  if (BatchName %in% colnames(Seurat@meta.data)){
    plot <- list(Batch_plot,DE_gene_P,nFeature_P,nCount_P)
    plot_size <- 500+100*ncluster
    png(paste0(outpath,"/",Sname,"_qc.png"),res = 300,height = 3000,width =1800+plot_size)
    qc_plot <- ggarrange(plotlist = plot, ncol = 2, nrow = 2, widths = c(1800,plot_size), heights = c(1600,1400))
    print(annotate_figure(qc_plot,top=text_grob(paste0(Sname,"_QC"),face = "bold", size = 18)))
    dev.off()
  } else {
    plot <- list(nFeature_P,nCount_P,DE_gene_P)
    plot_size <- 500+100*ncluster
    if (plot_size < 1500){
      pwidth <- 1500
    } else {
      pwidth <- plot_size
    }
    png(paste0(outpath,"/",Sname,"_qc.png"),res = 300,height = 3000,width = pwidth*2)
    qc_plot <- ggarrange(plotlist = plot, ncol = 2, nrow = 2, widths = c(pwidth,pwidth), heights = c(1500,1500))
    print(annotate_figure(qc_plot,top=text_grob(paste0(Sname,"_QC"),face = "bold", size = 18)))
  }

  return(c(Sname,ncell,ncluster,nfeature_median,nCount_median,DEG_up_mean,Ncell,Nfeature,Ncount,Data_quality))
}