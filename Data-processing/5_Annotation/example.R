suppressMessages(library(purrr))
suppressMessages(library(Seurat))
suppressMessages(library(sva))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(harmony))
suppressMessages(library(MAESTRO))
suppressMessages(library(stringr))

## 5.2 Generate signature list
Marker <- read.table(marker_path, sep = "\t", stringsAsFactors = F, header = T)
signature_list <- list()
for (i in 1:nrow(Marker)) {
  signature_list[Marker$Celltype[i]] <- list(unique(strsplit(Marker$Marker[i], ",")[[1]]))
}

## 5.3 Annotate with function
source('5_Annotate/selina_data_curate.R')
marker_res=MarkerAnnot(rds$RNA, rds$genes, clusterid=clusterid, signatures = signature_list, min.score = 0)
result_df$marker_res <- marker_res$cellscore_max_celltype[match(rownames(marker_res), result_df$cl)]

## 5.4 Use featureplot to check marker expression and curation. If some cluster's annotation are weird, replace it with the right one.

options(repr.plot.height=10,repr.plot.width=15)
DimPlot(rds$RNA, group.by='seurat_clusters', label = TRUE, pt.size = 0.1)

allmarker <- marker_list
allmarker <- unlist(allmarker, use.names=F)
allmarker <- unique(allmarker)
options(repr.plot.height=10,repr.plot.width=15)
FeaturePlot(rds$RNA, features = allmarker, slot = 'data')

# Compare seurat clusters plot with FeaturePlot
rds$RNA$manual.assign = 'assign'
rds$RNA$manual.assign[which(rds$RNA$seurat_clusters %in% c(cluster_id))] = 'newCellType'

options(repr.plot.height=6,repr.plot.width=8)
DimPlot(rds$RNA, group.by='manual.assign', label = TRUE, pt.size = 0.1)

# check cell number
dim(rds$RNA)
rds$RNA@project.name
folder_path <- dirname(file)

# name file with new name
rds$RNA@project.name <- 'GSEnumber_Tissue-Stage_cellnumber'

# save file
saveRDS(rds, paste0(folder_path, '/', rds$RNA@project.name, '.rds'))


## 5.5 Unify into different level
source("5_Annotation/selina_data_curate.R")

# all file in one tissue
files <- list.files("/fs/home/shixiaoying/Project/SELINA/Uterus/Adult", pattern = "rds$", full.names = T, recursive = T)
files

data_list <- map(files, readRDS)
for(i in 1:length(data_list)){
  if(length(data_list[[i]])==1){
    data_list[[i]]=data_list[[i]]
  }else{
    data_list[[i]]=data_list[[i]]$RNA
  }   
}

allLevel4CT <- c()
for (i in 1:length(data_list)){
  print(i)
  #     print(unique(data_list[[i]] $manual.assign))
  print(unique(sort(data_list[[i]]$curated.assign)))
  allLevel4CT <- c(allLevel4CT, unique(data_list[[i]]$curated.assign))
}
allLevel4CT <- data.frame(sort(unique(allLevel4CT)), stringsAsFactors = F)
colnames(allLevel4CT) <- 'manual'

# Filter out <10 cell type
for(i in 1:length(data_list)){
  data_list[[i]]$temp.assign = data_list[[i]]$curated.assign
  if(length(data_list[[i]])==1){
    tb=table(data_list[[i]]$curated.assign)
    for(ct in names(tb)){
      if(tb[ct]<=10){
        data_list[[i]]$curated.assign[data_list[[i]]$curated.assign == ct] = paste0(ct,'_Unknown')
      }
    }
  }else{
    tb=table(data_list[[i]]$RNA$curated.assign)
    for(ct in names(tb)){
      if(tb[ct]<=10){
        data_list[[i]]$RNA$curated.assign[data_list[[i]]$RNA$curated.assign == ct] = paste0(ct,'_Unknown')
      }
    }
  }   
}

# normal unification

for (i in 1:length(rdspath_sub)){
  if(length(grep('manual|Celltype',colnames(data_list[[i]]@meta.data)))>0){
    data_list[[i]]$curated.assign = gsub("_cell*", "", data_list[[i]]$curated.assign, ignore.case = T)
    data_list[[i]]$curated.assign = gsub('s$','',data_list[[i]]$curated.assign)
    data_list[[i]]$curated.assign = gsub(" cell$", "", data_list[[i]]$curated.assign, ignore.case = T)
    data_list[[i]]$curated.assign = sapply(data_list[[i]]$curated.assign, simpleCap)
    #and abbr->full name
    
  }
}

level <- read.table('5_Annotation/alllevel_file/Uterus_level.txt', 
                    sep='\t', stringsAsFactor=F, header=T)
level

# level1
for (i in 1:length(data_list)){
  data_list[[i]]$level1 = data_list[[i]]$curated.assign
  for (j in 1:nrow(level)){
    data_list[[i]]$level1[data_list[[i]]$curated.assign %in% level$manual[level$level1==level[j,'level1']]] = level[j,'level1']
  }
}

# level2
for (i in 1:length(data_list)){
  data_list[[i]]$level2 = data_list[[i]]$temp.assign
  for (j in 1:nrow(level)){
    data_list[[i]]$level2[data_list[[i]]$temp.assign %in% level$manual[level$level2==level[j,'level2']]] = level[j,'level2']
  }  
}

#T cell 
signatures <- read.table('5_Annotation/SELINA_Data_Curated/2_Annotation/Signature/T_level2.txt')
lineage = lineage_extract(lineage_cells = c('T'),level = 'level2')
lineage = suppressWarnings(harmony(lineage,'project'))
data_list = maestro_subtype(seurat_obj = lineage,signature = signatures,Parent = 'T',level = 'level2')

# stromal
lineage = lineage_extract(lineage_cells = c('Stromal'),level = 'level2')
lineage = suppressWarnings(harmony(lineage,'project'))
data_list = stromal_subtype(lineage, stromal_list, Parent='Stromal', 'level2') 

# lymphocyte
lineage = lineage_extract(lineage_cells = c('Lymphocyte'),level = 'level2')
lineage = suppressWarnings(harmony(lineage,'project'))
data_list = lymphocyte_subtype(lineage, Parent='Lymphocyte', 'level2') 

# Muscle
Muscle_list = list()
count = 0
for (i in 1:length(data_list)){
  if ('Muscle' %in% data_list[[i]]$level1){
    count = count+1
    Muscle_list[[count]] = data_list[[i]][,which(data_list[[i]]$level1 == 'Muscle')]
    Muscle_list[[count]]$project = data_list[[i]]@project.name
  }
}
Muscle_obj = Muscle_list %>% reduce(merge)
Muscle_obj_cca = suppressMessages(CCA(Muscle_obj,'project'))
options(repr.plot.height=5,repr.plot.width=8)
DimPlot(Muscle_obj_cca,group.by = 'project',label = T,repel = T)
DimPlot(Muscle_obj_cca,group.by = 'curated.assign',label = T)
DimPlot(Muscle_obj_cca,group.by = 'seurat_clusters',label = T)

# level3
for (i in 1:length(data_list)){
  data_list[[i]]$level3=data_list[[i]]$level2
}

# save results
data_list1 <- map(files, readRDS)
data_list2 <- list()
for (i in 1:length(data_list1)){
  if(length(data_list1[[i]])==1){
    data_list2[[i]]=data_list[[i]]
  }else{
    data_list2[[i]]=list()
    data_list2[[i]]$RNA=data_list[[i]]
    data_list2[[i]]$genes=data_list1[[i]]$genes
  }
  
}
for (i in 1:length(files)){
  saveRDS(data_list2[[i]],files[i])
}

