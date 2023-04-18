#load libraries

library(glmmSeq)
library(Seurat)
library(qvalue)
library(tidyverse)
library(tidyr)
library(lme4)
library(lmerTest)
library(parallel)
library(edgeR)
library(dplyr)
library(DESeq2)
library(glmGamPoi)
library(scran)

# load rds object
data <- readRDS("data.rds")

# define clusters
cluster_num <- length(unique(combined@meta.data$seuratclusters15))

# define clusters to skip based on low representation
skip_clusters=c(13,14,15)

#run DEG analysis

for (k in 1:15){
  #for (k in 1:1){
  print(paste0("starting cluster ",k-1))
  if(k %in% skip_clusters){
    next
  }
  data <- subset(x = combined, subset = seuratclusters15 == k-1)
  
  df_counts <- data.frame(t(as.matrix(data@assays$RNA@counts)))
  
  n_genes=ncol(df_counts)
  n_cells=nrow(df_counts)
  
  df_counts_meta = data.frame(rownames(df_counts))
  df_counts_meta$id <- df_counts_meta$rownames.df_counts.
  df_counts_meta$rownames.df_counts. = NULL
  df_counts_meta$subject = data$subject  
  df_counts_meta$cond = as.factor(data$cond)
  df_counts_meta$log_spawn_events = as.numeric(data$log_spawn_events)
  df_counts_meta$gsi = as.numeric(data$gsi)  
  df_counts_meta$run= as.factor(data$run)
  df_counts_meta$sample = as.factor(data$sample)
  df_counts_meta$pair= as.factor(data$pair)
  
  # cut out genes that were not detected in any cells from dataframe 
  df_counts_no_0 <- df_counts[,which(colSums(df_counts) != 0)]
  
  # store list of detected genes
  n_genes_no_0=ncol(df_counts_no_0)
  
  # add metadata back into the dataframe
  df_counts_no_0 <- cbind(df_counts_no_0,df_counts_meta)
  
  # split based on pair
  df_counts_no_0_split_by_pair <- split(df_counts_no_0, f = df_counts_no_0$pair) #split by pair
  
  #remove genes with 0 counts within each pair  
  for(l in 1:length(df_counts_no_0_split_by_pair)){
    print(paste0("finding good genes for pair ",l))
    temp_pair_l <- data.frame(df_counts_no_0_split_by_pair[[l]])
    temp_pair_l_counts <- temp_pair_l[,1:n_genes_no_0]
    temp_pair_l_counts_no_0 <- temp_pair_l_counts[,which(colSums(temp_pair_l_counts) != 0)]
    #temp <- cbind(good_genes, temp[20407:20410])
    out=data.frame(colnames(temp_pair_l_counts_no_0))
    assign(x=paste0("gene_list_pair_",l),value=get("out"))
  }
  
  # generate list of genes that were detected in all pairs
  good_gene_list <- gene_list_pair_1$colnames.temp_pair_l_counts_no_0. 
  for(m in 2:length(df_counts_no_0_split_by_pair)){
    print(paste0("integrating good genes from pair ",m," of 19"))
    temp_good_gene_list_m <- data.frame(value=get(paste0("gene_list_pair_",m)))
    temp_good_gene_list_m <- temp_good_gene_list_m$colnames.temp_pair_l_counts_no_0.
    good_gene_list <- intersect(good_gene_list,temp_good_gene_list_m)
  }
  
  p <- length(good_gene_list)
  
  # create a count dataframe for genes that were detected in all pairs
  df_counts_no_0_all_pairs <- df_counts_no_0[, good_gene_list]
  
  # create a count matrix for genes that were detected in all pairs  
  count_matrix_final <- as.matrix(df_counts_no_0_all_pairs)
  count_matrix_final <- as.data.frame(t(count_matrix_final))
  
  # add metadata into the dataframe for genes that were detected in all pairs   
  df_counts_no_0_all_pairs <- cbind(df_counts_no_0_all_pairs,df_counts_meta)
  
  #set up coldata for DESeq2
  coldata <- df_counts_meta
  coldata$subject <- as.factor(coldata$subject)
  coldata$cond <- as.factor(coldata$cond)
  coldata$log_spawn_events <- as.numeric(coldata$log_spawn_events)
  coldata$gsi <- as.numeric(coldata$gsi)
  coldata$sample <- as.factor(coldata$sample)
  coldata$run <- as.factor(coldata$run)
  coldata$pair <- as.factor(coldata$pair)
  
  cluster_size <- ncol(ncol(count_matrix_final))
  
  dds <- DESeqDataSetFromMatrix(countData = count_matrix_final,
                                colData = coldata,
                                design = ~ run + gsi + log_spawn_events + cond)
  
  size_factors <- calculateSumFactors(
    count_matrix_final,
    #sizes = seq(27, 101, 5),
    clusters = NULL,
    ref.clust = NULL,
    max.cluster.size = cluster_size,
    positive = TRUE,
    scaling = NULL,
    min.mean = NULL,
    subset.row = NULL,
    #BPPARAM = SerialParam()
  )                              
  
  sizeFactors(dds) <- size_factors
  
  dds <- estimateDispersions(
    dds,
    fitType = "glmGamPoi",
    #method = "runed",
    #sharingMode = 'maximum',
    useCR = TRUE,
    maxit = 100,
    weightThreshold = 0.01,
    quiet = FALSE,
    modelMatrix = NULL,
    minmu = 1e-06
  )
  
  disp <- as.matrix(mcols(dds))
  disp <- disp[,11]
  names(disp) <- good_gene_list
  
  glmmseq_counts <- data.frame(t(as.matrix(df_counts_no_0_all_pairs)))
  glmmseq_counts <- glmmseq_counts[1:p,]
  glmmseq_counts <- as.data.frame(sapply(glmmseq_counts, as.numeric))
  rownames(glmmseq_counts) <- colnames(df_counts_no_0_all_pairs[,1:p])
  
  results <- glmmSeq(~ cond + log_spawn_events + gsi + (1|run/sample/subject) + (1|pair/subject),
                     id = "subject",
                     countdata = glmmseq_counts,
                     metadata = coldata,
                     dispersion = disp,
                     removeDuplicatedMeasures = FALSE,
                     removeSingles=FALSE,
                     progress=TRUE,
                     cores = 3
  )
  
  out_final = data.frame(results@stats)
  out_final$cluster=k-1
  write.csv(out_final, paste0("results_cluster",k-1,".csv"))
} 