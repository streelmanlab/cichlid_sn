#load libraries
library(Seurat)
library(qvalue)
library(tidyverse)
library(lme4)
library(PROreg)

#load data
data <- readRDS("data.rds")

cluster_num <- length(unique(combined@meta.data$seuratclusters15))

#define the genes used to calculate a gene score (the number of genes in the list expressed in each nucleus or cell)
genes <- as.vector("gene_name_1", "gene_name_2", "gene_name_3", [...], "gene_name_n")

# define n as number of unique genes used to create gene score
n <- length(genes)

#add gene score column
createGeneScore = function(data, genes, new.col.name = "gene_score", normalize.by.n.features = F) {
  #' Create a gene score: the # of genes from the input expressed in each cell.
  #' By default, the gene score is not normalized by the # of features in each cell.
  #'
  #' @param data Seurat object
  #' @param genes vector of genes
  #' @param new.col.name name of meta.data column to store the score in
  #' @param normalize.by.n.features normalize gene score by # of features in the cell?
  
  # Change Matrix to Prescence/Absence
  mat = data@assays$RNA@counts[genes,]
  mat[which(mat > 1)] = 1
  
  # Add the new column
  data@meta.data[, new.col.name] = colSums(mat)
  
  # Normalize by # of Features
  if (normalize.by.n.features) { data@meta.data[, new.col.name] = data@meta.data[, new.col.name] / data$nFeature_RNA }
  
  return(data)
}



# define clusters to skip based on low representation
skip_clusters=c(13,14,15)

# run differential gene score analysis

bbmm_out=data.frame(est_bbmm_interc=integer(),
                    st_err_bbmm_interc=integer(),
                    t_bbmm_interc=integer(),
                    p_bbmm_interc=integer(),
                    est_bbmm_cond=integer(),
                    st_err_bbmm_cond=integer(),
                    t_bbmm_cond=integer(),
                    p_bbmm_cond=integer(),
                    est_bbmm_log_spawn_events=integer(),
                    st_err_bbmm_log_spawn_events=integer(),
                    t_bbmm_log_spawn_events=integer(),
                    p_bbmm_log_spawn_events=integer(),
                    est_bbmm_gsi=integer(),
                    st_err_bbmm_gsi=integer(),
                    t_bbmm_gsi=integer(),
                    p_bbmm_gsi=integer(),
                    cluster=integer(),
                    stringsAsFactors=FALSE)

for (k in 1:15){
  print(k)
  if(k %in% skip_clusters){
    next
  }
  
  data <- subset(x = combined, subset = seuratclusters15 == k-1)
  df <- data.frame(data$gene_score)
  df$cond = as.factor(data$cond)
  df$log_spawn_events = as.numeric(data$log_spawn_events)
  df$gsi = as.numeric(data$gsi)
  df$id = rownames(df)
  df$run = as.factor(data$run)
  df$sample = as.factor(data$sample)
  df$pair= as.factor(data$pair)
  df$subject = as.factor(data$subject)
  
  df_split <- split(df, f = df$subject)
  subject = as.factor(df$subject)
  cond = as.factor(df$cond)
  gsi = as.numeric(df$gsi)
  pair = as.factor(df$pair)
  run = as.factor(df$run)
  
  gene_score <- data$gene_score
  
  bbmm <- BBmm(fixed.formula = gene_score ~ as.numeric(cond) + as.numeric(log_spawn_events) + as.numeric(gsi), random.formula = ~ (subject %in% sample %in% run) + (subject %in% pair) , m=n, data = df, show = TRUE)
  bbmm_summary <- summary(bbmm)
  bbmm_summary <- data.frame(bbmm_summary$fixed.coefficients)
  
  bbmm_out[k,1]=bbmm_summary[1,1]
  bbmm_out[k,2]=bbmm_summary[1,2]
  bbmm_out[k,3]=bbmm_summary[1,3]
  bbmm_out[k,4]=bbmm_summary[1,4]
  bbmm_out[k,5]=bbmm_summary[2,1]
  bbmm_out[k,6]=bbmm_summary[2,2]
  bbmm_out[k,7]=bbmm_summary[2,3]
  bbmm_out[k,8]=bbmm_summary[2,4]
  bbmm_out[k,9]=bbmm_summary[3,1]
  bbmm_out[k,10]=bbmm_summary[3,2]
  bbmm_out[k,11]=bbmm_summary[3,3]
  bbmm_out[k,12]=bbmm_summary[3,4]
  bbmm_out[k,13]=bbmm_summary[4,1]
  bbmm_out[k,14]=bbmm_summary[4,2]
  bbmm_out[k,15]=bbmm_summary[4,3]
  bbmm_out[k,16]=bbmm_summary[4,4]
  bbmm_out[k,17]=cluster=paste0(k-1)
  
}

write.csv(bbmm_out, 
          "results.csv")

library(qvalue)
qobj <- qvalue(p = bbmm_out$p_bbmm_cond)
bbmm_out$bh <- p.adjust(bbmm_out$p_bbmm_cond, method = "BH")

bbmm_out$q <- qobj$qvalues
write.csv(bbmm_out, 
          "results_q.csv")
