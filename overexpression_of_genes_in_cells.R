
markerExpPerCellPerClusterQuick = function(obj, markers, include_stars = F, pt.alpha = 0.05) {
  #' This function tests for non-uniform expression of a list of genes across cells.
  #' The identity of the object are the groupings that are tested for non-uniform
  #' distribution. The function returns plots and dataframes with p-values and effect
  #' sizes for the overexpression of the list of genes across groups of cells.
  #' 
  #' @import BSDA
  #' @import ggplot2
  #' @import effsize
  #' 
  #' @param obj Seurat object
  #' @param markers vector of genes
  #' @param include_stars include stars denoting significance on the plots?
  #' @param pt.alpha alpha (transparancy) value of points on plots
  #' 
  #' @returns list of 4 elements:
  #' 1) boxplot for each group of cells on the x-axis and enrichment for the markers on the y-axis
  #' 2) plot of the effect size of the enrichment
  #' 3) dataframe for 1
  #' 4) dataframe for 2
  
  n_markers = T # use the number of markers expressed per cell as the score instead of sum of the expression of markers
  correct = T # correct for the total number of genes expressed per cell
  myslot = "counts"
  
  title_str = "Average"
  title_str = paste(title_str, if_else(myslot=="counts", "", "Normalized"))
  title_str = paste(title_str, "Expression of Marker Genes per Cell per Cluster")
  
  exp = GetAssayData(obj, assay = "RNA", slot=myslot)
  exp = exp[markers,]
  exp[which(exp > 0)] = 1
  per_cluster_df = data.frame()
  d_df = data.frame()
  clusters = levels(Idents(obj))
  if (! any(is.na(as.numeric( levels(Idents(obj)) )))) {
    clusters = sort(unique(as.numeric(as.vector(Idents(obj)))))
  }
  
  for (cluster in clusters) {
    cluster_cells <- WhichCells(obj, idents = cluster)
    
    avg_cluster_exp = colSums(exp[markers, cluster_cells])
    avg_cluster_exp = avg_cluster_exp/obj$nFeature_RNA[cluster_cells]
    other_cells = colnames(obj)[which(! colnames(obj) %in% cluster_cells)]
    other_avg = colSums(exp[markers, other_cells])
    other_avg = other_avg/obj$nFeature_RNA[other_cells]
    
    this.z.test = z.test(avg_cluster_exp, other_avg, sigma.x = sd(avg_cluster_exp), sigma.y = sd(other_avg), alternative = "greater")
    p    = this.z.test$p.value
    stat = this.z.test$statistic
    
    # Cohen's d
    all_exp = c(avg_cluster_exp, other_avg)
    test= effsize::cohen.d(all_exp, c(rep("cluster", length(avg_cluster_exp)),
                                      rep("other",   length(other_avg))))
    
    d=test$estimate
    up=test$conf.int[2]
    down = test$conf.int[1]
    mag=test$magnitude
    mag_pos=mag
    mag_pos[which(d < 0)] = "negligible"
    
    d_df = rbind(d_df, data.frame(cluster, mag, mag_pos, d, up, down, p, stat))
    per_cluster_df = rbind(per_cluster_df, data.frame(cluster, avg_cluster_exp, p, mag_pos))
  }
  
  d_df$p_val_adj = p.adjust(d_df$p, method = "bonferroni")
  d_df$cluster = factor(d_df$cluster, levels = clusters)
  p1  = ggplot(d_df, aes(cluster, d, color=mag, fill = mag)) + geom_pointrange(aes(ymin = down, ymax = up)) + ylab("Cohen's d") + ggtitle("Effect Size in Clusters")
  
  colnames(per_cluster_df) <- c("cluster", "avg_cluster_exp", "p", "mag_pos")
  per_cluster_df$avg_cluster_exp = as.numeric(as.vector(per_cluster_df$avg_cluster_exp))
  per_cluster_df$p_val_adj = unlist(sapply(1:length(clusters), function(x) rep(d_df$p_val_adj[which(d_df$cluster == clusters[x])], length(which(per_cluster_df$cluster == clusters[x]))) ))
  per_cluster_df$star = ""
  if (include_stars) {
    per_cluster_df$star = ifelse(per_cluster_df$p_val_adj < 0.001, "***",
                                 ifelse(per_cluster_df$p_val_adj < 0.01, "**",
                                        ifelse(per_cluster_df$p_val_adj < 0.05, "*", "")))
  }
  per_cluster_df$cluster = factor(per_cluster_df$cluster, levels = clusters)
  per_cluster_df$mag_pos = factor(per_cluster_df$mag_pos, levels = c("negligible", "small", "medium", "large"))
  p = ggplot(per_cluster_df, aes(cluster, avg_cluster_exp, fill=mag_pos, color=mag_pos)) + geom_text(aes(x= cluster, y = Inf, vjust = 1, label = star)) + geom_boxplot(alpha=0.6) +  geom_jitter(position=position_dodge2(width = 0.6), alpha = pt.alpha) + xlab("Cluster") + ylab("") + ggtitle(title_str) + scale_color_viridis(discrete = T,drop=TRUE, limits = levels(per_cluster_df$mag_pos), name = "Effect Size") + scale_fill_viridis(discrete = T, drop=TRUE, limits = levels(per_cluster_df$mag_pos), name = "Effect Size")
  
  return(list(p, p1, d_df, per_cluster_df))
}
