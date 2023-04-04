rna_path <-"<rna_path>" 

b1.data <- Read10X(data.dir = paste(rna_path, "/JTS07-B1/outs/filtered_feature_bc_matrix/", sep=""))
b2.data <- Read10X(data.dir = paste(rna_path, "/JTS07-B2/outs/filtered_feature_bc_matrix/", sep=""))
b3.data <- Read10X(data.dir = paste(rna_path, "/JTS07-B3/outs/filtered_feature_bc_matrix/", sep=""))
b4.data <- Read10X(data.dir = paste(rna_path, "/JTS07-B4/outs/filtered_feature_bc_matrix/", sep=""))
b5.data <- Read10X(data.dir = paste(rna_path, "/JTS07-B5/outs/filtered_feature_bc_matrix/", sep=""))
c1.data <- Read10X(data.dir = paste(rna_path, "/JTS07-C1/outs/filtered_feature_bc_matrix/", sep=""))
c2.data <- Read10X(data.dir = paste(rna_path, "/JTS07-C2/outs/filtered_feature_bc_matrix/", sep=""))
c3.data <- Read10X(data.dir = paste(rna_path, "/JTS07-C3/outs/filtered_feature_bc_matrix/", sep=""))
c4.data <- Read10X(data.dir = paste(rna_path, "/JTS07-C4/outs/filtered_feature_bc_matrix/", sep=""))
c5.data <- Read10X(data.dir = paste(rna_path, "/JTS07-C5/outs/filtered_feature_bc_matrix/", sep=""))

b1 <- CreateSeuratObject(counts = b1.data, project = "b1")
b2 <- CreateSeuratObject(counts = b2.data, project = "b2")
b3 <- CreateSeuratObject(counts = b3.data, project = "b3")
b4 <- CreateSeuratObject(counts = b4.data, project = "b4")
b5 <- CreateSeuratObject(counts = b5.data, project = "b5")

c1 <- CreateSeuratObject(counts = c1.data, project = "c1")
c2 <- CreateSeuratObject(counts = c2.data, project = "c2")
c3 <- CreateSeuratObject(counts = c3.data, project = "c3")
c4 <- CreateSeuratObject(counts = c4.data, project = "c4")
c5 <- CreateSeuratObject(counts = c5.data, project = "c5")


b1$sample <- "b1"
b2$sample <- "b2"
b3$sample <- "b3"
b4$sample <- "b4"
b5$sample <- "b5"

c1$sample <- "c1"
c2$sample <- "c2"
c3$sample <- "c3"
c4$sample <- "c4"
c5$sample <- "c5"


b1$cond <- "BHVE"
b2$cond <- "BHVE"
b3$cond <- "BHVE"
b4$cond <- "BHVE"
b5$cond <- "BHVE"

c1$cond <- "CTRL"
c2$cond <- "CTRL"
c3$cond <- "CTRL"
c4$cond <- "CTRL"
c5$cond <- "CTRL"


b1 <- NormalizeData(b1, normalization.method = "LogNormalize")
b2 <- NormalizeData(b2, normalization.method = "LogNormalize")
b3 <- NormalizeData(b3, normalization.method = "LogNormalize")
b4 <- NormalizeData(b4, normalization.method = "LogNormalize")
b5 <- NormalizeData(b5, normalization.method = "LogNormalize")

c1 <- NormalizeData(c1, normalization.method = "LogNormalize")
c2 <- NormalizeData(c2, normalization.method = "LogNormalize")
c3 <- NormalizeData(c3, normalization.method = "LogNormalize")
c4 <- NormalizeData(c4, normalization.method = "LogNormalize")
c5 <- NormalizeData(c5, normalization.method = "LogNormalize")


#i've already checked qc parameters for each sample individually, because they are all similar merging at this step for simplicity
combined <- merge(x=b1, y=c(b2,b3,b4,b5,c1,c2,c3,c4,c5), merge.data = TRUE, add.cell.ids = c("BHVE", "BHVE", "BHVE", "BHVE", "BHVE", "CTRL", "CTRL", "CTRL", "CTRL", "CTRL"))

# quality control
combined <- subset(combined, subset = nFeature_RNA > 300)
combined <- subset(combined, subset = nFeature_RNA < 3000)
combined <- subset(combined, subset = nCount_RNA > 500)
combined <- subset(combined, subset = nCount_RNA < 8000)
mt_genes = scan("<mt_gene_list.txt>", what = "character")
mt_genes = str_replace(mt_genes, "_", "-")
mt_genes = mt_genes[which(mt_genes %in% rownames(combined))]
combined$pct_mt = PercentageFeatureSet(combined, features = mt_genes)
combined <- subset(combined, subset = pct_mt < 5)

#find variable features
combined <- FindVariableFeatures(object = combined, selection.method = "vst", nfeatures = 4000, verbose = TRUE)

# dimensionality reduction
combined <- ScaleData(combined, verbose = TRUE)
combined <- RunPCA(combined, dim = 50, verbose = TRUE)
combined <- RunUMAP(combined, dims = 1:50, min.dist = 0.5, n.neighbors = 50, metric = "euclidean")
combined <- FindNeighbors(combined, reduction = "umap", k.param = 50 , dims = 1:2, prune.SNN = 0)

# primary clustering resolution
combined <- FindClusters(combined, resolution = 0.01, algorithm = 2) #15 clusters

# secondary clustering resolution
combined <- FindClusters(combined, resolution = 1.3, algorithm = 2) #53 clusters
