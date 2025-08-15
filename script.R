setwd("/Users/micoria/Documents/work/概普生物/测试/生信人测试数据/")

library(Seurat)
library(patchwork)
library(DoubletFinder)

mtx <- Read10X_h5(filename = "data/1_filtered_feature_bc_matrix.h5")

seu <- CreateSeuratObject(
  counts = mtx,
  project = "P1",       
  min.cells = 3,
  min.features = 200
)

gene_names <- rownames(seu)
has_human_mt    <- sum(grepl("^MT-", gene_names))      # 人：MT-ND1, MT-CO1...
has_mouse_mt_lo <- sum(grepl("^mt-", gene_names))      # 鼠：mt-Nd1, mt-Co1（常见小写）
has_mouse_mt_Mt <- sum(grepl("^Mt-", gene_names))      # 鼠：Mt-Nd1（另一种写法）
has_ENSG        <- sum(grepl("^ENSG", gene_names))     # 人的 Ensembl ID
has_ENSMUSG     <- sum(grepl("^ENSMUSG", gene_names))  # 鼠的 Ensembl ID

mt_pattern <- NULL
if (has_human_mt > 0) {
  mt_pattern <- "^MT-"
} else if (has_mouse_mt_lo > 0) {
  mt_pattern <- "^mt-"
} else if (has_mouse_mt_Mt > 0) {   
  mt_pattern <- "^Mt-"
} else if (has_ENSG > 100 || has_ENSMUSG > 100) {
  message("检测到 Ensembl ID。建议先把 ID 映射成基因符号再计算线粒体比例。")
} else {
  message("未检测到常见线粒体前缀（MT-/mt-/Mt-），请检查物种/命名。")
}
print(mt_pattern)

ribo_pattern <- NULL
if (sum(grepl("^RPL|^RPS", gene_names)) > 50) {
  ribo_pattern <- "^RPL|^RPS"  # 人：RPL*/RPS*
} else if (sum(grepl("^Rpl|^Rps", gene_names)) > 50) {
  ribo_pattern <- "^Rpl|^Rps"  # 鼠：Rpl*/Rps*
}
print(ribo_pattern)

if (!is.null(mt_pattern)) {
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = mt_pattern)
  cat("使用线粒体模式：", mt_pattern, "\n")
}
if (!is.null(ribo_pattern)) {
  seu[["percent.ribo"]] <- PercentageFeatureSet(seu, pattern = ribo_pattern)
  cat("使用核糖体模式：", ribo_pattern, "\n")
}

f1 <- VlnPlot(seu, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)
png("Figure/vio_QC.png", width = 800, height = 600)
print(f1)
dev.off()

seu <- subset(
  seu,
  subset = nFeature_RNA >= 200 &
    nFeature_RNA <= 6000 &  
    (is.null(percent.mt) | percent.mt <= 15)  
)
## 3) 归一化与降维聚类（SCT 推荐）
DefaultAssay(seu) <- "RNA"
seu <- SCTransform(seu, verbose = FALSE)           # 自动回归 nCount 等，省心
seu <- RunPCA(seu, verbose = FALSE)
ElbowPlot(seu, ndims = 50)

seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu, resolution = 0.4)         # 起步分辨率，后续可 0.3–0.8 探索
seu <- RunUMAP(seu, dims = 1:30)
DimPlot(seu, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + NoLegend()


## ====== 继续你的脚本（从这里往下粘贴即可） ======

## 0) 准备输出目录
if (!dir.exists("Figure")) dir.create("Figure", recursive = TRUE)
if (!dir.exists("Result")) dir.create("Result", recursive = TRUE)

## 1) 去双（用 scDblFinder）
suppressPackageStartupMessages({
  library(scDblFinder)
  library(SingleCellExperiment)
})

# 把 Seurat 转 SCE
sce <- as.SingleCellExperiment(seu)

# 设一个起始 doublet 率（7%），可按需要改 0.05–0.10
set.seed(1234)
sce <- scDblFinder(
  sce,
  dbr = 0.07   # 预估的双细胞比例
)

gc()
seu$scDblFinder.class <- colData(sce)$scDblFinder.class
table(seu$scDblFinder.class)

# 去掉 doublets
seu <- subset(seu, subset = scDblFinder.class == "singlet")

# 去双后重跑聚类/UMAP
seu <- FindNeighbors(seu, dims = 1:30, verbose = FALSE)
seu <- FindClusters(seu, resolution = 0.4, verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:30, verbose = FALSE)

p_umap_clusters <- DimPlot(seu, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + NoLegend()
library(ggplot2)
ggsave("Figure/umap_clusters_after_dbl.png", p_umap_clusters, width = 7, height = 5, dpi = 300)

## 2) 找各聚类 marker（用于手动注释大类）
DefaultAssay(seu) <- "SCT"
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers <- markers[order(markers$cluster, -markers$avg_log2FC), ]
write.csv(markers, "Result/cluster_markers.csv", row.names = FALSE)


## 3) 乳腺肿瘤组织常用 marker 面板（DotPlot + 一些 FeaturePlot）
features_panel <- c(
  # 上皮
  "EPCAM","KRT8","KRT18","KRT19","MUC1",
  # 单核/巨噬
  "LYZ","S100A8","S100A9","FCGR3A","MS4A7","CD68","LST1","MRC1","CD163",
  # T/NK
  "CD3D","CD3E","CD2","IL7R","CCR7","CD4","CD8A","NKG7","GNLY","TRAC",
  # B/浆细胞
  "MS4A1","CD79A","CD79B","IGHG1","MZB1","SDC1",
  # 内皮
  "PECAM1","VWF","KDR","ESAM","ENG",
  # 成纤维/周细胞/平滑肌
  "COL1A1","COL1A2","DCN","LUM","PDGFRA","PDGFRB","RGS5","ACTA2","TAGLN"
)

p_dot <- DotPlot(seu, features = features_panel) + RotatedAxis()
ggsave("Figure/dotpanel_main_types.png", p_dot, width = 10, height = 6, dpi = 300)

p_feats <- FeaturePlot(
  seu, reduction = "umap",
  features = c("EPCAM","KRT8","CD3D","MS4A1","PECAM1","COL1A1"),
  ncol = 3, order = TRUE
)
ggsave("Figure/feature_main_six.png", p_feats, width = 9, height = 6, dpi = 300)

## 4) 手动注释（请根据 DotPlot/markers 观察结果，把映射改成你的实际对应）
cluster_ids <- levels(seu$seurat_clusters)

# 初始化为 "UNK"
manual_labels <- setNames(rep("UNK", length(cluster_ids)), nm = cluster_ids)

# ====== 示例映射（务必按你的结果修改！）======
# 比如你看到 0 类是上皮（EPCAM/KRT* 高），1 类是单核巨噬（LYZ/CD68 高）……
manual_labels[c("0")] <- "Epithelial" ##上皮细胞
manual_labels[c("1")] <- "Epithelial" ##上皮细胞
manual_labels[c("2")] <- "Mononuclear macrophages" #单核/巨噬细胞
manual_labels[c("3")] <- "Epithelial" ##上皮细胞
manual_labels[c("4")] <- "Pericyte / Vascular Smooth Muscle Cell" #周细胞/血管平滑肌细胞
manual_labels[c("5")] <- "Endothelial" #内皮细胞
manual_labels[c("6")] <- "Epithelial" #上皮细胞，混合少量髓系细胞
manual_labels[c("7")] <- "Plasma cells" #浆细胞
manual_labels[c("8")] <- "Plasma cells" #浆细胞
manual_labels[c("9")] <- "Cancer-Associated Fibroblast" #成纤维细胞
manual_labels[c("10")] <- "T lymphocytes" #T细胞

# 写入元数据并画图
suppressPackageStartupMessages(library(plyr))
seu$celltype_main <- plyr::mapvalues(seu$seurat_clusters, from = names(manual_labels), to = manual_labels)

p_umap_ct <- DimPlot(seu, reduction = "umap", group.by = "celltype_main", label = TRUE)
ggsave("Figure/umap_celltype_main.png", p_umap_ct, width = 7, height = 5, dpi = 300)

# 各大类的细胞数
ct_counts <- as.data.frame(table(seu$celltype_main))
colnames(ct_counts) <- c("celltype_main","n_cells")
write.csv(ct_counts, "Result/celltype_counts.csv", row.names = FALSE)



#######确定肿瘤细胞
if (!dir.exists("Figure")) dir.create("Figure", recursive = TRUE)
if (!dir.exists("Result")) dir.create("Result", recursive = TRUE)

library(copykat)
library(ggplot2)

# 2) 取“上皮大类”细胞（你上面已写入 celltype_main）
Idents(seu) <- "celltype_main"
epi_cells <- WhichCells(seu, idents = "Epithelial")
epi <- subset(seu, cells = epi_cells)

# 3) 准备原始表达矩阵（基因行、细胞列；行名需是人类基因符号）
exp.raw <- GetAssayData(epi, assay = "RNA", slot = "counts")
# 可选：去掉全零基因
exp.raw <- exp.raw[rowSums(exp.raw) > 0, ]

# 4) 跑 CopyKAT
#    id.type='S' 表示 Symbol（符号），人类数据；cell.line='no' 表示组织样本而非细胞系
gc()
set.seed(1234)
ck <- copykat::copykat(
  rawmat = as.matrix(exp.raw),
  id.type = "S",
  cell.line = "no",
  ngene.chr = 5,
  win.size = 25,
  KS.cut = 0.1,
  distance = "euclidean",
  n.cores = max(1, parallel::detectCores() - 1)
)

# 5) 取结果：pred.annovar 里有每个细胞的 "aneuploid"/"diploid"
ck.pred <- data.frame(
  cell = rownames(ck$prediction),
  copykat.pred = ck$prediction$copykat.pred,
  stringsAsFactors = FALSE
)

# 6) 合并回原 Seurat 元数据（先给全体细胞填默认 NA，再写入上皮的预测）
seu$copykat.pred <- NA
seu$copykat.pred[match(ck.pred$cell, colnames(seu))] <- ck.pred$copykat.pred

# 7) 定义“肿瘤 vs 正常”
#    规则：上皮细胞 & copykat.pred == "aneuploid" → Tumor；其余 → Non-tumor（含正常上皮/其它细胞类型）
seu$tumor_call <- "Non-tumor"
seu$tumor_call[seu$celltype_main == "Epithelial" & seu$copykat.pred == "aneuploid"] <- "Tumor (Epithelial CNV+)"

# 8) 可视化 & 导出
p_tumor_umap <- DimPlot(seu, reduction = "umap", group.by = "tumor_call", label = FALSE) +
  ggtitle("Tumor vs Non-tumor (by CopyKAT on Epithelial)")
ggsave("Figure/umap_tumor_vs_nontumor_copykat.png", p_tumor_umap, width = 7, height = 5, dpi = 300)

# 统计数量
tumor_stats <- as.data.frame(table(seu$tumor_call, seu$celltype_main))
colnames(tumor_stats) <- c("tumor_call","celltype_main","n_cells")
write.csv(tumor_stats, "Result/tumor_call_counts.csv", row.names = FALSE)

# 导出上皮细胞中 CNV 调用的明细
epi_dt <- data.frame(
  cell = epi_cells,
  copykat.pred = seu$copykat.pred[epi_cells],
  tumor_call = seu$tumor_call[epi_cells],
  cluster = seu$seurat_clusters[epi_cells],
  stringsAsFactors = FALSE
)
write.csv(epi_dt, "Result/epithelial_copykat_calls.csv", row.names = FALSE)

saveRDS(seu, "Result/seu_with_tumor_call.rds")


