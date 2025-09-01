## ======================== Single-sample scRNA-seq pipeline ========================
options(stringsAsFactors = FALSE)
set.seed(1234)

## ---------------------------- 0) params & I/O -----------------------------------
params <- list(
  workdir     = "/mnt/DATA/home/Micoria/CA",
  h5_path     = "Data/1_filtered_feature_bc_matrix.h5",
  project     = "P1",
  min.cells   = 3,
  min.feats   = 200,
  max.feats   = 6000,
  max.mt      = 15,
  dims_use    = 1:30,
  res_main    = 0.4,
  ncores      = 16,        
  copykat     = list(ngene.chr = 5, win.size = 25, KS.cut = 0.1, dist = "euclidean"),
  caf_delta   = 0.15       
)

dirs <- list(
  fig   = file.path(params$workdir, "Figure"),
  res   = file.path(params$workdir, "Result")
)

## ------------------------------ 1) packages --------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(patchwork)
  library(ggplot2)
  library(DoubletFinder)          
  library(scDblFinder)
  library(SingleCellExperiment)
  library(copykat)
  library(GSVA)
  library(BiocParallel)
  library(plyr)
})

## ------------------------------ 2) helpers ---------------------------------------
logi <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "]", ..., "\n")

mk_dirs <- function() {
  for (d in unlist(dirs)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

guess_patterns <- function(genes) {
  list(
    mt   = if (any(grepl("^MT-", genes))) "^MT-" else if (any(grepl("^mt-", genes))) "^mt-" else if (any(grepl("^Mt-", genes))) "^Mt-" else NULL,
    ribo = if (sum(grepl("^RPL|^RPS", genes)) > 50) "^RPL|^RPS" else if (sum(grepl("^Rpl|^Rps", genes)) > 50) "^Rpl|^Rps" else NULL
  )
}

normalize_barcode <- function(x) {
  x <- sub("\\.1$", "-1", x)
  x <- sub("_1$",  "-1", x)
  x <- sub("^.*?([ACGT]{16}-\\d+)$", "\\1", x)
  x
}

ssgsea_scores <- function(mat, gsets, workers = 1) {
  present <- rownames(mat)
  gsets   <- lapply(gsets, function(gs) intersect(gs, present))
  lens    <- vapply(gsets, length, 1L)
  if (any(lens < 2)) {
    logi("Dropped gene sets (<2 genes):", paste(names(gsets)[lens < 2], collapse=", "))
    gsets <- gsets[lens >= 2]
  }
  stopifnot(length(gsets) >= 1)
  
  param <- GSVA::ssgseaParam(
    exprData  = mat, geneSets = gsets,
    minSize   = 1,   maxSize  = Inf,
    alpha     = 0.25, normalize = TRUE
  )
  bp <- if (.Platform$OS.type == "windows")
    SnowParam(workers = workers, progressbar = TRUE)
  else
    MulticoreParam(workers = workers, progressbar = TRUE)
  
  ssg <- GSVA::gsva(param, verbose = TRUE, BPPARAM = bp)  # rows=genesets, cols=cells
  
  z <- t(scale(t(ssg)))
  row_all_na <- function(v) all(is.na(v))
  for (i in seq_len(nrow(z))) {
    if (row_all_na(z[i, ])) {
      v <- as.numeric(ssg[i, ])
      if (sd(v) == 0) v <- v + rnorm(length(v), sd = 1e-8)
      z[i, ] <- scale(v)[,1]
    }
  }
  list(raw = ssg, z = z)
}

## ------------------------------ 3) I/O & QC -------------------------------------
setwd(params$workdir); mk_dirs()
logi("Reading 10x HDF5 ...")
mtx <- Read10X_h5(params$h5_path)

seu <- CreateSeuratObject(
  counts = mtx,
  project = params$project,
  min.cells = params$min.cells,
  min.features = params$min.feats
)

pat <- guess_patterns(rownames(seu))
if (!is.null(pat$mt))   seu[["percent.mt"]]   <- PercentageFeatureSet(seu, pattern = pat$mt)
if (!is.null(pat$ribo)) seu[["percent.ribo"]] <- PercentageFeatureSet(seu, pattern = pat$ribo)
logi("MT pattern:", ifelse(is.null(pat$mt),   "None", pat$mt),
     " | Ribo pattern:", ifelse(is.null(pat$ribo), "None", pat$ribo))

p_qc <- VlnPlot(seu, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)
ggsave(file.path(dirs$fig, "vio_QC.png"), p_qc, width = 8, height = 6, dpi = 300)

seu <- subset(seu, subset = nFeature_RNA >= params$min.feats &
                nFeature_RNA <= params$max.feats &
                (is.null(percent.mt) | percent.mt <= params$max.mt))

## ------------------------- 4) Normalize & dimension reduction --------------------
DefaultAssay(seu) <- "RNA"
logi("SCTransform (glmGamPoi=TRUE) ...")
seu <- SCTransform(seu, verbose = FALSE, method = "glmGamPoi")  # 更快、更稳
seu <- RunPCA(seu, verbose = FALSE)
seu <- FindNeighbors(seu, dims = params$dims_use)
seu <- FindClusters(seu, resolution = params$res_main)
seu <- RunUMAP(seu, dims = params$dims_use)
ggsave(file.path(dirs$fig, "umap_clusters_initial.png"),
       DimPlot(seu, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + NoLegend(),
       width = 7, height = 5, dpi = 300)

## ------------------------------- 5) Doublet removal ------------------------------
logi("Doublet removal via scDblFinder ...")
sce <- as.SingleCellExperiment(seu)
sce <- scDblFinder(sce, dbr = 0.07)
seu$scDblFinder.class <- colData(sce)$scDblFinder.class
logi("Doublet stats:"); print(table(seu$scDblFinder.class))

seu <- subset(seu, subset = scDblFinder.class == "singlet")
seu <- FindNeighbors(seu, dims = params$dims_use, verbose = FALSE)
seu <- FindClusters(seu, resolution = params$res_main,  verbose = FALSE)
seu <- RunUMAP(seu,   dims = params$dims_use,  verbose = FALSE)
ggsave(file.path(dirs$fig, "umap_clusters_after_dbl.png"),
       DimPlot(seu, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + NoLegend(),
       width = 7, height = 5, dpi = 300)

## --------------------------- 6) Marker discovery & panel -------------------------
DefaultAssay(seu) <- "SCT"
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers <- markers[order(markers$cluster, -markers$avg_log2FC), ]
write.csv(markers, file.path(dirs$res, "cluster_markers.csv"), row.names = FALSE)

panel <- c(
  "EPCAM","KRT8","KRT18","KRT19","MUC1",
  "LYZ","S100A8","S100A9","FCGR3A","MS4A7","CD68","LST1","MRC1","CD163",
  "CD3D","CD3E","CD2","IL7R","CCR7","CD4","CD8A","NKG7","GNLY","TRAC",
  "MS4A1","CD79A","CD79B","IGHG1","MZB1","SDC1",
  "PECAM1","VWF","KDR","ESAM","ENG",
  "COL1A1","COL1A2","DCN","LUM","PDGFRA","PDGFRB","RGS5","ACTA2","TAGLN"
)
ggsave(file.path(dirs$fig, "dotpanel_main_types.png"),
       DotPlot(seu, features = panel) + RotatedAxis(), width = 10, height = 6, dpi = 300)

## ------------------------------- 7) Manual annotation ----------------------------
cluster_ids <- levels(seu$seurat_clusters)
manual_labels <- setNames(rep("UNK", length(cluster_ids)), nm = cluster_ids)
manual_labels[c("0","1","3","6", "7")] <- "Epithelial"
manual_labels["2"] <- "TAMs"
manual_labels["4"] <- "VSMC"
manual_labels["5"] <- "Endothelial"
manual_labels["8"] <- "Plasma cells"
manual_labels["9"] <- "CAF"
manual_labels["10"] <- "T Cells"
seu$celltype_main <- plyr::mapvalues(seu$seurat_clusters, from = names(manual_labels), to = manual_labels)
ggsave(file.path(dirs$fig, "umap_celltype_main.png"),
       DimPlot(seu, reduction = "umap", group.by = "celltype_main", label = TRUE),
       width = 7, height = 5, dpi = 300)
write.csv(as.data.frame(table(seu$celltype_main)),
          file.path(dirs$res, "celltype_counts.csv"), row.names = FALSE)

## ------------------------------- 8) Tumor calling (CopyKAT) ----------------------
logi("CopyKAT on Epithelial ...")
Idents(seu) <- "celltype_main"
epi_cells <- WhichCells(seu, idents = "Epithelial")
epi <- subset(seu, cells = epi_cells)

exp.raw <- GetAssayData(epi, assay = "RNA", slot = "counts")
exp.raw <- exp.raw[rowSums(exp.raw) > 0, ]

ck <- copykat::copykat(
  rawmat  = as.matrix(exp.raw),
  id.type = "S", cell.line = "no",
  ngene.chr = params$copykat$ngene.chr,
  win.size  = params$copykat$win.size,
  KS.cut    = params$copykat$KS.cut,
  distance  = params$copykat$dist,
  n.cores   = max(1, params$ncores)
)

ck_cells_raw <- rownames(ck$prediction)
ck_pred_raw  <- ck$prediction$copykat.pred

ck_cells <- normalize_barcode(ck_cells_raw)
all_epi  <- normalize_barcode(colnames(epi))

pred_vec <- setNames(ck_pred_raw, ck_cells)
all_seu  <- normalize_barcode(colnames(seu))
vals <- pred_vec[all_seu]  
vals <- as.character(vals); names(vals) <- NULL
seu$copykat.pred <- vals

seu$tumor_call <- "Non-tumor"
seu$tumor_call[seu$celltype_main == "Epithelial" & seu$copykat.pred == "aneuploid"] <- "Tumor (Epithelial CNV+)"
logi("CopyKAT summary:"); print(table(seu$copykat.pred, useNA = "ifany"))
print(table(seu$tumor_call, useNA = "ifany"))

ggsave(file.path(dirs$fig, "umap_tumor_vs_nontumor_copykat.png"),
       DimPlot(seu, reduction = "umap", group.by = "tumor_call"),
       width = 7, height = 5, dpi = 300)
write.csv(as.data.frame(table(seu$tumor_call, seu$celltype_main)),
          file.path(dirs$res, "tumor_call_counts.csv"), row.names = FALSE)

write.csv(data.frame(
  cell    = colnames(epi),
  copykat = seu$copykat.pred[colnames(epi)],
  tumor   = seu$tumor_call[colnames(epi)],
  cluster = seu$seurat_clusters[colnames(epi)]
), file.path(dirs$res, "epithelial_copykat_calls.csv"), row.names = FALSE)

## ------------------------------- 9) CAF subtyping (ssGSEA) -----------------------
logi("CAF subtyping (ssGSEA) ...")
Idents(seu) <- "celltype_main"
caf_cells <- WhichCells(seu, idents = "CAF")
stopifnot(length(caf_cells) > 0)
caf <- subset(seu, cells = caf_cells)
assays_available <- tryCatch(Seurat::Assays(caf), error = function(e) names(caf@assays))
if ("SCT" %in% assays_available) {
  DefaultAssay(caf) <- "SCT"
} else if ("RNA" %in% assays_available) {
  DefaultAssay(caf) <- "RNA"
} else {
  DefaultAssay(caf) <- assays_available[1]
}
if (DefaultAssay(caf) == "RNA" && nrow(GetAssayData(caf, slot = "data")) == 0) {
  caf <- NormalizeData(caf, verbose = FALSE)
}
emat <- as.matrix(GetAssayData(caf, slot = "data"))
stopifnot(nrow(emat) > 0, ncol(emat) > 0)

gs.caf <- list(
  mCAF      = c("ACTA2","TAGLN","MYH11","CNN1","MYL9"),
  iCAF      = c("CXCL12","IL6","CCL2","CXCL14","LIF","PTGS2","IL11"),
  `ECM-CAF` = c("COL11A1","COL10A1","POSTN","THBS2","FN1","COL3A1","COL1A1","COL1A2"),
  apCAF     = c("HLA-DRA","HLA-DRB1","CD74","CIITA","HLA-DPA1","HLA-DPB1")
)

ssg <- ssgsea_scores(emat, gs.caf, workers = min(4, params$ncores))

for (gs in rownames(ssg$raw)) {
  caf[[paste0(gs, "_ssgsea")]]   <- as.numeric(ssg$raw[gs, ])
  caf[[paste0(gs, "_ssgsea_z")]] <- as.numeric(ssg$z[gs, ])
}

get_top     <- function(v) if (all(is.na(v))) NA_integer_ else which.max(replace(v, is.na(v), -Inf))
get_gap_col <- function(v) { vv <- sort(replace(v, is.na(v), -Inf), decreasing = TRUE); if (length(vv) < 2) NA_real_ else vv[1] - vv[2] }

zmat   <- ssg$z                             # rows=genesets, cols=cells
top_i  <- apply(zmat, 2, get_top)
top_lb <- ifelse(is.na(top_i), NA, rownames(zmat)[top_i])
gap    <- apply(zmat, 2, get_gap_col)

caf$CAF_subtype <- ifelse(is.na(top_lb), NA, ifelse(gap >= params$caf_delta, top_lb, "CAF-ambiguous"))
logi("CAF subtype (CAF subset):"); print(table(caf$CAF_subtype, useNA = "ifany"))

caf_idx <- match(normalize_barcode(colnames(caf)), normalize_barcode(colnames(seu)))
seu$CAF_subtype <- NA_character_
seu$CAF_subtype[caf_idx] <- as.character(caf$CAF_subtype)

logi("CAF subtype (in SEURAT, CAF only):"); print(table(seu$CAF_subtype[caf_cells], useNA = "ifany"))

cells_use <- intersect(caf_cells, colnames(seu)[!is.na(seu$CAF_subtype)])
if (length(cells_use) > 0) {
  caf_pl <- subset(seu, cells = cells_use)
  caf_pl$CAF_subtype <- factor(caf_pl$CAF_subtype, levels = c("mCAF","iCAF","ECM-CAF","apCAF","CAF-ambiguous"))
  ggsave(file.path(dirs$fig, "umap_CAF_subtypes_ssgsea.png"),
         DimPlot(caf_pl, reduction = "umap", group.by = "CAF_subtype") +
           ggtitle(sprintf("CAF subtypes (ssGSEA, assay=%s)", DefaultAssay(caf))),
         width = 7, height = 5, dpi = 300)
  
  panel_caf <- intersect(unique(c(
    "ACTA2","TAGLN","MYH11","CNN1","MYL9",
    "CXCL12","IL6","CCL2","CXCL14","LIF","PTGS2","IL11",
    "COL11A1","COL10A1","POSTN","THBS2","FN1","COL3A1","COL1A1","COL1A2",
    "HLA-DRA","HLA-DRB1","CD74","CIITA","HLA-DPA1","HLA-DPB1"
  )), rownames(caf_pl))
  if (length(panel_caf) > 0) {
    ggsave(file.path(dirs$fig, "dotpanel_CAF_subtypes_ssgsea.png"),
           DotPlot(caf_pl, features = panel_caf, group.by = "CAF_subtype") + RotatedAxis(),
           width = 10, height = 6, dpi = 300)
  }
  write.csv(as.data.frame(table(caf_pl$CAF_subtype, useNA = "ifany")),
            file.path(dirs$res, "CAF_subtype_counts_ssgsea.csv"), row.names = FALSE)
}

## ------------------------------- 10) Save ----------------------------------------
saveRDS(seu, file.path(dirs$res, "seu_with_tumor_and_CAF.rds"))
logi("All done.")
