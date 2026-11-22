# ============================================================
# Xenium Human Colon workflow -- adapted from the OSTA book chapter:
#   https://bioconductor.org/books/release/OSTA/pages/img-workflow-xenium.html
#   (OSTA book, CC BY 4.0)
#
# The code follows the book. One difference: the book uses scrapper's
# *.se wrappers (quickRnaQc.se, runPca.se, ...), which older scrapper
# versions lack; the steps below call the functions they wrap
# (computeRnaQcMetrics, normalizeCounts, runPca, buildSnnGraph,
# clusterGraph, scoreMarkers), so they run with both.
# ============================================================

library(BiocParallel)
library(DropletUtils)
library(ggplot2)
library(OSTA.data)
library(patchwork)
library(scater)
library(scrapper)
library(SingleR)
library(SpatialExperiment)
library(SpatialExperimentIO)

bp <- MulticoreParam(th <- 4)
set.seed(112358)

# ------------------------------------------------------------
# Load Xenium dataset
# ------------------------------------------------------------
id <- "Xenium_HumanColon_Oliveira"
pa <- OSTA.data_load(id, mol = FALSE)
dir.create(td <- tempfile())
unzip(pa, exdir = td)
# Book uses readXeniumSXE(td, addTx=FALSE) -- that parameter doesn't exist
# in this installed SpatialExperimentIO version (another version-gap vs.
# the book, same story as scrapper's .se wrappers). addParquetPaths=FALSE
# achieves the same thing: skip looking for transcripts.parquet, which
# isn't present since we downloaded with mol=FALSE.
(spe <- readXeniumSXE(td, addParquetPaths = FALSE))

# ------------------------------------------------------------
# Helper for spatial scatter plots (verbatim from the book)
# ------------------------------------------------------------
.plt_xy <- \(spe, col) {
  df <- data.frame(colData(spe), spatialCoords(spe))
  aes <- if (is.numeric(df[[col]])) {
    theme(
      legend.key.height = unit(1, "lines"),
      legend.key.width  = unit(0.5, "lines"))
  } else {
    list(
      theme(legend.key.size = unit(0, "lines")),
      guides(col = guide_legend(override.aes = list(alpha = 1, size = 2))))
  }
  ggplot(df, aes(x_centroid, y_centroid, col = .data[[col]])) +
    coord_equal() + theme_void() + aes +
    geom_point(stroke = 0, size = 1/3)
}

# ------------------------------------------------------------
# Quality control
#   Book:   spe <- quickRnaQc.se(spe, subsets=list())
#   Here:   compute -> suggest thresholds -> filter, by hand
# ------------------------------------------------------------
qc  <- computeRnaQcMetrics(counts(spe), subsets = list())
thr <- suggestRnaQcThresholds(qc)
spe$sum      <- qc$sum
spe$detected <- qc$detected
spe$keep     <- filterRnaQcMetrics(thr, qc)
round(100 * prop.table(table(keep = !spe$keep)), 2)

plot_dir <- "xenium_plots"
dir.create(plot_dir, showWarnings = FALSE)

p_qc <- .plt_xy(spe, "sum") +
  scale_color_viridis_c(
    "# counts", trans = "log1p",
    breaks = range(spe$sum), labels = c("low", "high")) +
  .plt_xy(spe, "detected") +
  scale_color_viridis_c(
    "# features", breaks = range(spe$detected), labels = c("low", "high")) +
  .plt_xy(spe[, order(!spe$keep)], "keep") +
  scale_color_manual(
    "high-quality", labels = c("no", "yes"), values = c("purple", "lavender"))
ggsave(file.path(plot_dir, "01_qc.png"), p_qc, width = 12, height = 4, dpi = 150, create.dir = TRUE)

ncol(spe <- spe[, spe$keep])

# ------------------------------------------------------------
# Define and crop to a spatial region (verbatim from the book)
# ------------------------------------------------------------
box <- list(xmin = 2e3, xmax = 5e3, ymin = 1e3, ymax = 4e3)

xy <- spatialCoords(spe)
i <-
  xy[, 1] > box$xmin &
  xy[, 1] < box$xmax &
  xy[, 2] > box$ymin &
  xy[, 2] < box$ymax
ncol(sub <- spe[, i])

df <- data.frame(xy, i)
p <- ggplot(df, aes(x_centroid, y_centroid)) +
  coord_equal() + theme_void() +
  theme(legend.position = "none")
p_crop <- p + geom_point(aes(col = i), stroke = 0, size = 0.1) |
  p + geom_point(data = df[i, ], stroke = 0, size = 0.2)
ggsave(file.path(plot_dir, "02_crop.png"), p_crop, width = 10, height = 5, dpi = 150, create.dir = TRUE)

# ------------------------------------------------------------
# Normalization and PCA
#   Book: sub <- normalizeRnaCounts.se(sub, size.factors=sfs)
#         sub <- runPca.se(sub, features=rownames(sub))
# ------------------------------------------------------------
sfs <- (sub$cell_area) / median(sub$cell_area)
logcounts(sub) <- normalizeCounts(counts(sub), size.factors = sfs)

pca <- runPca(as.matrix(logcounts(sub)[rownames(sub), ]), number = 25)
reducedDim(sub, "PCA") <- t(pca$components)

# ------------------------------------------------------------
# Visualize gene expression (verbatim from the book)
# ------------------------------------------------------------
gs <- c("PIGR", "IGHG3", "CEACAM6")
gs <- intersect(gs, rownames(sub))  # panel may not include all 3
es <- scale(t(as.matrix(logcounts(sub))))
es <- es[, gs, drop = FALSE]
colData(sub) <- cbind(colData(sub), es)
ps <- lapply(gs, \(.) .plt_xy(sub, .) + ggtitle(.))
p_gs <- wrap_plots(ps, nrow = 1) &
  scale_color_gradientn(NULL,
    labels = c("low", "high"),
    colors = rev(hcl.colors(9, "PuRd")),
    limits = rng, breaks = rng <- range(es)) &
  theme(plot.title = element_text(hjust = 0.5))
ggsave(file.path(plot_dir, "03_genes.png"), p_gs, width = 4 * length(gs), height = 4, dpi = 150, create.dir = TRUE)

# ------------------------------------------------------------
# Unsupervised clustering
#   Book: sub <- clusterGraph.se(sub, num.threads=th, resolution=0.5,
#                 method="leiden", output.name="Leiden",
#                 more.build.args=list(weight.scheme="jaccard"))
# ------------------------------------------------------------
snn <- buildSnnGraph(pca$components, weight.scheme = "jaccard", num.threads = th)
cl  <- clusterGraph(snn, method = "leiden", leiden.resolution = 0.5)
sub$Leiden <- cl$membership
table(sub$Leiden)

# ------------------------------------------------------------
# Supervised annotation: reference scRNA-seq + SingleR
# (unchanged from the book -- doesn't use scrapper's .se wrappers)
#
# Checkpointed to a local (non-OneDrive) path: this step is the most
# expensive in the whole script (SingleR against a ~280k-cell reference),
# and OneDrive has already caused two unrelated stalls/flakes today. If a
# later step (e.g. plotting) fails, re-running should not have to pay for
# this again.
# ------------------------------------------------------------
checkpoint <- "/tmp/xenium_osta_sub_annotated.rds"
if (file.exists(checkpoint)) {
  sub <- readRDS(checkpoint)
  cat("Loaded annotated `sub` from checkpoint:", checkpoint, "\n")
} else {
  id <- "Chromium_HumanColon_Oliveira"
  pa <- OSTA.data_load(id)
  dir.create(td <- tempfile())
  unzip(pa, exdir = td)

  sce <- read10xCounts(list.files(td, "h5$", full.names = TRUE))
  cd  <- read.csv(list.files(td, "cell_meta", full.names = TRUE))
  colData(sce) <- cbind(colData(sce), cd[, -1])
  table(sce$Level1)
  ncol(sce)

  sce <- sce[, sce$QCFilter == "Keep"]
  sce <- sce[, grepl("P2", sce$Patient)]
  assay(sce) <- as(assay(sce), "dgCMatrix")
  logcounts(sce) <- normalizeCounts(assay(sce), size.factors = rep(1, ncol(sce)))
  sce <- sce[rowData(sce)$ID %in% rowData(sub)$ID, ]
  rownames(sce) <- rowData(sce)$Symbol

  res <- SingleR(
    test = sub, ref = sce,
    labels = sce$Level2,
    de.method = "wilcox",
    aggr.ref = TRUE, BPPARAM = bp)

  sub$Level2 <- factor(res$pruned.labels)

  idx <- match(sub$Level2, sce$Level2)
  table(sub$Level1 <- factor(sce$Level1[idx]))

  lab <- list(
    tum = c("Tumor"),
    epi = c("Intestinal Epithelial"),
    imm = c("B cells", "T cells", "Myeloid"),
    str = c("Endothelial", "Fibroblast", "Smooth Muscle"))
  idx <- match(sub$Level1, unlist(lab))
  lab <- rep.int(names(lab), sapply(lab, length))
  table(sub$Level0 <- factor(lab[idx]))

  saveRDS(sub, checkpoint)
  cat("Saved annotated `sub` checkpoint to", checkpoint, "\n")
}

round(100 * prop.table(table(sub$Level1, sub$Leiden), 2), 1)
c(SingleR = nlevels(sub$Level1), Leiden = nlevels(sub$Leiden))

lapply(c("Leiden", "Level0", "Level1"), \(.) {
  pal <- if (. == "Level0") {
    c("gold", "cyan", "magenta", "black")
  } else {
    hcl.colors(nlevels(sub[[.]]), "Spectral")
  }
  .plt_xy(sub, .) + scale_color_manual(values = pal)
}) |> wrap_plots() -> p_annot
ggsave(file.path(plot_dir, "04_annotation.png"), p_annot, width = 14, height = 5, dpi = 150, create.dir = TRUE)

# ------------------------------------------------------------
# Marker genes
#   Book: de <- scoreMarkers.se(sub[, ok], sub$Level1[ok])
#         gs <- lapply(de, \(df) head(rownames(df), 5))
#   Here: scoreMarkers() + sort each group's Cohen's d table by
#   min.rank ourselves (verified this is what min.rank is for: a
#   single overall ranking combining all pairwise comparisons).
# ------------------------------------------------------------
ok <- !is.na(sub$Level1)
de_raw <- scoreMarkers(as.matrix(logcounts(sub))[, ok], groups = sub$Level1[ok])
de <- lapply(de_raw$cohens.d, \(df) df[order(df$min.rank), ])

gs <- lapply(de, \(df) head(rownames(df), 5))
gs <- unique(unlist(gs))

p_heat <- plotGroupedHeatmap(sub,
  features = gs, group = "Level1",
  scale = TRUE, center = TRUE, fontsize = 6)
ggsave(file.path(plot_dir, "05_markers.png"), p_heat, width = 8, height = 6, dpi = 150, create.dir = TRUE)

cat("\nDone. Plots saved to", plot_dir, "\n")
