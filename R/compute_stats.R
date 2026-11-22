#' @import SpatialFeatureExperiment
NULL

# Spatial statistics on a selected region. The neighbour graph is k-NN on
# cell centroids (available on every platform, unlike polygons), rebuilt on
# the selection so every cell has k neighbours inside the region.
# The @import above registers the coercion to SpatialFeatureExperiment.

.toSFEWithKNNNeighbors <- function(se, k = 6L) {
    sfe <- as(se, "SpatialFeatureExperiment")
    SpatialFeatureExperiment::colGraph(sfe, "knn") <-
        SpatialFeatureExperiment::findSpatialNeighbors(
            sfe, type = "spatialCoords", method = "knearneigh",
            k = as.integer(k),
            style = "W", zero.policy = TRUE)
    sfe
}

# Join counts need binary, symmetric weights: the test (Cliff and Ord 1981)
# assumes a symmetric matrix and counts joins, so the row-standardised k-NN
# graph used by the Moran statistics ("W", not symmetric) is not suitable;
# with it the variance can even be negative (?spdep::joincount.multi). Two
# cells are joined if either is among the other's k nearest neighbours.
.joincountWeights <- function(sfe) {
    nb <- spdep::make.sym.nb(SpatialFeatureExperiment::colGraph(sfe,
        "knn")$neighbours)
    spdep::nb2listw(nb, style = "B", zero.policy = TRUE)
}

.joincountDF <- function(sfe, group_col) {
    jc <- spdep::joincount.multi(factor(colData(sfe)[[group_col]]),
        .joincountWeights(sfe), zero.policy = TRUE)
    df <- as.data.frame(jc)
    colnames(df)[colnames(df) == "z-value"] <- "z_value"
    df$pair <- rownames(jc)
    df
}

.STATS_PLOTS <- c(
    "Spatial neighbour graph (k-NN)"  = "graph_plot",
    "Local Moran's I"                 = "localmoran_plot",
    "Moran scatter plot"              = "moran_plot",
    "Local cluster significance"      = "locclust_plot",
    "Join-count statistics"           = "joincount_plot",
    "Cell-type co-occurrence heatmap" = "joincount_heatmap_plot",
    "L-function (Lest)"               = "lest_plot",
    "L-cross function (Lcross)"       = "lcross_plot")

.computeSpatialStat <- function(se, plot_type, genes = character(0),
                                group_col = NA_character_, k = 6L,
                                alpha = 0.05) {
    has_groups <- !is.na(group_col) &&
        group_col %in% colnames(colData(se)) &&
        !all(is.na(colData(se)[[group_col]]))
    need_groups <- c("lest_plot", "lcross_plot", "joincount_plot",
        "joincount_heatmap_plot")
    # Unlabelled cells are left out only where the labels are used.
    if (has_groups && plot_type %in% need_groups)
        se <- se[, !is.na(colData(se)[[group_col]])]
    if (ncol(se) < max(5, k + 1)) {
        return(.notAvailablePlot(sprintf(
            "Too few cells (%d) in this selection for k = %d", ncol(se), k)))
    }
    if (plot_type %in% need_groups && !has_groups) {
        return(.notAvailablePlot(
            "No cell-type column chosen / available for this statistic"))
    }
    # L-functions work on point patterns; no neighbour graph involved.
    if (plot_type %in% c("lest_plot", "lcross_plot"))
        return(.lFunctionPlot(se, plot_type, group_col))

    sfe <- .toSFEWithKNNNeighbors(se, k)
    genes <- intersect(genes, rownames(sfe))
    gene_plots <- c("localmoran_plot", "moran_plot", "locclust_plot")
    if (plot_type %in% gene_plots && length(genes) == 0)
        return(.notAvailablePlot(
            "Choose at least one gene present in the data"))

    switch(plot_type,
        graph_plot = Voyager::plotColGraph(sfe, colGraphName = "knn",
            colGeometryName = 1) +
            ggplot2::theme_void(),
        localmoran_plot = .localMoranPlot(sfe, genes),
        moran_plot = .moranScatterPlot(sfe, genes[1]),
        locclust_plot = .localClusterPlot(sfe, genes[1], alpha, k),
        joincount_plot = .joincountBarPlot(sfe, group_col, k),
        joincount_heatmap_plot = .joincountHeatmap(sfe, group_col, k),
        .notAvailablePlot(paste("Unknown plot type:", plot_type)))
}

# Lest of the most frequent cell type, Lcross of the two most frequent ones
# (spatialFDA; types with >= 10 cells). The radius goes up to a quarter of
# the shorter side of the region, the usual limit for edge corrections.
.lFunctionPlot <- function(se, plot_type, group_col) {
    labs <- droplevels(factor(colData(se)[[group_col]]))
    counts <- sort(table(labs), decreasing = TRUE)
    top <- names(counts[counts >= 10])
    n_needed <- if (plot_type == "lcross_plot") 2 else 1
    if (length(top) < n_needed) {
        return(.notAvailablePlot(
            "Not enough cells per cell type in this selection"))
    }
    res <- spatialFDA::calcMetricPerFov(
        spe = se, fun = if (plot_type == "lest_plot") "Lest" else "Lcross",
        subsetby = "sample_id", marks = group_col,
        selection = top[seq_len(n_needed)],
        rSeq = .lRadii(SpatialExperiment::spatialCoords(se)),
        by = "sample_id", ncores = 1)
    spatialFDA::plotMetricPerFov(res, theo = TRUE, correction = "iso",
        x = "r", imageId = "sample_id")
}

.lRadii <- function(xy) {
    side <- min(diff(range(xy[, 1])), diff(range(xy[, 2])))
    seq(0, max(side, 1e-6) / 4, length.out = 100)
}

.runLocalMoran <- function(sfe, genes) {
    Voyager::runUnivariate(sfe, features = genes, type = "localmoran",
        zero.policy = TRUE, colGraphName = "knn",
        colGeometryName = 1)
}

.localMoranPlot <- function(sfe, genes) {
    sfe <- .runLocalMoran(sfe, genes)
    Voyager::plotLocalResult(sfe, name = "localmoran", features = genes,
        ncol = 2, colGeometryName = 1, divergent = TRUE,
        diverge_center = 0)
}

.moranScatterPlot <- function(sfe, gene) {
    sfe <- Voyager::runUnivariate(sfe, feature = gene, zero.policy = TRUE,
        type = "moran.plot", colGraphName = "knn",
        colGeometryName = 1)
    Voyager::moranPlot(sfe, feature = gene, graphName = "knn")
}

# Cells whose local Moran's I is significant (adjusted p < alpha), coloured
# by their cluster type (High-High, Low-Low, ...); the others are greyed.
.localClusterPlot <- function(sfe, gene, alpha, k) {
    sfe <- .runLocalMoran(sfe, gene)
    # localResults() is exported by SpatialFeatureExperiment, not Voyager.
    lm <- SpatialFeatureExperiment::localResults(sfe)[["localmoran"]][[1]]
    lm$locClust <- ifelse(lm[["-log10p_adj"]] > -log10(alpha),
        as.character(lm$mean), "non-significant")
    SpatialFeatureExperiment::localResults(sfe)[["localmoran"]][[1]] <- lm
    # The usual LISA colours: hot spots red, cold spots blue, the rest grey.
    cols <- c(`High-High` = "#D7191C", `Low-Low` = "#2C7BB6",
                `High-Low` = "#FDAE61", `Low-High` = "#ABD9E9",
                `non-significant` = "grey85")
    p <- Voyager::plotLocalResult(sfe, features = gene, name = "localmoran",
        attribute = "locClust", colGeometryName = 1)
    # Replace Voyager's colour scale (adding a second one prints a message).
    p$scales$scales <- Filter(function(s) !"colour" %in% s$aesthetics,
                                p$scales$scales)
    p + ggplot2::scale_colour_manual(values = cols, na.value = "grey85") +
        ggplot2::guides(colour = ggplot2::guide_legend(
            override.aes = list(size = 3))) +
        ggplot2::labs(caption = sprintf(
            "Significant: adjusted p < %g (local Moran's I, k = %d)",
            alpha, k))
}

.joincountBarPlot <- function(sfe, group_col, k) {
    df <- .joincountDF(sfe, group_col)
    df <- df[order(-abs(df$z_value)), ]
    df$pair <- factor(df$pair, levels = rev(df$pair))
    ggplot2::ggplot(df, ggplot2::aes(x = pair, y = z_value,
        fill = z_value > 0)) +
        ggplot2::geom_col() +
        ggplot2::coord_flip() +
        ggplot2::scale_fill_manual(
            values = c(`TRUE` = "#F7A11A", `FALSE` = "#4ECDC4"),
            guide = "none") +
        ggplot2::labs(x = NULL, y = "z-value", title = sprintf(
            "Join-count statistics (%s, k = %d)", group_col, k)) +
        ggplot2::theme_minimal()
}

.joincountHeatmap <- function(sfe, group_col, k) {
    df <- .joincountDF(sfe, group_col)
    lv <- levels(droplevels(factor(colData(sfe)[[group_col]])))
    z <- matrix(NA_real_, length(lv), length(lv), dimnames = list(lv, lv))
    parts <- strsplit(df$pair, ":")
    for (i in seq_along(parts)) {
        if (length(parts[[i]]) != 2) next          # skip the "Jtot" total row
        a <- parts[[i]][1]
        b <- parts[[i]][2]
        z[a, b] <- z[b, a] <- df$z_value[i]
    }
    zl <- max(abs(z), na.rm = TRUE)
    if (!is.finite(zl) || zl == 0) {
        return(.notAvailablePlot(
            "Not enough cell-type pairs in this selection"))
    }
    pal <- grDevices::colorRampPalette(c("#4ECDC4", "white", "#F7A11A"))
    pheatmap::pheatmap(z, cluster_rows = FALSE, cluster_cols = FALSE,
        color = pal(101), breaks = seq(-zl, zl, length.out = 102),
        # Values are unreadable beyond ~15 types; the colours remain.
        display_numbers = nrow(z) <= 15, number_format = "%.1f",
        na_col = "grey90",
        main = sprintf("Cell-type co-occurrence (join-count z, k = %d)", k),
        silent = TRUE)
}

# Moran's I scan over the most expressed genes (or metadata(se)$gene_ranks).
.computeSVGRanking <- function(se, top_n = 200L, k = 6L) {
    # An empty table carries the reason, shown in the panel.
    empty <- function(why) structure(data.frame(gene = character(0),
                                                moran_i = numeric(0)),
                                                reason = why)
    if (ncol(se) < max(5, k + 1))
        return(empty(sprintf(
            "Too few cells in this selection (%d) for k = %d.", ncol(se), k)))
    cand <- S4Vectors::metadata(se)$gene_ranks
    if (is.null(cand)) {
        a <- if ("counts" %in% SummarizedExperiment::assayNames(se))
            "counts" else "logcounts"
        tot <- Matrix::rowSums(SummarizedExperiment::assay(se, a))
        cand <- rownames(se)[order(tot, decreasing = TRUE)]
    }
    cand <- intersect(utils::head(cand, top_n), rownames(se))
    if (length(cand) == 0) return(empty("No genes to rank in this dataset."))
    sfe <- .toSFEWithKNNNeighbors(se, k)
    m <- Voyager::calculateMoransI(sfe, features = cand, colGraphName = "knn",
        colGeometryName = 1, zero.policy = TRUE)
    tbl <- data.frame(gene = rownames(m), moran_i = m$moran, row.names = NULL)
    tbl <- tbl[!is.na(tbl$moran_i), , drop = FALSE]
    if (nrow(tbl) == 0)
        return(empty(paste("Moran's I could not be computed for any gene",
            "(no variation in this selection?).")))
    tbl[order(-tbl$moran_i), , drop = FALSE]
}
