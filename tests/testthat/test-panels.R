test_that("panel constructors set their defaults", {
    expect_identical(SpatialMapPanel()@PendingIDs, "")
    p <- SpatialPlotPanel()
    expect_identical(p@ColorBy, "")
    expect_identical(p@Selection, "")
    expect_identical(SpatialQCPanel()@UseFences, TRUE)
    s <- SpatialStatsPanel()
    expect_identical(s@PlotType, "graph_plot")
    expect_identical(s@K, 6L)
    v <- SVGPanel(K = 12)
    expect_identical(v@K, 12L)
    expect_identical(v@TopN, 200L)
    n <- SpatialNichePanel(K = 30)
    expect_identical(n@K, 30L)
})

test_that("every panel reports a name and a valid hex colour", {
    panels <- list(SpatialMapPanel(), SpatialPlotPanel(), FOVZoomPanel(),
                   SpatialQCPanel(),
        FOVStatsPanel(), SpatialStatsPanel(), SVGPanel(), SpatialNichePanel())
    for (p in panels) {
        expect_type(.fullName(p), "character")
        expect_match(.panelColor(p), "^#[0-9A-Fa-f]{6}$")
    }
})

test_that("all selection-following panels share the Selection slot", {
    for (p in list(SpatialPlotPanel(), FOVZoomPanel(), SpatialQCPanel(),
                   FOVStatsPanel(),
        SpatialStatsPanel(), SVGPanel(), SpatialNichePanel()))
        expect_true(is(p, "SpatialPanelBase"))
    expect_false(is(SpatialMapPanel(), "SpatialPanelBase"))
})

test_that("default layouts adapt to the platform", {
    cos_qc <- .defaultPanels(.makeTestSPE(), "qc")
    expect_true(any(vapply(cos_qc, is, logical(1), "FOVZoomPanel")))
    xen_qc <- .defaultPanels(.makeTestXeniumSPE(), "qc")
    expect_false(any(vapply(xen_qc, is, logical(1), "FOVZoomPanel")))
    an <- .defaultPanels(.makeTestVisiumHDSPE(), "analysis")
    expect_true(any(vapply(an, is, logical(1), "SpatialNichePanel")))
    expect_identical(an[[2]]@ColorBy, "Level1")
})

test_that("plot panels render on synthetic data without the app", {
    skip_if_not_installed("SpaceTrooper")
    spe <- .makeTestVisiumHDSPE()
    expect_s3_class(.generateOutput(SpatialPlotPanel(ColorBy = "sum"), spe),
                    "ggplot")
    expect_s3_class(.generateOutput(SpatialPlotPanel(ColorBy = "GENE1"), spe),
                    "ggplot")
    p <- .generateOutput(SpatialPlotPanel(UsePolygons = TRUE), spe)
    expect_s3_class(p, "ggplot")    # graceful message, no polygons loaded
    # not a tile platform: a message, not an error
    expect_s3_class(.generateOutput(FOVZoomPanel(), spe), "ggplot")
})

test_that("stats degrade gracefully on too few cells and honour k", {
    spe <- .makeTestSPE(n_cells = 3)
    expect_s3_class(.computeSpatialStat(spe, "graph_plot", k = 6), "ggplot")
    tbl <- .computeSVGRanking(spe)
    expect_equal(nrow(tbl), 0)
    expect_named(tbl, c("gene", "moran_i"))
    expect_s3_class(.computeSpatialStat(.makeTestSPE(), "joincount_plot",
                                        group_col = NA), "ggplot")
})

# ── Integration on real data, when present locally ──
test_that("stats and SVG ranking run on the real CosMx data", {
    path <- .realCosmxPath()
    skip_if(is.na(path), "real CosMx dataset not present on this machine")
    spe <- setPlatform(readRDS(path), "cosmx")
    sub <- .subsetBySelection(spe, "72,73,74")
    expect_s3_class(.computeSpatialStat(sub, "graph_plot", k = 10), "ggplot")
    tbl <- .computeSVGRanking(sub, top_n = 20, k = 10)
    expect_gt(nrow(tbl), 0)
    expect_true(all(diff(tbl$moran_i) <= 0))
})

test_that("Xenium test bundle: loader, k-NN graph and known top SVG", {
    path <- .xeniumTestDataPath()
    skip_if(is.na(path), "bundled Xenium test data not found")
    skip_if_not_installed("SpatialExperimentIO")
    spe <- loadXenium(path)
    expect_identical(getPlatform(spe), "xenium")
    expect_equal(ncol(spe), 36)
    expect_s3_class(.computeSpatialStat(spe, "graph_plot"), "ggplot")
    tbl <- .computeSVGRanking(spe, top_n = 50)
    # verified on this bundle with k = 6
    expect_identical(tbl$gene[1], "Defa5")
})

test_that("genes with a hyphen in their name can be plotted", {
    skip_if_not_installed("SpaceTrooper")
    spe <- .makeTestVisiumHDSPE()
    rownames(spe)[1] <- "MT-ATP6"
    p <- .generateOutput(SpatialPlotPanel(ColorBy = "MT-ATP6"), spe)
    expect_s3_class(p, "ggplot")
    expect_no_error(ggplot2::ggplot_build(p))
})

test_that("SVG tables are kept per session, not shared", {
    a <- new.env()
    b <- new.env()
    assign("p", data.frame(gene = "A"), envir = .svgTables(a))
    assign("p", data.frame(gene = "B"), envir = .svgTables(b))
    expect_identical(get("p", envir = .svgTables(a))$gene, "A")
    expect_identical(get("p", envir = .svgTables(b))$gene, "B")
})

test_that("the launcher offers the spatial panels in 'Organize panels'", {
    expect_s3_class(spiSEE(.makeTestVisiumHDSPE(), mode = "qc",
                                cache = FALSE), "shiny.appobj")
    expect_false(any(vapply(.extraPanels(.makeTestVisiumHDSPE()), is,
                            logical(1), "FOVZoomPanel")))
    expect_true(any(vapply(.extraPanels(.makeTestSPE()), is, logical(1),
                           "FOVZoomPanel")))
})

test_that("join counts use binary symmetric weights", {
    sfe <- .toSFEWithKNNNeighbors(.makeTestSPE(n_cells = 120), k = 6L)
    lw <- .joincountWeights(sfe)
    expect_identical(lw$style, "B")
    m <- spdep::listw2mat(lw)
    expect_true(isSymmetric(unname(m)))
    expect_true(all(m %in% c(0, 1)))
    df <- .joincountDF(sfe, "SingleR_label")
    expect_true(all(df$Variance[!is.na(df$Variance)] >= 0))
})

test_that("local cluster significance level is a panel parameter", {
    expect_identical(SpatialStatsPanel()@Alpha, 0.05)
    expect_identical(SpatialStatsPanel(Alpha = 0.01)@Alpha, 0.01)
    expect_error(SpatialStatsPanel(Alpha = 1.5), "Alpha")
    expect_true("alpha" %in% names(formals(.computeSpatialStat)))
})

test_that("QC histogram: readable names and flagged-cell count", {
    spe <- .makeTestSPE(n_cells = 200)
    colData(spe)$total_outlier_mc <- scuttle::isOutlier(colData(spe)$total,
        nmads = 1)
    n <- sum(colData(spe)$total_outlier_mc)
    p <- .generateOutput(SpatialQCPanel(Metric = "total", UseFences = TRUE),
                         spe)
    expect_s3_class(p, "ggplot")
    expect_identical(p$labels$colour, "Outlier fences (medcouple)")
    expect_match(p$labels$subtitle, sprintf("^%d of 200 cells flagged", n))
    expect_identical(.prettyMetric("transcript_counts"), "Transcript counts")
    mc <- structure(c("NO", "HIGH", "LOW", "NO"),
                    class = c("outlier.filter", "logical"))
    expect_identical(.isFlagged(mc), c(FALSE, TRUE, TRUE, FALSE))
    expect_identical(.isFlagged(c(TRUE, FALSE)), c(TRUE, FALSE))
    expect_identical(.flaggedText(c(TRUE, FALSE, NA, FALSE)),
                     "1 of 4 cells flagged as outliers (25.0%)")
})

test_that("placeholders show their message in the panel", {
    p <- .generateOutput(FOVZoomPanel(), .makeTestSPE())
    expect_identical(p$labels$title,
                     "Select FOVs on the map, then press Validate.")
    expect_true(any(vapply(p$layers, function(l) inherits(l$geom, "GeomText"),
                           logical(1))))
    f <- tempfile(fileext = ".png")
    grDevices::png(f, 500, 400)
    print(p)
    grDevices::dev.off()
    expect_gt(file.size(f), 0)
})

test_that("FOV zoom explains a dataset without cell polygons", {
    spe <- .makeTestSPE()
    spe$polygons <- NULL
    fov <- as.character(spe$fov[1])
    p <- .generateOutput(FOVZoomPanel(Selection = fov), spe)
    expect_match(p$labels$title, "has none")
})

test_that("panel outputs turn on Shiny's busy indicators", {
    for (x in list(SpatialPlotPanel(), SpatialMapPanel(), SVGPanel(),
                   SpatialNichePanel())) {
        x@PanelId <- 1L
        html <- htmltools::renderTags(.defineOutput(x))$html
        expect_match(html, "shinyBusySpinners")
        expect_match(html, .getEncodedName(x), fixed = TRUE)
    }
})

test_that("scale bar is turned off only when the technology is unknown", {
    f_new <- function(spe, scaleBar = TRUE) NULL      # SpaceTrooper >= 1.3
    f_old <- function(spe) NULL                       # SpaceTrooper 1.0
    spe <- .makeTestSPE()
    S4Vectors::metadata(spe)$technology <- NULL
    expect_identical(.scaleBarArgs(f_new, spe), list(scaleBar = FALSE))
    expect_identical(.scaleBarArgs(f_old, spe), list())
    S4Vectors::metadata(spe)$technology <- "10X_Xenium"
    expect_identical(.scaleBarArgs(f_new, spe), list())
})

test_that("an empty SVG table says why", {
    spe <- .makeTestSPE()
    expect_match(attr(.generateOutput(SVGPanel(), spe), "reason"),
                 "Press Compute")
    few <- .computeSVGRanking(spe[, 1:4], k = 6L)
    expect_equal(nrow(few), 0L)
    expect_match(attr(few, "reason"), "Too few cells")
})

test_that("an initial selection is shown on the map", {
    spe <- .makeTestSPE()
    panels <- .defaultPanels(spe, "qc", selection = "1,2")
    expect_identical(panels[[1]]@CommittedIDs, "1,2")
    expect_match(.MAP_JS, "knownOnly(data.committed)", fixed = TRUE)
})

test_that("gene menus put the most expressed genes first", {
    spe <- .makeTestSPE()
    S4Vectors::metadata(spe)$gene_ranks <- NULL
    top <- rownames(spe)[which.max(Matrix::rowMeans(logcounts(spe)))]
    S4Vectors::metadata(spe)$spiSEE_gene_order <- .geneOrder(spe)
    expect_identical(.geneChoices(spe)[1], top)
})

test_that("plotSpatialCells draws polygons, points, genes and panels", {
    spe <- .makeTestSPE(n_cells = 60)
    spe$grp <- rep(c("a", "b", "c"), 20)
    p <- plotSpatialCells(spe, "grp")
    expect_s3_class(p, "ggplot")
    expect_true(inherits(p$layers[[1]]$geom, "GeomPoint"))   # no polygons
    g <- rownames(spe)[1]
    p2 <- plotSpatialCells(spe, c(g, rownames(spe)[2]), ncol = 2)
    expect_s3_class(p2$facet, "FacetWrap")
    expect_error(plotSpatialCells(spe, c("grp", g)), "mixes")
    expect_error(plotSpatialCells(spe, "nope"), "neither")
    sq <- function(x, y) sf::st_polygon(list(rbind(c(x, y), c(x + 1, y),
        c(x + 1, y + 1), c(x, y + 1), c(x, y))))
    xy <- SpatialExperiment::spatialCoords(spe)
    spe$polygons <- sf::st_sf(cell_id = seq_len(ncol(spe)),
        global = sf::st_sfc(lapply(seq_len(ncol(spe)),
                                   function(i) sq(xy[i, 1], xy[i, 2]))))
    p3 <- plotSpatialCells(spe, "grp", borders = TRUE)
    expect_true(inherits(p3$layers[[1]]$geom, "GeomSf"))
})

test_that("the selection given at launch is checked", {
    spe <- .makeTestSPE()
    expect_identical(.checkSelection(" 1, 2 ,", spe), "1,2")
    expect_identical(.checkSelection(c("2", "1"), spe), "2,1")
    expect_warning(sel <- .checkSelection("1,99", spe), "99")
    expect_identical(sel, "1")
    expect_identical(.checkSelection("", spe), "")
})

test_that("one sample per app; cell names made unique", {
    a <- .makeTestXeniumSPE()
    b <- .makeTestXeniumSPE(seed = 2)
    b$sample_id <- "s2"
    both <- cbind(a, b)
    expect_error(spiSEE(both, cache = FALSE), "sample =")
    expect_identical(ncol(.oneSample(both, "s2")), ncol(b))
    expect_error(.oneSample(both, "nope"), "one of")
    expect_message(app <- spiSEE(both, sample = "s2", cache = FALSE),
                    NA)
    dup <- a
    colnames(dup)[2] <- colnames(dup)[1]
    expect_message(spiSEE(dup, cache = FALSE), "not unique")
})

test_that("the app starts on CosMx-like data without a FOV table", {
    app <- spiSEE(.makeTestSPE(), mode = "qc", cache = FALSE)
    expect_s3_class(app, "shiny.appobj")
    shiny::testServer(app, {
        session$flushReact()
        expect_true(TRUE)
    })
})

test_that("a panel added while the app runs takes the map selection", {
    spe <- .makeTestSPE()
    shiny::testServer(function(input, output, session) {
        pObjects <- new.env()
        rObjects <- shiny::reactiveValues(modified = list())
        map <- SpatialMapPanel(CommittedIDs = "1,2", PanelId = 1L)
        new_panel <- SpatialQCPanel(PanelId = 1L)
        pObjects$memory <- list(SpatialMapPanel1 = map,
                                SpatialQCPanel1 = new_panel)
        .createObservers(new_panel, spe, input, session, pObjects, rObjects)
    }, {
        expect_identical(pObjects$memory$SpatialQCPanel1@Selection, "1,2")
    })
})

test_that("text and logical label columns are drawn as categories", {
    spe <- .makeTestSPE()
    spe$lab_chr <- as.character(spe$SingleR_label)
    spe$flag <- spe$total > 50
    expect_true(is.factor(.forceCategorical(spe, "lab_chr")$lab_chr))
    expect_true(is.factor(.forceCategorical(spe, "flag")$flag))
    expect_true(is.numeric(.forceCategorical(spe, "total")$total))
    p <- .generateOutput(SpatialPlotPanel(ColorBy = "lab_chr",
                                            UsePolygons = FALSE), spe)
    expect_s3_class(p, "ggplot")
})

test_that("L-function radii scale with the region", {
    r <- .lRadii(cbind(c(0, 400), c(0, 1000)))
    expect_equal(max(r), 100)
    expect_length(r, 100)
})

test_that("gene statistics keep the cells without a label", {
    spe <- .makeTestSPE(n_cells = 60)
    spe$SingleR_label[1:10] <- NA
    seen <- new.env()
    local_mocked_bindings(.toSFEWithKNNNeighbors = function(se, k) {
        seen$n <- ncol(se)
        stop("stop here")
    })
    try(.computeSpatialStat(spe, "graph_plot", group_col = "SingleR_label"),
        silent = TRUE)
    expect_identical(seen$n, 60L)
    try(.computeSpatialStat(spe, "joincount_plot",
                            group_col = "SingleR_label"), silent = TRUE)
    expect_identical(seen$n, 50L)
})

test_that("loaders read the example data of the reader packages", {
    dir <- system.file("extdata", "CosMx_DBKero_Tiny", package = "SpaceTrooper")
    skip_if(!nzchar(dir))
    f <- withr::local_tempfile(fileext = ".rds")
    saveRDS(SpaceTrooper::readCosmxSPE(dir, sampleName = "DBKero"), f)
    spe <- suppressMessages(loadCosmx(f))
    expect_identical(getPlatform(spe), "cosmx")
    skip_if_not_installed("SpatialExperimentIO")
    d <- system.file("extdata", "Xenium_small",
                        package = "SpatialExperimentIO")
    x <- suppressMessages(loadXenium(d))
    expect_identical(getPlatform(x), "xenium")
    expect_true("logcounts" %in% SummarizedExperiment::assayNames(x))
    skip_if_not_installed("VisiumIO")
    skip_if_not_installed("BiocIO")
    v <- suppressMessages(suppressWarnings(loadVisiumHD(
        system.file("extdata", "segmented_outputs", package = "VisiumIO"))))
    expect_identical(getPlatform(v), "visiumhd")
})

test_that("per-FOV summary works from colData, fov_stats or nothing", {
    spe <- .makeTestSPE()
    p <- .generateOutput(FOVStatsPanel(), spe)
    expect_s3_class(p, "ggplot")
    expect_setequal(unique(p$data$name), utils::head(.qcMetrics(spe), 6))
    S4Vectors::metadata(spe)$fov_stats <- data.frame(fov = spe$fov,
        name = "log total", value = log(spe$total))
    p2 <- .generateOutput(FOVStatsPanel(Selection = "1,2"), spe)
    expect_identical(unique(p2$data$name), "log total")
    expect_setequal(levels(p2$data$fov), c("1", "2"))
    p3 <- .generateOutput(FOVStatsPanel(Selection = "99"), spe)
    expect_match(p3$labels$title, "No per-FOV values")
    expect_match(.generateOutput(FOVStatsPanel(),
        .makeTestXeniumSPE())$labels$title, "FOV-based")
})

test_that("panel parameters are checked", {
    expect_error(SpatialStatsPanel(K = 1), "'K'")
    expect_error(SVGPanel(TopN = 0), "'TopN'")
    expect_error(SpatialNichePanel(Lambda = 1.5), "'Lambda'")
    expect_error(SpatialNichePanel(NNiches = 1), "'NNiches'")
    expect_s4_class(SpatialNichePanel(), "SpatialNichePanel")
})

test_that("local clusters and L-functions are drawn", {
    skip_on_cran()
    spe <- .makeTestVisiumHDSPE(n_cells = 90)
    for (pt in c("locclust_plot", "lest_plot", "lcross_plot")) {
        p <- suppressWarnings(suppressMessages(.computeSpatialStat(spe, pt,
            genes = rownames(spe)[1], group_col = "Level1", k = 6L)))
        expect_true(inherits(p, "ggplot"), info = pt)
        title <- c(as.character(p$labels$title), "")[1]
        expect_false(grepl("^(No |Not |Too |Choose)", title), info = title)
    }
})

test_that("local clusters need a gene, not a cell-type column", {
    skip_on_cran()
    spe <- .makeTestVisiumHDSPE(n_cells = 90)
    spe$Level1[1:5] <- NA
    for (grp in list(NA_character_, "Level1")) {
        p <- suppressWarnings(suppressMessages(.computeSpatialStat(spe,
            "locclust_plot", genes = rownames(spe)[1], group_col = grp,
            k = 6L)))
        expect_true(inherits(p, "ggplot"))
        title <- c(as.character(p$labels$title), "")[1]
        expect_false(grepl("^(No |Not |Too |Choose)", title), info = title)
    }
})

test_that("FOV tiles come from the FOV table, else from the cells", {
    spe <- .makeTestSPE()
    r <- .fovRects(spe)
    expect_identical(r$id, .allRegionIDs(spe))
    S4Vectors::metadata(spe)$fov_positions <- data.frame(fov = 4:1,
        x_global_px = c(300, 200, 100, 0), y_global_px = 0)
    S4Vectors::metadata(spe)$fov_dim <- c(xdim = 100, ydim = 50)
    r <- .fovRects(spe)
    expect_identical(r$id, c("1", "2", "3", "4"))
    expect_equal(r$xmin, c(0, 100, 200, 300))
    expect_equal(r$ymax, rep(50, 4))
})
