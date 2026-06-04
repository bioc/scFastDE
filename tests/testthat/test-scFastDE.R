library(SingleCellExperiment)
library(scFastDE)

# ── Helper: minimal 2-condition 6-donor SCE ───────────────────────────────────
make_sce <- function(seed = 42, n_genes = 100, n_cells = 60) {
    set.seed(seed)
    counts <- matrix(rpois(n_genes * n_cells, 8L),
                     nrow = n_genes, ncol = n_cells)
    rownames(counts) <- paste0("Gene", seq_len(n_genes))
    colnames(counts) <- paste0("Cell", seq_len(n_cells))
    # Inject DE signal into first 10 genes for treatment
    counts[1:10, 31:n_cells] <- counts[1:10, 31:n_cells] * 3L
    sce <- SingleCellExperiment(assays = list(counts = counts))
    sce$donor     <- rep(paste0("D", 1:6), each = 10)
    sce$cell_type <- "Tcell"
    sce$condition <- rep(c("ctrl", "treat"), each = 30)
    sce
}

# ══════════════════════════════════════════════════════════════════════════════
# filterSparseDonors
# ══════════════════════════════════════════════════════════════════════════════

test_that("filterSparseDonors returns a SingleCellExperiment", {
    sce <- make_sce()
    res <- filterSparseDonors(sce, donor = "donor",
                               cell_type = "cell_type", min_cells = 5)
    testthat::expect_s4_class(res, "SingleCellExperiment")
})

test_that("filterSparseDonors removes donors below threshold", {
    sce <- make_sce()
    # All donors have 10 cells; require 15 to force removal
    res <- filterSparseDonors(sce, donor = "donor",
                               cell_type = "cell_type", min_cells = 15)
    testthat::expect_lt(ncol(res), ncol(sce))
})

test_that("filterSparseDonors with action='flag' adds colData column", {
    sce <- make_sce()
    res <- filterSparseDonors(sce, donor = "donor",
                               cell_type = "cell_type",
                               min_cells = 15, action = "flag")
    testthat::expect_true("scFastDE_sparse" %in% names(colData(res)))
    testthat::expect_type(res$scFastDE_sparse, "logical")
})

test_that("filterSparseDonors errors on bad column name", {
    sce <- make_sce()
    testthat::expect_error(
        filterSparseDonors(sce, donor = "nonexistent",
                            cell_type = "cell_type"),
        regexp = "not found in colData"
    )
})

test_that("filterSparseDonors returns unchanged SCE when no sparse donors", {
    sce <- make_sce()
    res <- filterSparseDonors(sce, donor = "donor",
                               cell_type = "cell_type", min_cells = 1)
    testthat::expect_equal(ncol(res), ncol(sce))
})

# ══════════════════════════════════════════════════════════════════════════════
# fastPseudobulk
# ══════════════════════════════════════════════════════════════════════════════

test_that("fastPseudobulk returns a list with correct elements", {
    sce <- make_sce()
    pb <- fastPseudobulk(sce, donor = "donor",
                          cell_type = "cell_type", target_type = "Tcell")
    testthat::expect_type(pb, "list")
    testthat::expect_true(all(c("pseudobulk", "donor_weights",
                       "donor_ncells") %in% names(pb)))
})

test_that("fastPseudobulk pseudobulk has correct dimensions", {
    sce <- make_sce()
    pb <- fastPseudobulk(sce, donor = "donor",
                          cell_type = "cell_type", target_type = "Tcell")
    testthat::expect_equal(nrow(pb$pseudobulk), nrow(sce))
    testthat::expect_equal(ncol(pb$pseudobulk), length(unique(sce$donor)))
})

test_that("fastPseudobulk donor_weights are sqrt of cell counts", {
    sce <- make_sce()
    pb <- fastPseudobulk(sce, donor = "donor",
                          cell_type = "cell_type", target_type = "Tcell")
    expected <- sqrt(pb$donor_ncells)
    testthat::expect_equal(pb$donor_weights, expected)
})

test_that("fastPseudobulk errors on invalid target_type", {
    sce <- make_sce()
    testthat::expect_error(
        fastPseudobulk(sce, donor = "donor",
                        cell_type = "cell_type", target_type = "NK"),
        regexp = "not found"
    )
})

test_that("fastPseudobulk errors when fewer than 2 donors", {
    sce <- make_sce()
    sce_1d <- sce[, sce$donor == "D1"]
    testthat::expect_error(
        fastPseudobulk(sce_1d, donor = "donor",
                        cell_type = "cell_type", target_type = "Tcell"),
        regexp = "At least 2"
    )
})

# ══════════════════════════════════════════════════════════════════════════════
# fastDE
# ══════════════════════════════════════════════════════════════════════════════

test_that("fastDE returns a FDEResult", {
    sce <- make_sce()
    res <- fastDE(sce, donor = "donor", cell_type = "cell_type",
                  condition = "condition", target_type = "Tcell",
                  min_cells = 5)
    testthat::expect_s4_class(res, "FDEResult")
})

test_that("fastDE deTable has expected columns", {
    sce <- make_sce()
    res <- fastDE(sce, donor = "donor", cell_type = "cell_type",
                  condition = "condition", target_type = "Tcell",
                  min_cells = 5)
    testthat::expect_true(all(c("logFC", "P.Value", "adj.P.Val") %in%
                    names(deTable(res))))
})

test_that("fastDE detects injected DE signal", {
    sce <- make_sce(seed = 42)
    res <- fastDE(sce, donor = "donor", cell_type = "cell_type",
                  condition = "condition", target_type = "Tcell",
                  min_cells = 5)
    dt <- as.data.frame(deTable(res))
    sig_genes <- rownames(dt)[dt$adj.P.Val < 0.05]
    # At least some of the first 10 injected genes should be detected
    detected <- intersect(sig_genes, paste0("Gene", 1:10))
    testthat::expect_gt(length(detected), 0)
})

test_that("fastDE errors on non-SCE input", {
    testthat::expect_error(
        fastDE(list(), donor = "donor", cell_type = "cell_type",
               condition = "condition", target_type = "Tcell"),
        regexp = "SingleCellExperiment"
    )
})

test_that("fastDE errors when condition has more than 2 levels", {
    sce <- make_sce()
    sce$condition <- rep(c("a", "b", "c"), length.out = ncol(sce))
    testthat::expect_error(
        fastDE(sce, donor = "donor", cell_type = "cell_type",
               condition = "condition", target_type = "Tcell",
               min_cells = 1),
        regexp = "exactly 2 levels"
    )
})

# ══════════════════════════════════════════════════════════════════════════════
# plotDEResults
# ══════════════════════════════════════════════════════════════════════════════

test_that("plotDEResults returns a ggplot object", {
    sce <- make_sce()
    res <- fastDE(sce, donor = "donor", cell_type = "cell_type",
                  condition = "condition", target_type = "Tcell",
                  min_cells = 5)
    p <- plotDEResults(res)
    testthat::expect_s3_class(p, "ggplot")
})

test_that("plotDEResults errors on non-FDEResult input", {
    testthat::expect_error(plotDEResults(list()), regexp = "FDEResult")
})

# ══════════════════════════════════════════════════════════════════════════════
# FDEResult S4 class and accessors
# ══════════════════════════════════════════════════════════════════════════════

test_that("FDEResult constructor and accessors work", {
    library(S4Vectors)
    de <- DataFrame(
        gene      = c("G1", "G2"),
        logFC     = c(1.2, -0.5),
        P.Value   = c(0.001, 0.5),
        adj.P.Val = c(0.01, 0.8),
        AveExpr   = c(3.1, 2.2),
        t         = c(4.1, 0.8),
        B         = c(2.1, -1.2)
    )
    pb <- matrix(rpois(20, 10), nrow = 2, ncol = 10)
    rownames(pb) <- c("G1", "G2")
    colnames(pb) <- paste0("D", 1:10)
    wts <- setNames(sqrt(1:10), paste0("D", 1:10))

    obj <- FDEResult(deTable = de, pseudobulk = pb,
                     donorWeights = wts,
                     params = list(cell_type = "Tcell"))

    testthat::expect_s4_class(obj, "FDEResult")
    testthat::expect_equal(nrow(deTable(obj)), 2)
    testthat::expect_equal(dim(pseudobulk(obj)), c(2, 10))
    testthat::expect_length(donorWeights(obj), 10)
})

test_that("show method for FDEResult prints without error", {
    library(S4Vectors)
    de <- DataFrame(
        gene = "G1", logFC = 1.2, P.Value = 0.001,
        adj.P.Val = 0.01, AveExpr = 3.1, t = 4.1, B = 2.1
    )
    pb <- matrix(rpois(10, 10), nrow = 1, ncol = 10)
    rownames(pb) <- "G1"
    colnames(pb) <- paste0("D", 1:10)
    obj <- FDEResult(de, pb, setNames(sqrt(1:10), paste0("D", 1:10)))
    testthat::expect_output(show(obj), "FDEResult")
})

# ══════════════════════════════════════════════════════════════════════════════
# Paired design support
# ══════════════════════════════════════════════════════════════════════════════

# Helper: paired SCE where SAME donors appear in BOTH conditions
make_paired_sce <- function(seed = 99, n_genes = 100, n_donors = 6,
                            cells_per_donor_cond = 10) {
    set.seed(seed)
    n_cells <- n_donors * 2 * cells_per_donor_cond
    counts <- matrix(rpois(n_genes * n_cells, 8L),
                     nrow = n_genes, ncol = n_cells)
    rownames(counts) <- paste0("Gene", seq_len(n_genes))
    colnames(counts) <- paste0("Cell", seq_len(n_cells))

    # Inject DE signal: first 15 genes upregulated in "treat"
    treat_cols <- seq(n_donors * cells_per_donor_cond + 1, n_cells)
    counts[1:15, treat_cols] <- counts[1:15, treat_cols] * 4L

    sce <- SingleCellExperiment(assays = list(counts = counts))
    # Same donors in both conditions
    sce$donor     <- rep(rep(paste0("D", seq_len(n_donors)),
                            each = cells_per_donor_cond), 2)
    sce$cell_type <- "Tcell"
    sce$condition <- rep(c("ctrl", "treat"),
                         each = n_donors * cells_per_donor_cond)
    sce
}

test_that("fastPseudobulk with condition creates donor x cond samples", {
    sce <- make_paired_sce()
    pb <- fastPseudobulk(sce, donor = "donor",
                          cell_type = "cell_type",
                          target_type = "Tcell",
                          condition = "condition")
    # 6 donors x 2 conditions = 12 samples
    testthat::expect_equal(ncol(pb$pseudobulk), 12)
    testthat::expect_true("sample_info" %in% names(pb))
    testthat::expect_equal(nrow(pb$sample_info), 12)
    testthat::expect_true(all(c("donor", "condition") %in% names(pb$sample_info)))
})

test_that("fastDE detects paired design automatically", {
    sce <- make_paired_sce()
    res <- fastDE(sce, donor = "donor", cell_type = "cell_type",
                  condition = "condition", target_type = "Tcell",
                  min_cells = 5)
    testthat::expect_s4_class(res, "FDEResult")
    testthat::expect_true(params(res)$is_paired)
    # Should produce 12 samples (6 donors x 2 conditions)
    testthat::expect_equal(ncol(pseudobulk(res)), 12)
})

test_that("fastDE detects injected DE signal in paired design", {
    sce <- make_paired_sce(seed = 42)
    res <- fastDE(sce, donor = "donor", cell_type = "cell_type",
                  condition = "condition", target_type = "Tcell",
                  min_cells = 5)
    dt <- as.data.frame(deTable(res))
    sig_genes <- rownames(dt)[dt$adj.P.Val < 0.05]
    # At least some of the 15 injected DE genes should be detected
    detected <- intersect(sig_genes, paste0("Gene", 1:15))
    testthat::expect_gt(length(detected), 0)
})

test_that("unpaired design still works after paired fix", {
    sce <- make_sce()
    res <- fastDE(sce, donor = "donor", cell_type = "cell_type",
                  condition = "condition", target_type = "Tcell",
                  min_cells = 5)
    testthat::expect_s4_class(res, "FDEResult")
    testthat::expect_false(params(res)$is_paired)
    dt <- as.data.frame(deTable(res))
    sig_genes <- rownames(dt)[dt$adj.P.Val < 0.05]
    detected <- intersect(sig_genes, paste0("Gene", 1:10))
    testthat::expect_gt(length(detected), 0)
})
