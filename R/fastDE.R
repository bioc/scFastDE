#' @title Fast Vectorised Differential Expression
#'
#' @description
#' Runs differential expression analysis across all genes simultaneously
#' using a weighted limma-voom linear model. Unlike tools that loop over
#' genes serially, \code{fastDE} passes the entire pseudo-bulk matrix to
#' limma, which uses LAPACK routines to fit all gene models in one
#' vectorised call.
#'
#' The function automatically detects whether the experimental design is
#' \strong{paired} (same donors in multiple conditions) or
#' \strong{unpaired} (each donor in one condition only) and builds the
#' appropriate linear model:
#' \itemize{
#'   \item \strong{Paired:} pseudo-bulk is aggregated per donor-condition
#'     pair, and a \code{~ 0 + condition + donor} model accounts for
#'     donor-level variation.
#'   \item \strong{Unpaired:} pseudo-bulk is aggregated per donor, and a
#'     \code{~ 0 + condition} model is used.
#' }
#'
#' @details
#' The analysis pipeline is:
#' \enumerate{
#'   \item Filter lowly-expressed genes (CPM > \code{min_cpm} in at least
#'     \code{min_donors} samples).
#'   \item Apply \code{limma::voom} with sample cell-count weights.
#'   \item Fit a weighted linear model with the appropriate design.
#'   \item Extract the contrast of interest with \code{limma::makeContrasts}.
#'   \item Apply empirical Bayes moderation with \code{limma::eBayes}.
#'   \item Return all results as a \code{\link{FDEResult}} object.
#' }
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'   object with raw counts in \code{assay(sce, "counts")}.
#' @param donor A \code{character(1)} naming the donor column in
#'   \code{colData(sce)}.
#' @param cell_type A \code{character(1)} naming the cell type column in
#'   \code{colData(sce)}.
#' @param condition A \code{character(1)} naming the condition column in
#'   \code{colData(sce)} (e.g. \code{"treatment"}, \code{"disease"}).
#'   Must have exactly two unique values.
#' @param target_type A \code{character(1)} specifying which cell type
#'   to test. Must be present in \code{colData(sce)[[cell_type]]}.
#' @param contrast A \code{character(1)} specifying the contrast in
#'   limma syntax, e.g. \code{"treat - ctrl"}. If \code{NULL}, the
#'   second level minus the first level is used. Default: \code{NULL}.
#' @param min_cpm A \code{numeric(1)} minimum CPM for a gene to be
#'   considered expressed. Default: \code{1}.
#' @param min_donors An \code{integer(1)} minimum number of samples in
#'   which a gene must exceed \code{min_cpm}. Default: \code{2}.
#' @param min_cells An \code{integer(1)} minimum cells per donor before
#'   pseudo-bulk aggregation. Passed to \code{\link{filterSparseDonors}}.
#'   Default: \code{10}.
#' @param BPPARAM A \code{\link[BiocParallel]{BiocParallelParam}} object.
#'   Default: \code{SerialParam()}.
#'
#' @return A \code{\link{FDEResult}} object containing:
#'   \itemize{
#'     \item \code{deTable}: per-gene statistics (\code{logFC},
#'       \code{AveExpr}, \code{t}, \code{P.Value}, \code{adj.P.Val},
#'       \code{B}).
#'     \item \code{pseudobulk}: the aggregated count matrix.
#'     \item \code{donorWeights}: the per-sample \code{sqrt(n_cells)}
#'       weights.
#'     \item \code{params}: the analysis parameters.
#'   }
#'
#' @examples
#' library(SingleCellExperiment)
#'
#' set.seed(42)
#' n_genes <- 200
#' counts <- matrix(rpois(n_genes * 60, 8), nrow = n_genes, ncol = 60)
#' rownames(counts) <- paste0("Gene", seq_len(n_genes))
#' colnames(counts) <- paste0("Cell", seq_len(60))
#'
#' # Inject DE signal into first 10 genes for treatment group
#' counts[1:10, 31:60] <- counts[1:10, 31:60] * 3L
#'
#' sce <- SingleCellExperiment(assays = list(counts = counts))
#' sce$donor     <- rep(paste0("D", 1:6), each = 10)
#' sce$cell_type <- "Tcell"
#' sce$condition <- rep(c("ctrl", "treat"), each = 30)
#'
#' result <- fastDE(sce,
#'                  donor       = "donor",
#'                  cell_type   = "cell_type",
#'                  condition   = "condition",
#'                  target_type = "Tcell",
#'                  min_cells   = 5)
#' result
#' head(deTable(result))
#'
#' @seealso \code{\link{fastPseudobulk}}, \code{\link{filterSparseDonors}},
#'   \code{\link{plotDEResults}}
#'
#' @importFrom limma voom lmFit makeContrasts contrasts.fit eBayes topTable
#' @importFrom SummarizedExperiment colData
#' @importFrom S4Vectors DataFrame
#' @importFrom BiocParallel SerialParam
#' @importFrom methods is
#' @importFrom stats model.matrix
#' @export
fastDE <- function(sce,
                    donor       = NULL,
                    cell_type   = NULL,
                    condition   = NULL,
                    target_type = NULL,
                    contrast    = NULL,
                    min_cpm     = 1,
                    min_donors  = 2L,
                    min_cells   = 10L,
                    BPPARAM     = SerialParam()) {
    # ── Validation ────────────────────────────────────────────────────────────
    if (!is(sce, "SingleCellExperiment"))
        stop("'sce' must be a SingleCellExperiment object.")
    for (col in c(donor, cell_type, condition)) {
        if (is.null(col) || !col %in% names(colData(sce)))
            stop("Column '", col, "' not found in colData(sce).")
    }
    if (is.null(target_type))
        stop("'target_type' must be specified.")

    # ── Step 1: filter sparse donors ─────────────────────────────────────────
    sce <- filterSparseDonors(sce,
                               donor     = donor,
                               cell_type = cell_type,
                               min_cells = min_cells,
                               action    = "remove")

    # ── Step 2: detect paired vs unpaired design ─────────────────────────────
    sce_sub   <- sce[, as.character(colData(sce)[[cell_type]]) ==
                     target_type]
    donor_sub <- as.character(colData(sce_sub)[[donor]])
    cond_sub  <- as.character(colData(sce_sub)[[condition]])

    # Validate exactly 2 condition levels
    all_cond_levels <- unique(cond_sub)
    if (length(all_cond_levels) != 2)
        stop("'condition' must have exactly 2 levels. Found: ",
             paste(sort(all_cond_levels), collapse = ", "))

    # For each donor, find which conditions they appear in
    donor_conds <- tapply(cond_sub, donor_sub, function(x) {
        unique(x)
    })
    is_paired <- any(vapply(donor_conds, length, integer(1)) > 1L)

    if (is_paired) {
        message("Detected paired design (same donors in multiple ",
                "conditions). Aggregating by donor x condition.")
    } else {
        message("Detected unpaired design (each donor in one ",
                "condition).")
    }

    # ── Step 3: build pseudo-bulk ─────────────────────────────────────────────
    if (is_paired) {
        pb <- fastPseudobulk(sce,
                              donor       = donor,
                              cell_type   = cell_type,
                              target_type = target_type,
                              condition   = condition)
    } else {
        pb <- fastPseudobulk(sce,
                              donor       = donor,
                              cell_type   = cell_type,
                              target_type = target_type)
    }

    pb_mat      <- pb$pseudobulk
    weights     <- pb$sample_weights
    sample_info <- pb$sample_info
    n_samples   <- ncol(pb_mat)

    # ── Step 4: build design matrix ───────────────────────────────────────────
    if (is_paired) {
        # Paired: use ~ 0 + condition + donor to block on donor
        cond_f <- factor(sample_info$condition)
        cond_levels <- sort(levels(cond_f))
        cond_f <- factor(sample_info$condition, levels = cond_levels)

        donor_f <- factor(sample_info$donor)

        # Check we have enough residual df
        n_params <- length(cond_levels) + length(levels(donor_f)) - 1L
        if (n_samples <= n_params)
            stop("Not enough samples (", n_samples,
                 ") for the paired model with ",
                 length(levels(donor_f)), " donors and ",
                 length(cond_levels), " conditions. ",
                 "Need at least ", n_params + 1L, " samples.")

        design <- model.matrix(~ 0 + cond_f + donor_f)
        cond_levels <- make.names(cond_levels, unique = TRUE)
        colnames(design)[seq_along(cond_levels)] <- cond_levels
        # Clean up donor column names — must be valid R names
        donor_cols <- seq(length(cond_levels) + 1L, ncol(design))
        colnames(design)[donor_cols] <- make.names(
            sub("^donor_f", "", colnames(design)[donor_cols]),
            unique = TRUE
        )
    } else {
        # Unpaired: ~ 0 + condition (original approach)
        # Get condition per donor from original SCE data
        donor_condition <- tapply(cond_sub, donor_sub, function(x) {
            names(sort(table(x), decreasing = TRUE))[1]
        })
        donor_condition <- donor_condition[colnames(pb_mat)]

        cond_levels <- sort(unique(donor_condition))
        if (length(cond_levels) != 2)
            stop("'condition' must have exactly 2 levels. Found: ",
                 paste(cond_levels, collapse = ", "))

        cond_f <- factor(donor_condition, levels = cond_levels)
        design <- model.matrix(~ 0 + cond_f)
        colnames(design) <- cond_levels
    }

    # ── Step 5: filter lowly expressed genes ─────────────────────────────────
    lib_sizes <- colSums(pb_mat)
    cpm_mat   <- sweep(pb_mat, 2, lib_sizes / 1e6, FUN = "/")
    keep      <- rowSums(cpm_mat >= min_cpm) >= min_donors

    if (sum(keep) == 0)
        stop("No genes passed the expression filter. ",
             "Try lowering 'min_cpm' or 'min_donors'.")

    pb_filt <- pb_mat[keep, , drop = FALSE]
    message(sprintf("Testing %d / %d genes after expression filter.",
                    sum(keep), nrow(pb_mat)))

    # ── Step 6: voom + weighted linear model ─────────────────────────────────
    v <- voom(pb_filt, design,
              weights  = matrix(rep(weights, each = nrow(pb_filt)),
                                nrow = nrow(pb_filt)),
              plot     = FALSE)

    fit <- lmFit(v, design, weights = v$weights)

    # ── Step 7: contrast ──────────────────────────────────────────────────────
    if (is.null(contrast))
        contrast <- paste(cond_levels[2], "-", cond_levels[1])

    contrast_mat <- makeContrasts(contrasts = contrast, levels = design)
    fit2         <- contrasts.fit(fit, contrast_mat)
    fit2         <- eBayes(fit2)

    # ── Step 8: collect results ───────────────────────────────────────────────
    tt <- topTable(fit2, number = Inf, sort.by = "none")
    de_df <- DataFrame(
        gene      = rownames(tt),
        logFC     = tt[["logFC"]],
        AveExpr   = tt[["AveExpr"]],
        t         = tt[["t"]],
        P.Value   = tt[["P.Value"]],
        adj.P.Val = tt[["adj.P.Val"]],
        B         = tt[["B"]],
        row.names = rownames(tt)
    )

    FDEResult(
        deTable      = de_df,
        pseudobulk   = pb_filt,
        donorWeights = weights,
        params       = list(
            cell_type   = target_type,
            condition   = condition,
            donor       = donor,
            contrast    = contrast,
            min_cpm     = min_cpm,
            min_donors  = min_donors,
            min_cells   = min_cells,
            is_paired   = is_paired
        )
    )
}
