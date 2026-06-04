#' @title Build Donor-Weighted Pseudo-Bulk Profiles
#'
#' @description
#' Aggregates single-cell counts into pseudo-bulk profiles for a specified
#' cell type, using vectorised sparse matrix operations.
#'
#' When \code{condition} is provided, aggregation is performed per
#' donor-condition pair (one column per sample, e.g. \code{D1_ctrl},
#' \code{D1_stim}). This is the correct approach for paired experimental
#' designs where the same donors contribute cells under multiple conditions.
#'
#' When \code{condition} is \code{NULL} (default), aggregation is performed
#' per donor only (backward-compatible behaviour for unpaired designs).
#'
#' Each sample's weight equals \code{sqrt(n_cells)}, giving more influence
#' to well-represented samples while not discarding those with fewer cells.
#'
#' @details
#' Pseudo-bulk aggregation sums the raw counts for all cells belonging to
#' each donor (or donor-condition pair) within a given cell type:
#'
#' \deqn{PB_{g,s} = \sum_{c \in \text{sample}_s, \text{type}_t} X_{g,c}}
#'
#' Sample weights are then computed as:
#'
#' \deqn{w_s = \sqrt{n_{s,t}}}
#'
#' where \eqn{n_{s,t}} is the number of cells for sample \eqn{s} in cell
#' type \eqn{t}.
#'
#' The aggregation uses sparse matrix column-sums grouped by sample,
#' avoiding a slow R-level loop.
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'   object with raw counts in \code{assay(sce, "counts")}.
#' @param donor A \code{character(1)} naming the donor column in
#'   \code{colData(sce)}.
#' @param cell_type A \code{character(1)} naming the cell type column in
#'   \code{colData(sce)}.
#' @param target_type A \code{character(1)} specifying which cell type to
#'   aggregate. Must be a value present in \code{colData(sce)[[cell_type]]}.
#' @param condition A \code{character(1)} naming the condition column in
#'   \code{colData(sce)}. When provided, aggregation is done per
#'   donor-condition pair (required for paired designs). Default:
#'   \code{NULL} (aggregate per donor only).
#' @param assay_name A \code{character(1)} name of the assay to aggregate.
#'   Default: \code{"counts"}.
#'
#' @return A \code{list} with elements:
#'   \describe{
#'     \item{\code{pseudobulk}}{A \code{matrix} of aggregated counts,
#'       rows = genes, columns = samples.}
#'     \item{\code{sample_weights}}{A named \code{numeric} vector of
#'       per-sample weights (\code{sqrt(n_cells)}).}
#'     \item{\code{sample_ncells}}{A named \code{integer} vector giving
#'       the raw cell count per sample.}
#'     \item{\code{sample_info}}{A \code{data.frame} with columns
#'       \code{sample_id}, \code{donor}, and (if applicable)
#'       \code{condition} — one row per pseudo-bulk sample.}
#'   }
#'   For backward compatibility, \code{donor_weights} and
#'   \code{donor_ncells} are also provided as aliases when
#'   \code{condition} is \code{NULL}.
#'
#' @examples
#' library(SingleCellExperiment)
#'
#' set.seed(42)
#' counts <- matrix(rpois(500 * 60, 5), nrow = 500, ncol = 60)
#' rownames(counts) <- paste0("Gene", seq_len(500))
#' colnames(counts) <- paste0("Cell", seq_len(60))
#' sce <- SingleCellExperiment(assays = list(counts = counts))
#' sce$donor     <- rep(paste0("D", 1:6), each = 10)
#' sce$cell_type <- "Tcell"
#' sce$condition <- rep(c("ctrl", "treat"), each = 30)
#'
#' # Unpaired (aggregate by donor only)
#' pb <- fastPseudobulk(sce, donor = "donor",
#'                       cell_type = "cell_type",
#'                       target_type = "Tcell")
#' dim(pb$pseudobulk)
#'
#' # Paired (aggregate by donor x condition)
#' pb2 <- fastPseudobulk(sce, donor = "donor",
#'                        cell_type = "cell_type",
#'                        target_type = "Tcell",
#'                        condition = "condition")
#' dim(pb2$pseudobulk)
#' pb2$sample_info
#'
#' @seealso \code{\link{filterSparseDonors}}, \code{\link{fastDE}}
#'
#' @importFrom SummarizedExperiment assay colData
#' @importFrom methods is
#' @importFrom stats setNames
#' @export
fastPseudobulk <- function(sce,
                            donor       = NULL,
                            cell_type   = NULL,
                            target_type = NULL,
                            condition   = NULL,
                            assay_name  = "counts") {
    # ── Validation ────────────────────────────────────────────────────────────
    .validate_sce_coldata(sce, list(donor = donor, cell_type = cell_type))
    if (is.null(target_type))
        stop("'target_type' must be specified.")
    if (!assay_name %in% names(assays(sce)))
        stop("Assay '", assay_name, "' not found in sce.")
    if (!is.null(condition))
        .validate_sce_coldata(sce, list(condition = condition))

    ct_vals <- as.character(colData(sce)[[cell_type]])
    if (!target_type %in% ct_vals)
        stop("'target_type' = '", target_type, "' not found in ",
             "'", cell_type, "' column.")

    # ── Subset to target cell type ────────────────────────────────────────────
    sce_sub    <- sce[, ct_vals == target_type]
    donors_sub <- as.character(colData(sce_sub)[[donor]])

    # ── Build sample labels ───────────────────────────────────────────────────
    if (!is.null(condition)) {
        # Paired: aggregate by donor × condition
        cond_sub    <- as.character(colData(sce_sub)[[condition]])
        sample_ids  <- paste(donors_sub, cond_sub, sep = "___")
        sample_lvls <- unique(sample_ids)
        n_samples   <- length(sample_lvls)

        # Build sample_info
        parts <- strsplit(sample_lvls, "___", fixed = TRUE)
        sample_info <- data.frame(
            sample_id = sample_lvls,
            donor     = vapply(parts, `[`, character(1), 1L),
            condition = vapply(parts, `[`, character(1), 2L),
            stringsAsFactors = FALSE
        )
    } else {
        # Unpaired: aggregate by donor only
        sample_ids  <- donors_sub
        sample_lvls <- unique(donors_sub)
        n_samples   <- length(sample_lvls)

        sample_info <- data.frame(
            sample_id = sample_lvls,
            donor     = sample_lvls,
            stringsAsFactors = FALSE
        )
    }

    if (n_samples < 2)
        stop("At least 2 samples required for DE analysis. Found: ",
             n_samples, " sample(s) for cell type '", target_type, "'.")

    # ── Vectorised pseudo-bulk aggregation ────────────────────────────────────
    counts_mat <- assay(sce_sub, assay_name)
    sample_idx <- match(sample_ids, sample_lvls)

    # Sparse indicator matrix (cells x samples)
    indicator <- Matrix::sparseMatrix(
        i = seq_along(sample_idx),
        j = sample_idx,
        x = 1,
        dims = c(ncol(sce_sub), n_samples),
        dimnames = list(colnames(sce_sub), sample_lvls)
    )

    # Vectorised aggregation: genes x samples — one matrix multiply
    pb_matrix <- as.matrix(counts_mat %*% indicator)

    # ── Sample weights: sqrt(n_cells) ─────────────────────────────────────────
    n_cells_per_sample <- tabulate(sample_idx, nbins = n_samples)
    names(n_cells_per_sample) <- sample_lvls
    sample_weights <- sqrt(n_cells_per_sample)

    result <- list(
        pseudobulk     = pb_matrix,
        sample_weights = sample_weights,
        sample_ncells  = n_cells_per_sample,
        sample_info    = sample_info
    )

    # Backward-compatible aliases when condition is NULL
    if (is.null(condition)) {
        result$donor_weights <- sample_weights
        result$donor_ncells  <- n_cells_per_sample
    }

    result
}


# ── Internal helper ───────────────────────────────────────────────────────────

#' Access assays from SCE (internal helper)
#' @keywords internal
#' @noRd
assays <- function(x) SummarizedExperiment::assays(x)
