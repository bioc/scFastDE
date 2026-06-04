#' @title Filter Donors with Too Few Cells per Cell Type
#'
#' @description
#' Removes or flags donors that have fewer than \code{min_cells} cells
#' for a given cell type. When forming pseudo-bulk profiles, donors with
#' very few cells produce highly variable, unreliable aggregated counts
#' that can inflate false-positive DE calls. This function provides a
#' principled pre-filtering step before calling \code{\link{fastPseudobulk}}.
#'
#' @details
#' For each combination of \code{cell_type} level and \code{donor} level,
#' the function counts the number of cells in the \code{sce} and compares
#' it to \code{min_cells}. Donors below the threshold are either removed
#' (when \code{action = "remove"}) or flagged in a new \code{colData}
#' column (when \code{action = "flag"}).
#'
#' As a rule of thumb, \code{min_cells = 10} is a reasonable default for
#' common cell types. For rare cell types (< 1\% frequency), consider
#' lowering to \code{min_cells = 5}.
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'   object.
#' @param donor A \code{character(1)} naming the donor column in
#'   \code{colData(sce)}.
#' @param cell_type A \code{character(1)} naming the cell type column in
#'   \code{colData(sce)}.
#' @param min_cells A \code{integer(1)} minimum number of cells a donor
#'   must have for a given cell type to be retained. Default: \code{10}.
#' @param action A \code{character(1)}, either \code{"remove"} (drop
#'   sparse cells from the SCE) or \code{"flag"} (add a logical column
#'   \code{scFastDE_sparse} to \code{colData}). Default: \code{"remove"}.
#'
#' @return The input \code{sce} with sparse donor-cell type combinations
#'   either removed or flagged depending on \code{action}.
#'   When \code{action = "flag"}, a column \code{scFastDE_sparse} is
#'   added to \code{colData}: \code{TRUE} means the cell belongs to a
#'   sparse donor-cell type combination.
#'
#' @examples
#' library(SingleCellExperiment)
#'
#' set.seed(42)
#' counts <- matrix(rpois(500 * 80, 5), nrow = 500, ncol = 80)
#' rownames(counts) <- paste0("Gene", seq_len(500))
#' colnames(counts) <- paste0("Cell", seq_len(80))
#' sce <- SingleCellExperiment(assays = list(counts = counts))
#' sce$donor     <- rep(paste0("D", 1:8), each = 10)
#' sce$cell_type <- rep(c("Tcell", "Bcell"), times = 40)
#'
#' # Remove donor-cell type combos with fewer than 5 cells
#' sce_filtered <- filterSparseDonors(sce, donor = "donor",
#'                                     cell_type = "cell_type",
#'                                     min_cells = 5)
#' ncol(sce_filtered)
#'
#' @seealso \code{\link{fastPseudobulk}}, \code{\link{fastDE}}
#'
#' @importFrom SummarizedExperiment colData "colData<-"
#' @importFrom methods is
#' @export
filterSparseDonors <- function(sce,
                                donor     = NULL,
                                cell_type = NULL,
                                min_cells = 10L,
                                action    = c("remove", "flag")) {
    action <- match.arg(action)

    .validate_sce_coldata(sce, list(donor = donor, cell_type = cell_type))
    if (!is.numeric(min_cells) || min_cells < 1)
        stop("'min_cells' must be a positive integer.")

    donors     <- as.character(colData(sce)[[donor]])
    cell_types <- as.character(colData(sce)[[cell_type]])
    combo      <- paste(donors, cell_types, sep = "___")

    # Count cells per donor-cell_type combination
    combo_counts <- table(combo)

    # Flag sparse combinations
    is_sparse <- combo_counts[combo] < min_cells
    n_sparse  <- sum(is_sparse)

    if (n_sparse == 0) {
        message("No sparse donor-cell type combinations found at ",
                "min_cells = ", min_cells, ".")
        return(sce)
    }

    message(sprintf(
        "%d cells in %d donor-cell type combinations below min_cells=%d.",
        n_sparse,
        sum(combo_counts < min_cells),
        min_cells
    ))

    if (action == "flag") {
        colData(sce)[["scFastDE_sparse"]] <- as.logical(is_sparse)
        return(sce)
    }

    # action == "remove"
    sce[, !is_sparse]
}
