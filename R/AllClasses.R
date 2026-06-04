#' @title FDEResult: Fast DE Result Container
#'
#' @description
#' An S4 class to store the output of \code{\link{fastDE}}. Slots hold
#' per-gene DE statistics, the pseudo-bulk matrix used, donor weights,
#' and the analysis parameters.
#'
#' @slot deTable A \code{DataFrame} with one row per gene containing
#'   columns \code{logFC}, \code{AveExpr}, \code{t}, \code{P.Value},
#'   \code{adj.P.Val}, and \code{B}.
#' @slot pseudobulk A \code{matrix} of pseudo-bulk counts (genes x donors)
#'   used as input to the DE model.
#' @slot donorWeights A named \code{numeric} vector of per-donor weights
#'   (sqrt of cell count) used in the linear model.
#' @slot params A \code{list} of analysis parameters including
#'   \code{cell_type}, \code{condition}, \code{donor}, \code{min_cells}.
#'
#' @exportClass FDEResult
#' @importFrom methods new is
#' @importFrom S4Vectors DataFrame
setClass(
    "FDEResult",
    representation(
        deTable      = "DataFrame",
        pseudobulk   = "matrix",
        donorWeights = "numeric",
        params       = "list"
    )
)

#' @title Constructor for FDEResult
#'
#' @description Create a new \code{FDEResult} object.
#'
#' @param deTable A \code{DataFrame} of per-gene DE statistics.
#' @param pseudobulk A \code{matrix} of pseudo-bulk counts.
#' @param donorWeights A named \code{numeric} vector of donor weights.
#' @param params A \code{list} of analysis parameters.
#'
#' @return A \code{FDEResult} object.
#'
#' @examples
#' library(S4Vectors)
#' de <- DataFrame(
#'     logFC     = c(1.2, -0.5, 0.8),
#'     P.Value   = c(0.001, 0.5, 0.02),
#'     adj.P.Val = c(0.01, 0.8, 0.05)
#' )
#' pb <- matrix(rpois(30, 10), nrow = 3, ncol = 10)
#' rownames(pb) <- paste0("Gene", 1:3)
#' colnames(pb) <- paste0("Donor", 1:10)
#' obj <- FDEResult(
#'     deTable      = de,
#'     pseudobulk   = pb,
#'     donorWeights = setNames(sqrt(1:10), paste0("Donor", 1:10)),
#'     params       = list(cell_type = "T_cells", condition = "group")
#' )
#' obj
#'
#' @export
FDEResult <- function(deTable, pseudobulk, donorWeights, params = list()) {
    new("FDEResult",
        deTable      = deTable,
        pseudobulk   = pseudobulk,
        donorWeights = donorWeights,
        params       = params)
}

# ── Generics ──────────────────────────────────────────────────────────────────

#' @title Accessor for DE table in a FDEResult
#'
#' @description Returns the per-gene differential expression statistics
#'   \code{DataFrame} from a \code{FDEResult} object.
#'
#' @param x A \code{FDEResult} object.
#' @param ... Additional arguments (not used).
#'
#' @return A \code{DataFrame} with columns \code{logFC}, \code{AveExpr},
#'   \code{t}, \code{P.Value}, \code{adj.P.Val}, and \code{B}.
#'
#' @examples
#' library(S4Vectors)
#' de <- DataFrame(logFC = c(1.2, -0.5), P.Value = c(0.001, 0.5),
#'                 adj.P.Val = c(0.01, 0.8))
#' pb <- matrix(rpois(20, 10), nrow = 2, ncol = 10)
#' rownames(pb) <- paste0("Gene", 1:2)
#' colnames(pb) <- paste0("Donor", 1:10)
#' obj <- FDEResult(de, pb, setNames(sqrt(1:10), paste0("Donor", 1:10)))
#' deTable(obj)
#'
#' @export
setGeneric("deTable", function(x, ...) standardGeneric("deTable"))

#' @title Accessor for pseudo-bulk matrix in a FDEResult
#'
#' @description Returns the pseudo-bulk count matrix from a
#'   \code{FDEResult} object.
#'
#' @param x A \code{FDEResult} object.
#' @param ... Additional arguments (not used).
#'
#' @return A \code{matrix} of pseudo-bulk counts (genes x donors).
#'
#' @examples
#' library(S4Vectors)
#' de <- DataFrame(logFC = c(1.2, -0.5), P.Value = c(0.001, 0.5),
#'                 adj.P.Val = c(0.01, 0.8))
#' pb <- matrix(rpois(20, 10), nrow = 2, ncol = 10)
#' rownames(pb) <- paste0("Gene", 1:2)
#' colnames(pb) <- paste0("Donor", 1:10)
#' obj <- FDEResult(de, pb, setNames(sqrt(1:10), paste0("Donor", 1:10)))
#' pseudobulk(obj)
#'
#' @export
setGeneric("pseudobulk", function(x, ...) standardGeneric("pseudobulk"))

#' @title Accessor for donor weights in a FDEResult
#'
#' @description Returns the per-donor weight vector from a
#'   \code{FDEResult} object.
#'
#' @param x A \code{FDEResult} object.
#' @param ... Additional arguments (not used).
#'
#' @return A named \code{numeric} vector of donor weights.
#'
#' @examples
#' library(S4Vectors)
#' de <- DataFrame(logFC = c(1.2, -0.5), P.Value = c(0.001, 0.5),
#'                 adj.P.Val = c(0.01, 0.8))
#' pb <- matrix(rpois(20, 10), nrow = 2, ncol = 10)
#' rownames(pb) <- paste0("Gene", 1:2)
#' colnames(pb) <- paste0("Donor", 1:10)
#' obj <- FDEResult(de, pb, setNames(sqrt(1:10), paste0("Donor", 1:10)))
#' donorWeights(obj)
#'
#' @export
setGeneric("donorWeights", function(x, ...) standardGeneric("donorWeights"))



# ── Methods ───────────────────────────────────────────────────────────────────

#' @title Show method for FDEResult
#'
#' @description Prints a compact summary of a \code{FDEResult} object.
#'
#' @param object A \code{FDEResult} object.
#'
#' @return Invisibly returns \code{object}.
#'
#' @importFrom methods show
#' @export
setMethod("show", "FDEResult", function(object) {
    dt <- deTable(object)
    pb <- pseudobulk(object)
    object_params <- params(object)

    sig <- sum(dt[["adj.P.Val"]] < 0.05, na.rm = TRUE)
    is_paired <- isTRUE(object_params$is_paired)
    cat("FDEResult\n")
    cat("  Genes tested   :", nrow(dt), "\n")
    cat("  Samples        :", ncol(pb), "\n")
    cat("  Significant    :", sig, "(adj.P.Val < 0.05)\n")
    if (!is.null(object_params$cell_type))
        cat("  Cell type      :", object_params$cell_type, "\n")
    if (!is.null(object_params$condition))
        cat("  Condition      :", object_params$condition, "\n")
    cat("  Design         :", if (is_paired) "paired" else "unpaired",
        "\n")
    invisible(object)
})

#' @rdname deTable
#' @export
setMethod("deTable", "FDEResult", function(x, ...) x@deTable)

#' @rdname pseudobulk
#' @export
setMethod("pseudobulk", "FDEResult", function(x, ...) x@pseudobulk)

#' @rdname donorWeights
#' @export
setMethod("donorWeights", "FDEResult", function(x, ...) x@donorWeights)

#' @title Accessor for analysis parameters in a FDEResult
#'
#' @description Returns the analysis parameter list from a
#'   \code{FDEResult} object.
#'
#' @param x A \code{FDEResult} object.
#' @param ... Additional arguments (not used).
#'
#' @return A \code{list} of analysis parameters.
#'
#' @examples
#' library(S4Vectors)
#' de <- DataFrame(logFC = c(1.2, -0.5), P.Value = c(0.001, 0.5),
#'                 adj.P.Val = c(0.01, 0.8))
#' pb <- matrix(rpois(20, 10), nrow = 2, ncol = 10)
#' rownames(pb) <- paste0("Gene", 1:2)
#' colnames(pb) <- paste0("Donor", 1:10)
#' obj <- FDEResult(
#'     de,
#'     pb,
#'     setNames(sqrt(1:10), paste0("Donor", 1:10)),
#'     params = list(cell_type = "T_cells", condition = "group")
#' )
#' params(obj)
#'
#' @importFrom S4Vectors params
#' @export
setMethod("params", "FDEResult", function(x, ...) x@params)
