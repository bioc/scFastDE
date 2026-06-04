.validate_sce_coldata <- function(sce, columns) {
    if (!is(sce, "SingleCellExperiment"))
        stop("'sce' must be a SingleCellExperiment object.")

    missing_cols <- vapply(columns, function(col) {
        is.null(col) || length(col) != 1L || !col %in% names(colData(sce))
    }, logical(1))

    if (any(missing_cols)) {
        arg_name <- names(columns)[which(missing_cols)[1L]]
        col <- columns[[arg_name]]
        if (is.null(col) || length(col) != 1L)
            col <- ""
        stop("'", arg_name, "' column '", col,
             "' not found in colData(sce).")
    }

    invisible(TRUE)
}
