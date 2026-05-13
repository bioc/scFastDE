#' @title Volcano Plot of Differential Expression Results
#'
#' @description
#' Produces a volcano plot from a \code{\link{FDEResult}} object, with
#' genes coloured by significance and optionally labelled. Points are
#' coloured by whether they exceed the \code{lfc_thresh} and
#' \code{fdr_thresh} thresholds.
#'
#' @param result A \code{\link{FDEResult}} object from \code{\link{fastDE}}.
#' @param fdr_thresh A \code{numeric(1)} FDR threshold for significance.
#'   Default: \code{0.05}.
#' @param lfc_thresh A \code{numeric(1)} absolute log2 fold-change
#'   threshold. Default: \code{1}.
#' @param top_n An \code{integer(1)} number of top significant genes to
#'   label by name. Default: \code{10}.
#' @param point_size A \code{numeric(1)} point size. Default: \code{1}.
#' @param point_alpha A \code{numeric(1)} point transparency. Default:
#'   \code{0.6}.
#' @param colours A named \code{character} vector of colours for
#'   \code{"up"}, \code{"down"}, and \code{"ns"} (not significant).
#'   Default: \code{c(up = "#E24B4A", down = "#378ADD", ns = "#888780")}.
#'
#' @return A \code{ggplot2} object.
#'
#' @examples
#' library(SingleCellExperiment)
#'
#' set.seed(42)
#' n_genes <- 200
#' counts <- matrix(rpois(n_genes * 60, 8), nrow = n_genes, ncol = 60)
#' rownames(counts) <- paste0("Gene", seq_len(n_genes))
#' colnames(counts) <- paste0("Cell", seq_len(60))
#' counts[1:10, 31:60] <- counts[1:10, 31:60] * 3L
#'
#' sce <- SingleCellExperiment(assays = list(counts = counts))
#' sce$donor     <- rep(paste0("D", 1:6), each = 10)
#' sce$cell_type <- "Tcell"
#' sce$condition <- rep(c("ctrl", "treat"), each = 30)
#'
#' result <- fastDE(sce, donor = "donor", cell_type = "cell_type",
#'                  condition = "condition", target_type = "Tcell",
#'                  min_cells = 5)
#' plotDEResults(result)
#'
#' @seealso \code{\link{fastDE}}, \code{\link{FDEResult}}
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_vline geom_hline scale_color_manual labs theme_bw theme element_text element_blank element_rect element_line geom_text
#' @importFrom methods is
#' @importFrom utils head
#' @export
plotDEResults <- function(result,
                           fdr_thresh  = 0.05,
                           lfc_thresh  = 1,
                           top_n       = 10L,
                           point_size  = 1,
                           point_alpha = 0.6,
                           colours     = c(up   = "#E24B4A",
                                           down = "#378ADD",
                                           ns   = "#888780")) {
    if (!is(result, "FDEResult"))
        stop("'result' must be a FDEResult object.")

    dt <- as.data.frame(deTable(result))

    # Significance classification
    dt$sig <- "ns"
    dt$sig[dt$adj.P.Val < fdr_thresh & dt$logFC >  lfc_thresh] <- "up"
    dt$sig[dt$adj.P.Val < fdr_thresh & dt$logFC < -lfc_thresh] <- "down"
    dt$sig <- factor(dt$sig, levels = c("up", "down", "ns"))

    dt$neg_log10_p <- -log10(dt$P.Value + 1e-300)

    # Top genes to label
    sig_genes <- dt[dt$sig != "ns", ]
    sig_genes <- sig_genes[order(sig_genes$adj.P.Val), ]
    label_genes <- head(sig_genes, top_n)

    p <- ggplot(dt, aes(x = .data[["logFC"]],
                        y = .data[["neg_log10_p"]],
                        colour = .data[["sig"]])) +
        geom_point(size = point_size, alpha = point_alpha) +
        geom_vline(xintercept = c(-lfc_thresh, lfc_thresh),
                   linetype = "dashed", colour = "grey50",
                   linewidth = 0.4) +
        geom_hline(yintercept = -log10(fdr_thresh),
                   linetype = "dashed", colour = "grey50",
                   linewidth = 0.4) +
        scale_color_manual(
            values = colours,
            labels = c(
                up   = sprintf("Up (%d)",   sum(dt$sig == "up")),
                down = sprintf("Down (%d)", sum(dt$sig == "down")),
                ns   = sprintf("NS (%d)",   sum(dt$sig == "ns"))
            )
        ) +
        labs(
            x      = expression(log[2]~"fold change"),
            y      = expression(-log[10]~"p-value"),
            colour = "Regulation",
            title  = paste0("DE: ", result@params$cell_type,
                            " | ", result@params$contrast)
        ) +
        theme_bw(base_size = 11) +
        theme(
            panel.grid.minor = element_blank(),
            axis.text        = element_text(size = 10),
            plot.title       = element_text(size = 12)
        )

    # Label top genes if any are significant
    if (nrow(label_genes) > 0) {
        p <- p + geom_text(
            data        = label_genes,
            aes(label   = .data[["gene"]]),
            size        = 3,
            vjust       = -0.5,
            show.legend = FALSE
        )
    }

    p
}
