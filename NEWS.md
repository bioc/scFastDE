# scFastDE 0.99.3

## Bug fixes

* Improve overall stability.

# scFastDE 0.99.2

## Bug fixes

* Improve overall stability.

# scFastDE 0.99.1

## Bug fixes

* Removed `inst/CITATION` file lacking a DOI (Bioconductor BiocCheck WARNING).
* Updated R version dependency from 4.5.0 to 4.6.0 to match Bioconductor 3.23.

# scFastDE 0.99.0

## New features

* Initial release submitted to Bioconductor.
* `filterSparseDonors()`: remove donors below a minimum cell count threshold
  per cell type before pseudo-bulk aggregation.
* `fastPseudobulk()`: build donor-weighted pseudo-bulk profiles from a
  SingleCellExperiment using vectorised sparse matrix operations. Supports
  aggregation per donor (unpaired) or per donor × condition pair (paired).
* `fastDE()`: run vectorised differential expression across all genes
  simultaneously using limma-voom with donor cell-count weights.
  Automatically detects paired designs (same donors in multiple conditions)
  and uses a `~ 0 + condition + donor` blocking model.
* `plotDEResults()`: volcano plot with significance colouring and
  top-gene labelling.
* `FDEResult` S4 class for structured storage of DE results.
