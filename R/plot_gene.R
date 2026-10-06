# Check that `gene` is a single valid symbol present in the exon annotation of
# `x`, otherwise stop with a message suggesting similar or available symbols.
validate_gene_symbol <- function(x, gene) {
    if (length(gene) != 1 || is.na(gene) || gene == "") {
        stop("Gene symbol cannot be empty or NA. Please provide a valid gene symbol.")
    }

    if (nrow(exons(x)) == 0) {
        stop(glue::glue(
            "No exon annotation found in the data object.\n",
            "Gene '{gene}' cannot be plotted without exon information.\n",
            "Please add exon annotation using: exons(your_object) <- your_exons"
        ))
    }

    available_genes <- unique(stats::na.omit(exons(x)$symbol))

    if (gene %in% available_genes) {
        return(invisible(TRUE))
    }

    stop(gene_not_found_message(gene, available_genes))
}

gene_not_found_message <- function(gene, available_genes, max_similar = 10, max_available = 20) {
    prefix <- tolower(substr(gene, 1, 3))
    similar_genes <- available_genes[startsWith(tolower(available_genes), prefix)]

    error_msg <- c(
        glue::glue("Gene '{gene}' not found in exon annotation."),
        "Please check the gene symbol spelling."
    )

    if (length(similar_genes) > 0) {
        shown_genes <- paste(utils::head(similar_genes, max_similar), collapse = ", ")
        if (length(similar_genes) > max_similar) {
            shown_genes <- glue::glue("{shown_genes}, ... ({length(similar_genes) - max_similar} more)")
        }
        error_msg <- c(error_msg, glue::glue("Similar genes found: {shown_genes}"))
    }

    if (length(available_genes) <= max_available) {
        error_msg <- c(
            error_msg,
            glue::glue("Available genes: {paste(available_genes, collapse = ', ')}")
        )
    } else {
        error_msg <- c(
            error_msg,
            glue::glue("Use unique(exons(your_object)$symbol) to see all {length(available_genes)} available genes.")
        )
    }

    paste(error_msg, collapse = "\n")
}

plot_gene_impl <- function(
        x,
        gene,
        window_prop = 0.3,
        anno_regions = NULL,
        binary_threshold = NULL,
        avg_method = c("mean", "median"),
        spaghetti = FALSE,
        heatmap = TRUE,
        heatmap_subsample = 50,
        smoothing_window = 2000,
        gene_anno = TRUE,
        palette = ggplot2::scale_colour_brewer(palette = "Set1"),
        line_size = 1,
        mod_scale = c(0, 1),
        span = NULL
) {
    if (!missing(span)) {
        warning("the 'span' argument has been deprecated, please use 'smoothing_window' instead")
    }

    avg_method <- match.arg(avg_method)

    validate_gene_symbol(x, gene)

    if (length(window_prop) == 1) {
        # convert to two sided window_prop
        window_prop <- c(window_prop, window_prop)
    }

    exons_anno <- query_exons_symbol(x, symbol = gene)

    feature_chr <- unique(exons_anno$chr)
    feature_start <- min(exons_anno$start)
    feature_end <- max(exons_anno$end)

    plot_region(
        x = x,
        chr = feature_chr,
        start = feature_start,
        end = feature_end,
        anno_regions = anno_regions,
        binary_threshold = binary_threshold,
        avg_method = avg_method,
        spaghetti = spaghetti,
        heatmap = heatmap,
        heatmap_subsample = heatmap_subsample,
        smoothing_window = smoothing_window,
        gene_anno = gene_anno,
        window_prop = window_prop,
        palette = palette,
        line_size = line_size,
        mod_scale = mod_scale
    )
}

#' @rdname plot_gene
#'
#' @inheritParams plot_region
#'
#' @details
#' This function plots the methylation data for a given gene. The main trendline plot shows the average methylation
#' probability across the gene. The heatmap plot shows the methylation probability for each read across the gene. The
#' gene annotation plot shows the exons of the gene. In the heatmap, each row represents one or more non-overlapping reads
#' where the coloured segments represent the methylation probability at each position. Data along a read is connected by
#' a grey line. The gene annotation plot shows the isoforms and exons of the gene, with arrows indicating the direction
#' of transcription.
#'
#' Since V3.0.0 NanoMethViz has changed the smoothing strategy from a loess smoothing to a weighted moving average. This
#' is because the loess smoothing was too computationally expensive for large datasets and had a span parameter that was
#' difficult to tune. The new smoothing strategy is controlled by the smoothing_window argument.
#'
#' @examples
#' nmr <- load_example_nanomethresult()
#' plot_gene(nmr, "Peg3")
#'
#' @export
setMethod("plot_gene", signature(x = "NanoMethResult", gene = "character"),
    plot_gene_impl
)

#' @describeIn plot_gene S4 method for ModBamResult
setMethod("plot_gene", signature(x = "ModBamResult", gene = "character"),
    plot_gene_impl
)
