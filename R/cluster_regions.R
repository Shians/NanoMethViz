#' Cluster regions by K-means
#'
#' Cluster regions by k-means based on their methylation profiles. In order to
#' cluster using k-means the methylation profile of each region is interpolated
#' and sampled at fixed points. The first 10 principal components are used for
#' the k-means clustering. The clustering is best behaved in regions of similar
#' width and CpG density.
#'
#' @param x the NanoMethResult object.
#' @param regions a table of regions containing at least columns chr, strand,
#'   start and end.
#' @param centers number of centers for k-means, identical to the number of
#'   output clusters.
#' @param grid_method the method for generating the sampling grid. The default
#'   option "density" attempts to create a grid with similar density as the
#'   data, "uniform" creates a grid of uniform density.
#'
#' @return the table of regions given by the 'regions' argument with the column
#'   'cluster' added.
#'
#' @importFrom stats quantile approxfun prcomp kmeans sd
#' @importFrom graphics hist
#' @export
#'
#' @examples
#' nmr <- load_example_nanomethresult()
#' gene_anno <- exons_to_genes(NanoMethViz::exons(nmr))
#' # uniform grid due to low number of input features
#' gene_anno_clustered <- cluster_regions(nmr, gene_anno, centers = 2, grid_method = "uniform")
#' plot_agg_regions(nmr, gene_anno_clustered, group_col = "cluster")
cluster_regions <- function(x, regions, centers = 2, grid_method = c("density", "uniform")) {
    grid_method <- match.arg(grid_method)

    methy_df <- tibble(
        regions,
        methy = query_methy(x, regions$chr, regions$start, regions$end, simplify = FALSE)
    )

    # rescale position values
    for (i in 1:nrow(methy_df)) {
        start <- methy_df$start[i]
        end <- methy_df$end[i]
        methy_df$methy[[i]] <- if (methy_df$strand[i] == "-") {
            # re-sort by position due to flip
            methy_df$methy[[i]] %>%
                mutate(pos = 1 - ((.data$pos - start) / (end - start))) %>%
                arrange(.data$pos)
        } else {
            # include * case
            methy_df$methy[[i]] %>%
                mutate(pos = (.data$pos - start) / (end - start))
        }
    }

    # take grid size as twice the 90% percentile
    grid_length <- 2 * floor(quantile(purrr::map_int(methy_df$methy, ~length(unique(.x$pos))), 0.9))

    if (grid_method == "uniform") {
        grid <- seq(0, 1, length.out = grid_length)
    } else if (grid_method == "density") {
        unique_pos <- unlist(purrr::map(methy_df$methy, ~unique(.x$pos)))
        hist_out <- hist(unique_pos, plot = FALSE, breaks = grid_length)
        epdf <- approxfun(hist_out$mids, hist_out$density/sum(hist_out$density), rule = 2)
        epdf_vals <- epdf(seq(0, 1, length.out = grid_length))

        ecdf_vals <- cumsum(epdf_vals)/max(cumsum(epdf_vals))
        inv_cdf <- approxfun(ecdf_vals, seq(0, 1, length.out = grid_length), rule = 2)

        grid <- inv_cdf(seq(0, 1, length.out = grid_length))
    }

    sample_methy_grid <- function(x, grid) {
        summaried_methy <- x %>%
            group_by(.data$pos) %>%
            summarise(methy_prob = median(e1071::sigmoid(.data$statistic))) %>%
            ungroup()

        smoothed_f <- approxfun(
            summaried_methy$pos,
            summaried_methy$methy_prob,
            rule = 2
        )

        smoothed_f(grid)
    }

    methy_matrix <- map(methy_df$methy, sample_methy_grid, grid = grid) %>%
        unlist() %>%
        matrix(nrow = length(methy_df$methy), byrow = TRUE)

    pca <- prcomp(methy_matrix)
    pc10 <- as.matrix(pca$x[, 1:min(ncol(pca$x), 10)])
    km_pca <- kmeans(methy_matrix, centers = centers)

    row_means_split <- function(x, f) {
        out <- list()
        unique_f <- unique(f)

        for (i in seq_along(unique_f)) {
            out[[i]] <- colMeans(x[which(f == i), , drop = FALSE])
        }

        out
    }

    cluster_means <- row_means_split(methy_matrix, km_pca$cluster)

    out <- regions %>%
        mutate(cluster = km_pca$cluster)

    out$cluster_dev <- numeric(nrow(out))
    for (i in 1:nrow(out)) {
        out$cluster_dev[i] <- sum((methy_matrix[i, ] - cluster_means[[out$cluster[i]]])^2) / ncol(methy_matrix)
    }

    out %>%
        mutate(cluster = factor(.data$cluster))
}
