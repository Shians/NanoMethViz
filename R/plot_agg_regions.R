#' Plot aggregate regions
#'
#' @param x the NanoMethResult object.
#' @param regions a table of regions or GRanges, or a list of such objects. The
#'   table of regions must contain chr, strand, start and end columns.
#' @param groups_feature if 'features' is a list, a vector of characters of the same
#'   length as the list containing names for each member.
#' @param flank the number of flanking bases to add to each side of each region.
#' @param stranded TRUE if negative strand features should have coordinates
#'   flipped to reflect features like transcription start sites.
#' @param span the span for loess smoothing.
#' @param palette the colour palette used for groups.
#'
#' @return a ggplot object containing the aggregate methylation trend.
#'
#' @examples
#' nmr <- load_example_nanomethresult()
#' gene_anno <- exons_to_genes(NanoMethViz::exons(nmr))
#' plot_agg_regions(nmr, gene_anno)
#'
#' @export
plot_agg_regions <- function(
    x,
    regions,
    groups_feature = NULL,
    flank = 2000,
    stranded = TRUE,
    span = 0.05,
    palette = ggplot2::scale_colour_brewer(palette = "Set1")
) {
    regions <- .normalise_regions(regions, groups_feature)
    has_groups <- !is.null(names(regions))

    # query methylation data
    methy_data <- .get_agg_methy_data(regions, x, flank, stranded, has_groups)

    # take binned means
    grid_size <- 2^10
    binned_pos_df <- .get_grid(grid_size = grid_size, rel_pos = methy_data$rel_pos)

    methy_data <- methy_data %>%
        dplyr::mutate(interval = cut(.data$rel_pos, breaks = grid_size)) %>%
        dplyr::group_by(.data$interval, .data$group, .drop = FALSE) %>%
        dplyr::summarise(methy_prob = mean(.data$methy_prob), .groups = "drop") %>%
        dplyr::left_join(binned_pos_df, by = "interval")

    # set up plot
    p <- ggplot2::ggplot() +
        ggplot2::ylim(c(0, 1)) +
        ggplot2::theme_minimal() +
        .agg_geom_smooth(methy_data, span = span, group = !is.null(groups_feature)) +
        palette

    if (flank == 0) {
        labels <- c("start", "end")
        breaks = c(0, 1)
        limits = c(0, 1)
    } else {
        breaks = c(-.33, 0, 1, 1.33)
        limits = c(-0.33, 1.33)
        kb_marker <- round(flank / 1000, 1)
        labels <- c(glue::glue("-{kb_marker}kb"), "start", "end", glue::glue("+{kb_marker}kb"))
        p <- p +
            ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey80") +
            ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "grey80")
    }

    p + ggplot2::coord_cartesian(clip = "off") +
        ggplot2::scale_x_continuous(
            name = "Relative Position",
            breaks = breaks,
            limits = limits,
            labels = labels) +
        ggplot2::ylab("Average Methylation Proportion")
}

#' Plot aggregate regions with grouped samples
#'
#' @inheritParams plot_agg_regions
#' @param palette the colour palette used for samples.
#'
#' @return a ggplot object containing the aggregate methylation trend.
#'
#' @examples
#' nmr <- load_example_nanomethresult()
#' gene_anno <- exons_to_genes(NanoMethViz::exons(nmr))
#' plot_agg_regions_sample_grouped(nmr, gene_anno)
#'
#' @export
plot_agg_regions_sample_grouped <- function(
    x,
    regions,
    flank = 2000,
    stranded = TRUE,
    span = 0.05,
    palette = ggplot2::scale_colour_brewer(palette = "Set1")
) {
    regions <- .normalise_regions(regions)

    # query methylation data per region group
    methy_data <- purrr::map(
        regions,
        function(features) {
            .get_features_with_rel_pos(
                    features,
                    methy = methy(x),
                    flank = flank,
                    stranded = stranded) %>%
                dplyr::bind_rows()
        }
    )
    methy_data <- methy_data %>%
        dplyr::bind_rows()

    # take binned means
    grid_size <- 2^10
    binned_pos_df <- .get_grid(grid_size = grid_size, rel_pos = methy_data$rel_pos)

    # take binned means
    methy_data <- methy_data %>%
        dplyr::mutate(interval = cut(.data$rel_pos, breaks = grid_size)) %>%
        dplyr::group_by(.data$sample, .data$interval) %>%
        dplyr::summarise(methy_prob = mean(.data$methy_prob), .groups = "drop")

    # add bin intervals
    methy_data <- dplyr::left_join(methy_data, binned_pos_df, by = "interval") %>%
        dplyr::mutate(rel_pos = .data$binned_pos) %>%
        dplyr::left_join(samples(x), by = "sample")

    methy_data <- dplyr::left_join(methy_data, samples(x), by = c("sample", "group"))

    kb_marker <- round(flank / 1000, 1)

    labels <- c(glue::glue("-{kb_marker}kb"), "start", "end", glue::glue("+{kb_marker}kb"))

    # set up plot
    ggplot2::ggplot() +
        ggplot2::ylim(c(0, 1)) +
        ggplot2::theme_minimal() +
        # smoothed line
        ggplot2::stat_smooth(
            aes(x = .data$rel_pos, y = .data$methy_prob, group = .data$group, col = .data$group),
            method = "loess",
            formula = "y ~ x",
            span = span,
            na.rm = TRUE,
            se = FALSE,
            data = methy_data) +
        palette +
        ggplot2::coord_cartesian(clip = "off") +
        # start and end vertical dashes
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey80") +
        ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "grey80") +
        ggplot2::scale_x_continuous(
            name = "Relative Position",
            breaks = c(-.33, 0, 1, 1.33),
            limits = c(-0.33, 1.33),
            labels = labels) +
        ggplot2::ylab("Average Methylation Proportion")
}

.split_rows <- function(x) {
    split(x, seq_len(nrow(x)))
}

# get features with flanking region and relative positions
.get_features_with_rel_pos <- function(features, methy, flank, stranded, feature_ids = NULL) {
    .scale_to_feature <- function(x, feature, stranded = stranded) {
        within <- x >= feature$start & x <= feature$end
        upstream <- x < feature$start
        downstream <- x > feature$end

        out <- numeric(length(x))
        out[within] <- (x[within] - feature$start) / (feature$end - feature$start)
        out[upstream] <- (x[upstream] - feature$start) / flank / 3
        out[downstream] <- 1 + ((x[downstream] - feature$end) / flank / 3)
        if (stranded && length(feature$strand) > 0 && feature$strand == "-") {
            out <- 1 - out
        }
        out
    }

    out <- query_methy(
        methy,
        features$chr,
        features$start - flank * 1.10,
        features$end + flank * 1.10,
        simplify = FALSE)

    # remove empty queries
    empty <- purrr::map_lgl(out, function(x) nrow(x) == 0)
    out <- out[!empty]
    features <- features[!empty, ]

    out <- purrr::map2(out, .split_rows(features),
        function(x, y) {
            x <- x %>%
                dplyr::group_by(.data$sample, .data$pos) %>%
                dplyr::summarise(methy_prob = mean(.data$statistic > 0), .groups = "drop")

            x$rel_pos <- .scale_to_feature(x$pos, y, stranded = stranded)

            x
        }
    )

    if (!is.null(feature_ids)) {
        assertthat::assert_that(length(feature_ids) == nrow(out))
        names(out) <- feature_ids
    }

    out
}

.agg_geom_smooth <- function(data, span, group) {
    if (group) {
        ggplot2::stat_smooth(
            aes(x = .data$binned_pos,
                y = .data$methy_prob,
                group = .data$group,
                col = .data$group),
            method = "loess",
            formula = "y ~ x",
            span = span,
            na.rm = TRUE,
            se = FALSE,
            data = data)
    } else {
        ggplot2::stat_smooth(
            aes(x = .data$binned_pos, y = .data$methy_prob),
            method = "loess",
            formula = "y ~ x",
            span = span,
            na.rm = TRUE,
            se = FALSE,
            data = data)
    }
}

.get_grid <- function(grid_size, rel_pos) {
    binned_pos <- seq(-1.1/3, 1 + 1.1/3, length.out = grid_size + 2)[-c(1, grid_size + 2)]
    binned_intervals <- levels(cut(rel_pos, breaks = grid_size))

    tibble::tibble(
        interval = levels(cut(rel_pos, breaks = grid_size)),
        binned_pos = seq(-1.1/3, 1 + 1.1/3, length.out = grid_size+ 2)[2:(1 + grid_size)]
    )
}

.is_df_or_granges <- function(x) {
    is.data.frame(x) || is(x, "GRanges")
}

# check that regions have required columns and are of correct class
.validate_regions <- function(regions) {
    assertthat::assert_that(.is_df_or_granges(regions) || is.list(regions))
    required <- c("chr", "strand", "start", "end")
    if (.is_df_or_granges(regions)) {
        if (!.has_required_columns(regions, required)) {
            missing <- setdiff(required, colnames(regions))
            missing <- paste(missing, collapse = ", ")
            stop(glue::glue("columns missing from 'regions': {missing}"))
        }
    } else if (is.list(regions)) {
        invalids <- !purrr::map_lgl(
            regions,
            .has_required_columns,
            required = required
        )

        if (any(invalids)) {
            invalids_id <- paste(which(invalids), collapse = ", ")
            stop(glue::glue("elements of 'region' do not have all required columns: {invalids_id}"))
        }
    }
}

.has_required_columns <- function(regions, required) {
    all(required %in% colnames(regions))
}

.normalise_regions <- function(regions, groups_feature = NULL) {
    .validate_regions(regions)

    # convert GRanges to tibble
    if (is(regions, "GRanges")) {
        regions <- tibble::as_tibble(regions) %>%
            dplyr::rename(chr = "seqnames")
    }

    # convert data.frame regions to single element list
    if (!is.null(dim(regions))) {
        regions <- list(regions)
    }

    if (!is.null(groups_feature)) {
        names(regions) <- groups_feature
    }

    regions
}

.get_agg_methy_data <- function(regions, x, flank, stranded, has_groups) {
    methy_data <- purrr::map(
        regions,
        function(features) {
            .get_features_with_rel_pos(
                    features,
                    methy = methy(x),
                    flank = flank,
                    stranded = stranded) %>%
                dplyr::bind_rows()
        }
    )

    methy_data <- methy_data %>%
        dplyr::bind_rows(.id = "group")

    # average methylation across relative coordinates
    methy_data
}
