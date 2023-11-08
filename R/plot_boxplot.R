#' Plot boxplots for regions
#'
#' @param x the NanoMethResult object.
#' @param regions a table of regions containing at least columns chr, strand,
#'   start and end. Any additional columns can be used for grouping.
#' @param binary_threshold the modification probability such that calls with
#'   modification probability above the threshold are considered methylated, and
#'   those with probability equal or below are considered unmethylated.
#' @param group_col the column to group aggregated trends by. This column can
#'   be in from the regions table or samples(x).
#' @param palette the ggplot colour palette used for groups.
#'
#' @return a ggplot object containing the methylation boxplots.
#'
#' @examples
#' nmr <- load_example_nanomethresult()
#' gene_anno <- exons_to_genes(NanoMethViz::exons(nmr))
#' plot_boxplot(nmr, gene_anno)
#' plot_boxplot(nmr, gene_anno, group_col = "sample")
#' plot_boxplot(nmr, gene_anno, group_col = "group")
#'
#' @export
# plot_boxplot <- function(
#     x,
#     regions,
#     binary_threshold = 0.5,
#     group_col = NULL,
#     palette = ggplot2::scale_colour_brewer(palette = "Set1")
# ) {
#     if (!is.null(group_col)) {
#         avail_columns <- c(colnames(samples(x)), colnames(regions))
#         assertthat::assert_that(
#             group_col %in% avail_columns,
#             msg = glue::glue("'{group_col}' could not be found in columns of 'regions' or samples(x)")
#         )
#     }
#
#     # grouped regions crashes downstream operations
#     regions <- ungroup(regions)
#
#
# }
