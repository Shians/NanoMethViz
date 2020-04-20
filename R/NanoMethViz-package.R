#' @details The main plotting functions in this package are [plot_gene()] and
#'   [plot_region()].
#' @docType package
#'
#' @importFrom magrittr %>%
#' @importFrom methods new .valueClassTest
#' @importFrom rlang .data
#' @importFrom ggplot2 aes geom_rect geom_segment geom_text ggplot theme_void
#'   xlim ylim ggplot_build rel unit
#' @importFrom dplyr filter group_by inner_join mutate n select summarise
#' @importFrom tidyr unnest
#' @importFrom glue glue
#' @importFrom assertthat assert_that is.readable is.writeable
#' @importFrom stringr str_extract
#' @importFrom readr cols col_character col_integer col_logical col_double
#' @importFrom tibble tibble as_tibble
#' @import patchwork
#' @import assertthat
#' @keywords internal
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
