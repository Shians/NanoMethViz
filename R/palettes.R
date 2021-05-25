colour_blind_palettes <- list(
    # http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
    "rcookbook1" = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
    "rcookbook2" = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),

    "rcookbook_remix" = c("#009E73", "#D55E00", "#0072B2", "#E69F00", "#56B4E9", "#CC79A7", "#F0E442", "#999999"),

    # https://www.ibm.com/design/language/resources/color-library
    "ibm" = c("#648FFF", "#785EF0", "#DC267F", "#FE6100", "#FFB000"),

    # https://www.nature.com/articles/nmeth.1618
    "wong" = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),

    # https://personal.sron.nl/~pault/
    "tol" = c("#332288", "#117733", "#44AA99", "#88CCEE", "#DDCC77", "#CC6677", "#AA4499", "#882255"),
    "tol_remixed" = c("#332288", "#882255", "#117733", "#AA4499", "#44AA99", "#CC6677", "#88CCEE", "#DDCC77")
)

scale_cb_d <- function(palette) {
    ggplot2::scale_colour_manual(values = colour_blind_palettes[[palette]])
}
