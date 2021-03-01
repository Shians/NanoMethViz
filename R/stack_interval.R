stacked_interval_inds <- function(reads) {
    reads <- arrange(reads, start)
    reads$group <- 1:nrow(reads)
    reads$.index <- 1:nrow(reads)

    current <- 1
    merged <- numeric()

    for (i in 1:nrow(reads)) {
        if (i == current || i %in% merged) {
            next
        }

        curr_end <- reads$end[current]
        candidates <- reads %>%
            filter(
                start >= curr_end,
                !.index %in% merged
            )

        while (nrow(candidates) > 0) {
            cand_ind <- 1

            merged <- c(merged, candidates$.index[cand_ind])
            reads$group[candidates$.index[cand_ind]] <- reads$group[current]

            curr_end <- candidates$end[cand_ind]
            candidates <- candidates %>%
                filter(
                    start >= curr_end,
                    !.index %in% merged
               )
        }

        merged <- c(merged, current)
        current <- current + 1
        while (current %in% merged) {
            current <- current + 1
        }

        curr_end <- reads$end[current]
    }

    reads %>%
        pull(group) %>%
        factor()
}

# intervals <- tibble(
#     start = sample(1:2000, 100),
#     end = start + sample(5:200, 100),
#     group = 1:100
# )
#
# ggplot(intervals, aes(x = group)) +
#     geom_linerange(aes(ymin = start, ymax = end))
#
# ggplot(stacked_interval_inds(intervals), aes(x = group)) +
#     geom_linerange(aes(ymin = start, ymax = end))
#
