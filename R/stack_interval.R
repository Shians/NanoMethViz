stacked_intervals <- function(reads) {
    reads$old_order <- 1:nrow(reads)
    reads <- arrange(reads, start)
    reads$index <- 1:nrow(reads)
    reads$group <- 1:nrow(reads)

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
                !index %in% merged
            )

        while (nrow(candidates) > 0) {
            cand_ind <- 1

            merged <- c(merged, candidates$index[cand_ind])
            reads$group[candidates$index[cand_ind]] <- reads$group[current]

            curr_end <- candidates$end[cand_ind]
            candidates <- candidates %>%
                filter(
                    start >= curr_end,
                    !index %in% merged
               )
        }

        merged <- c(merged, current)
        current <- current + 1
        while (current %in% merged) {
            current <- current + 1
        }

        curr_end <- reads$end[current]
    }

    reads$group <- factor(reads$group)
    reads
}

stacked_interval_inds <-  function(reads) {
    stacked_intervals(reads) %>%
        arrange(old_order) %>%
        pull(group)
}
