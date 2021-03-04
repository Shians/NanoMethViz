stacked_intervals <- function(reads) {
    reads$old_order <- 1:nrow(reads)
    reads <- dplyr::arrange(reads, start)
    reads$index <- 1:nrow(reads)
    reads$group <- 1:nrow(reads)

    current <- 1
    merged <- numeric()

    for (i in 1:nrow(reads)) {
        if (i == current || i %in% merged) {
            # skip if already merged
            next
        }

        # update current group end
        curr_end <- reads$end[current]

        # list reads after current end
        candidates <- reads %>%
            dplyr::filter(
                .data$start >= curr_end,
                !.data$index %in% merged
            )

        while (nrow(candidates) > 0) {
            cand_ind <- 1

            # merge in candidate
            merged <- c(merged, candidates$index[cand_ind])
            reads$group[candidates$index[cand_ind]] <- reads$group[current]

            # update current group end
            curr_end <- candidates$end[cand_ind]

            # update candidates
            candidates <- candidates %>%
                dplyr::filter(
                    .data$start >= curr_end,
                    !.data$index %in% merged
               )
        }

        # consider current read merged
        merged <- c(merged, current)

        current <- current + 1
        while (current %in% merged) {
            # skip forward if next read is already merged
            current <- current + 1
        }

        # update end for next read group
        curr_end <- reads$end[current]
    }

    reads$group <- factor(reads$group)
    reads
}

stacked_interval_inds <-  function(reads) {
    stacked_intervals(reads) %>%
        dplyr::arrange(.data$old_order) %>%
        dplyr::pull(.data$group)
}
