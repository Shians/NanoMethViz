stacked_intervals <- function(reads) {
    # record original order of reads
    reads$old_order <- seq_len(nrow(reads))

    # sort by start and assign index number
    reads <- dplyr::arrange(reads, .data$start)
    reads$index <- seq_len(nrow(reads))

    # set initial groups as one per read
    reads$group <- seq_len(nrow(reads))

    # initialize 'merged' as an empty numeric vector
    merged <- numeric()

    # loop over rows
    for (i in seq_len(nrow(reads))) {
        # skip if row is already merged
        if (i %in% merged) {
            next
        }

        # set group end to 'end' of current row
        curr_end <- reads$end[i]

        # find rows for merging
        candidates <- reads[!reads$index %in% merged & reads$start >= curr_end, ]

        # while candidates exist, merge left-most candidate into group
        while (nrow(candidates) > 0) {
            # merge first candidate into group
            merged <- c(merged, candidates$index[1])
            reads$group[candidates$index[1]] <- reads$group[i]

            # update end position of current group
            curr_end <- candidates$end[1]

            # find new candidates
            candidates <- reads[!reads$index %in% merged & reads$start >= curr_end, ]
        }

        # add current read into list of merged reads
        merged <- c(merged, i)
    }

    reads$group <- factor(reads$group)
    reads
}

stacked_interval_inds <-  function(reads) {
    stacked_intervals(reads) %>%
        dplyr::arrange(.data$old_order) %>%
        dplyr::pull(.data$group)
}
