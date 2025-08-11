#' @importFrom purrr map_lgl
assert_readable <- function(x) {
    readable <- purrr::map_lgl(x, fs::file_exists)

    if (any(!readable)) {
        unreadables <- x[!readable]
        if (length(unreadables) == 1) {
            stop(glue::glue(
                "File path '{unreadables}' does not exist.\n",
                "Please check that the file path is correct and the file exists."
            ))
        } else {
            unreadables <- paste(glue::glue("'{unreadables}'"), collapse = ", ")
            stop(glue::glue(
                "File paths {unreadables} do not exist.\n",
                "Please check that all file paths are correct and the files exist."
            ))
        }
    }
}

#' @importFrom purrr map_lgl
assert_has_index <- function(x) {
    has_index <- purrr::map_lgl(paste0(x, ".bai"), fs::file_exists)

    if (any(!has_index)) {
        no_index <- x[!has_index]
        if (length(no_index) == 1) {
            stop(glue::glue(
                "BAM file '{no_index}' is missing its index file (.bai).\n",
                "To fix this, run: samtools index {no_index}\n",
                "Or ensure the .bai file is in the same directory as the BAM file."
            ))
        } else {
            no_index_list <- paste(glue::glue("'{no_index}'"), collapse = ", ")
            stop(glue::glue(
                "BAM files {no_index_list} are missing their index files (.bai).\n",
                "To fix this, run 'samtools index' on each file:\n",
                "{paste0('samtools index ', no_index, collapse = '\n')}"
            ))
        }
    }
}

#' Assert that file paths are readable
#'
#' This function checks whether all provided file paths exist and are readable.
#' If any file paths do not exist, it throws an informative error message.
#'
#' @param x A character vector of file paths to check for existence
#'
#' @return Nothing if all files exist, otherwise throws an error
#'
#' @keywords internal
assert_valid_genomic_coords <- function(chr, start, end, allow_equal = FALSE) {
    # Check for missing values
    if (any(is.na(chr))) {
        stop("Chromosome cannot be NA. Please provide valid chromosome identifiers.")
    }

    if (any(is.na(start)) || any(is.na(end))) {
        stop("Start and end positions cannot be NA. Please provide valid genomic coordinates.")
    }

    # Convert to numeric if needed
    if (!is.numeric(start)) {
        start_numeric <- suppressWarnings(as.numeric(start))
        if (any(is.na(start_numeric))) {
            stop(glue::glue(
                "Start position must be numeric. Got: {paste(start[is.na(start_numeric)], collapse = ', ')}\n",
                "Please provide numeric start positions."
            ))
        }
        start <- start_numeric
    }

    if (!is.numeric(end)) {
        end_numeric <- suppressWarnings(as.numeric(end))
        if (any(is.na(end_numeric))) {
            stop(glue::glue(
                "End position must be numeric. Got: {paste(end[is.na(end_numeric)], collapse = ', ')}\n",
                "Please provide numeric end positions."
            ))
        }
        end <- end_numeric
    }

    # Check for positive coordinates
    if (any(start < 1)) {
        stop(glue::glue(
            "Start positions must be positive (>= 1). Got: {paste(start[start < 1], collapse = ', ')}\n",
            "Genomic coordinates are 1-based in this package."
        ))
    }

    if (any(end < 1)) {
        stop(glue::glue(
            "End positions must be positive (>= 1). Got: {paste(end[end < 1], collapse = ', ')}\n",
            "Genomic coordinates are 1-based in this package."
        ))
    }

    # Check start <= end
    comparison <- if (allow_equal) start <= end else start < end
    if (any(!comparison)) {
        invalid_idx <- which(!comparison)
        error_msg <- if (allow_equal) {
            "Start positions must be <= end positions."
        } else {
            "Start positions must be < end positions."
        }

        stop(glue::glue(
            "{error_msg}\n",
            "Invalid coordinates found:\n",
            "{paste(glue::glue('  {chr[invalid_idx]}:{start[invalid_idx]}-{end[invalid_idx]}'), collapse = '\n')}\n",
            "Please check your genomic coordinates."
        ))
    }
}

#' Validate sample annotation data frame
#' @param samples sample annotation data.frame
#' @param required_cols required column names
#' @param context description of where this is being used for error messages
#'
#' @keywords internal
assert_valid_samples <- function(samples, required_cols = c("sample", "group"), context = "sample annotation") {
    if (!is.data.frame(samples)) {
        stop(glue::glue(
            "The {context} must be a data.frame.\n",
            "Got object of class: {class(samples)[1]}\n",
            "Please provide a data.frame with columns: {paste(required_cols, collapse = ', ')}"
        ))
    }

    if (nrow(samples) == 0) {
        stop(glue::glue(
            "The {context} data.frame is empty.\n",
            "Please provide sample information with columns: {paste(required_cols, collapse = ', ')}"
        ))
    }

    # Check required columns
    assert_has_columns(samples, required_cols)

    # Check for duplicate sample names
    if ("sample" %in% required_cols && any(duplicated(samples$sample))) {
        duplicates <- samples$sample[duplicated(samples$sample)]
        stop(glue::glue(
            "Found duplicate sample names in {context}: {paste(unique(duplicates), collapse = ', ')}\n",
            "Each sample must have a unique identifier."
        ))
    }

    # Check for empty/NA sample names
    if ("sample" %in% required_cols) {
        if (any(is.na(samples$sample) | samples$sample == "")) {
            stop(glue::glue(
                "Sample names cannot be empty or NA in {context}.\n",
                "Please provide valid sample identifiers."
            ))
        }
    }
}

#' Validate exon annotation data frame
#' @param exons exon annotation data.frame
#'
#' @keywords internal
assert_valid_exons <- function(exons) {
    if (!is.data.frame(exons)) {
        stop(glue::glue(
            "Exon annotation must be a data.frame.\n",
            "Got object of class: {class(exons)[1]}\n",
            "Please provide a data.frame with required exon columns."
        ))
    }

    required_cols <- c("gene_id", "chr", "strand", "start", "end", "transcript_id", "symbol")
    assert_has_columns(exons, required_cols)

    if (nrow(exons) > 0) {
        # Validate genomic coordinates in exons
        tryCatch({
            assert_valid_genomic_coords(exons$chr, exons$start, exons$end, allow_equal = TRUE)
        }, error = function(e) {
            stop(glue::glue(
                "Invalid genomic coordinates found in exon annotation:\n",
                "{e$message}"
            ))
        })

        # Check strand values
        valid_strands <- c("+", "-", "*")
        if (any(!exons$strand %in% valid_strands)) {
            invalid_strands <- unique(exons$strand[!exons$strand %in% valid_strands])
            stop(glue::glue(
                "Invalid strand values in exon annotation: {paste(invalid_strands, collapse = ', ')}\n",
                "Valid strand values are: {paste(valid_strands, collapse = ', ')}"
            ))
        }
    }
}

#' Enhanced validation with better error messages
#' @param data_path path to methylation data file
#' @param samples sample annotation data.frame
#' @param sample_check_limit number of entries to check for sample matching
#'
#' @keywords internal
assert_valid_methy_samples <- function(data_path, samples, sample_check_limit = 1000) {
    tryCatch({
        # Check file format by reading header
        head_values <- read.table(
            gzfile(data_path),
            col.names = methy_col_names(),
            nrows = sample_check_limit,
            comment.char = ""
        )
    }, error = function(e) {
        stop(glue::glue(
            "Failed to read methylation data file '{data_path}'.\n",
            "Error: {e$message}\n",
            "Please check that the file is a valid tabix-indexed methylation file."
        ))
    })

    if (nrow(head_values) == 0) {
        stop(glue::glue(
            "Methylation data file '{data_path}' appears to be empty.\n",
            "Please check that the file contains methylation data."
        ))
    }

    # Check sample matching
    data_samples <- unique(head_values$sample)
    anno_samples <- samples$sample

    matched_samples <- intersect(data_samples, anno_samples)

    if (length(matched_samples) == 0) {
        stop(glue::glue(
            "No sample names from the data file match the sample annotation.\n",
            "Data samples (first {length(data_samples)}): {paste(data_samples, collapse = ', ')}\n",
            "Annotation samples: {paste(anno_samples, collapse = ', ')}\n",
            "Please ensure sample names match between your data and annotation."
        ))
    }

    # Warn about unmatched samples
    unmatched_data <- setdiff(data_samples, anno_samples)
    unmatched_anno <- setdiff(anno_samples, data_samples)

    if (length(unmatched_data) > 0) {
        warning(glue::glue(
            "Found {length(unmatched_data)} samples in data file not present in annotation: {paste(unmatched_data, collapse = ', ')}\n",
            "These samples will be ignored in analysis."
        ))
    }

    if (length(unmatched_anno) > 0) {
        warning(glue::glue(
            "Found {length(unmatched_anno)} samples in annotation not present in data file: {paste(unmatched_anno, collapse = ', ')}\n",
            "These samples will have no data for analysis."
        ))
    }

    message(glue::glue(
        "Successfully matched {length(matched_samples)} samples between data and annotation."
    ))
}
