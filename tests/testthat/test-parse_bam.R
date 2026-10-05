test_that("parse_bam_cpp handles multi-code MM tags", {
    # C at read positions 2, 6, 10, 14 (1-based)
    seq <- "ACGTCACGTCACGTCACGTC"
    cigar <- "20M"
    mm <- "C+mh,0,2;"
    # one ML value per declared code per listed position (m, h)
    ml <- "200,100,180,90"
    logit <- function(p) log(p / (1 - p))

    # m is the first declared code, so takes the first ML value of each pair
    out_m <- parse_bam_cpp(seq, cigar, mm, ml, 100, "+", "m")
    expect_equal(out_m$pos, c(101L, 109L))
    expect_equal(out_m$statistic, c(logit(200 / 255), logit(180 / 255)))

    # h is the second declared code, so takes the second ML value of each pair
    out_h <- parse_bam_cpp(seq, cigar, mm, ml, 100, "+", "h")
    expect_equal(out_h$pos, c(101L, 109L))
    expect_equal(out_h$statistic, c(logit(100 / 255), logit(90 / 255)))

    # mod codes not present in the declaration return no rows
    out_a <- parse_bam_cpp(seq, cigar, mm, ml, 100, "+", "a")
    expect_equal(nrow(out_a), 0)
})

test_that("parse_bam_cpp keeps ML stream in sync across multi-code groups", {
    seq <- "ACGTCACGTCACGTCACGTC"
    cigar <- "20M"
    # first group declares two codes, second group a single code
    mm <- "C+mh,0;G+m,3;"
    ml <- "200,100,250"
    logit <- function(p) log(p / (1 - p))

    # only one ML value is consumed for the single-code G group
    out_m <- parse_bam_cpp(seq, cigar, mm, ml, 100, "+", "m")
    expect_equal(out_m$pos, c(101L, 117L))
    expect_equal(out_m$statistic, c(logit(200 / 255), logit(250 / 255)))

    # group dropped when mod code absent, later groups still parsed correctly
    out_m2 <- parse_bam_cpp(seq, cigar, "C+ah,0;G+m,3;", "200,100,250", 100, "+", "m")
    expect_equal(out_m2$pos, 117L)
    expect_equal(out_m2$statistic, logit(250 / 255))
})

test_that("parse_bam_cpp handles multi-code MM tags with '.' flag", {
    seq <- "ACGTCACGTCACGTCACGTC"
    cigar <- "20M"
    mm <- "C+mh.,0,1;"
    ml <- "200,100,180,90"
    logit <- function(p) log(p / (1 - p))

    out_m <- parse_bam_cpp(seq, cigar, mm, ml, 100, "+", "m")
    # unlisted C at read position 5 skipped, recorded with probability 0
    expect_equal(out_m$pos, c(101L, 104L, 106L))
    expect_equal(out_m$statistic, c(logit(200 / 255), logit(0), logit(180 / 255)))

    out_h <- parse_bam_cpp(seq, cigar, mm, ml, 100, "+", "h")
    expect_equal(out_h$pos, c(101L, 104L, 106L))
    expect_equal(out_h$statistic, c(logit(100 / 255), logit(0), logit(90 / 255)))
})
