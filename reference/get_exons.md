# Get exon annotations

Helper functions are provided for obtaining exon annotations from
relevant TxDb packages on Bioconductor for the construction of
NanoMethResult objects.

## Usage

``` r
get_cgi_mm10()

get_cgi_grcm39()

get_cgi_t2t()

get_cgi_hg19()

get_cgi_hg38()

get_exons_mm10()

get_exons_grcm39()

get_exons_hg19()

get_exons_hg38()

get_exons_t2t()
```

## Value

tibble (data.frame) object containing exon annotation.

## Examples

``` r
cgi_mm10 <- get_cgi_mm10()

cgi_GRCm39 <- get_cgi_grcm39()

cgi_t2t <- get_cgi_t2t()

cgi_hg19 <- get_cgi_hg19()

cgi_hg38 <- get_cgi_hg38()

mm10_exons <- get_exons_mm10()

grcm39_exons <- get_exons_grcm39()

hg19_exons <- get_exons_hg19()

hg38_exons <- get_exons_hg38()

t2t_exons <- get_exons_t2t()
```
