#pragma once
#include <Rcpp.h>

#include "utils.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <exception>
#include <map>
#include <unordered_map>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

struct entry {
    std::string sample;
    std::string chr;
    int pos;
    double stat;
};

struct MethyCount {
    uint N=0;
    uint X=0;
};

struct GenomicPos {
    std::string chr;
    int pos;
};

bool operator < (GenomicPos const &a, GenomicPos const &b) {
    return (a.pos < b.pos  || a.chr < b.chr) && (a.chr <= b.chr);
}

typedef std::map<GenomicPos, MethyCount> MethyData;

// [[Rcpp::export]]
std::vector<std::string> 
convert_methy_to_dss_cpp(
    std::string input,
    std::string output_dir
);
