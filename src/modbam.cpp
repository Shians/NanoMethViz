#include <string>
#include <sstream>

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector get_char_pos_cpp(CharacterVector x, CharacterVector c) {

    NumericVector out;
    for (auto i = 0; i < x[0].size(); i++) {
        if (x[0][i] == c[0][0]) {
            out.push_back(i+1);
        }
    }

    return out;
}

// [[Rcpp::export]]
DataFrame cigar_tokeniser_cpp(CharacterVector x) {
    std::istringstream sstream(Rcpp::as<std::string>(x[0]));

    int num;
    char ch;

    NumericVector count;
    CharacterVector state;

    while (sstream >> num && sstream >> ch) {
        count.push_back(num);
        state.push_back(ch);
    }

    return DataFrame::create(_["state"] = state, _["count"] = count);
}
