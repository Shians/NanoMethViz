#include <string>
#include <Rcpp.h>

using namespace std;
using namespace Rcpp;

// function to count "CG" dinucleotides in a DNA string
// [[Rcpp::export]]
int count_cg_cpp(std::string str) {
    int count = 0;
    int length = str.length();

    // loop through string
    for (int i = 0; i < length - 1; i++) {
        // if "CG" found, increment count
        if (str[i] == 'C' && str[i + 1] == 'G') {
            count++;
        }
    }

    return count;
}
