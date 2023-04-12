#include <string>
#include <sstream>
#include <algorithm>
#include <ranges>
#include <string_view>

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

// [[Rcpp::export]]
DataFrame cigar_tokeniser_cpp2(CharacterVector x) {
    std::istringstream sstream(Rcpp::as<std::string>(x[0]));

    int num;
    char ch;

    std::vector<int> count;
    std::vector<char> state;

    while (sstream >> num && sstream >> ch) {
        count.push_back(num);
        state.push_back(ch);
    }

    return DataFrame::create(_["state"] = state, _["count"] = count);
}


// [[Rcpp::export]]
DataFrame mod_tokeniser_cpp(std::string string, std::string scores) {
    // Remove semicolons from string
    string.erase(std::remove(string.begin(), string.end(), ';'), string.end());

    // Split mod_string and convert to numeric
    std::vector<std::string> mod_string_split;
    std::stringstream mod_string_stream(string);
    std::string mod_string_segment;
    while (std::getline(mod_string_stream, mod_string_segment, ',')) {
        mod_string_split.push_back(mod_string_segment);
    }
    // Remove leading annotation
    mod_string_split.erase(mod_string_split.begin());

    std::vector<int> mod_offsets(mod_string_split.size());
    std::transform(mod_string_split.begin(), mod_string_split.end(), mod_offsets.begin(),
        [](const std::string& s){ return std::stoi(s); }
    );

    // Calculate mod_pos
    std::vector<int> mod_pos(mod_offsets.size());
    std::transform(mod_offsets.begin(), mod_offsets.end(), mod_offsets.begin(),
        [](int i){ return i + 1; }
    );
    std::partial_sum(mod_offsets.begin(), mod_offsets.end(), mod_pos.begin());

    // Split mod_scores and convert to numeric
    std::vector<std::string> mod_scores_split;
    std::stringstream mod_scores_stream(scores);
    std::string mod_score_segment;
    while (std::getline(mod_scores_stream, mod_score_segment, ',')) {
        mod_scores_split.push_back(mod_score_segment);
    }
    std::vector<double> mod_prob(mod_scores_split.size());
    std::transform(mod_scores_split.begin(), mod_scores_split.end(), mod_prob.begin(),
        [](const std::string& s){ return std::stod(s) / 255.0; }
    );

    // Create data frame and return
    return DataFrame::create(_["mod_pos"] = mod_pos, _["mod_prob"] = NumericVector(mod_prob.begin(), mod_prob.end()));
}

