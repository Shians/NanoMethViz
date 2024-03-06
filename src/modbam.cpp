#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

#include <Rcpp.h>
using namespace Rcpp;

// get the position of character c occurring in string x
// [[Rcpp::export]]
NumericVector
get_char_pos_cpp(CharacterVector x, CharacterVector c) {
    std::vector<int> out;
    out.reserve(x[0].size());

    // for each character in x, if c matches x[i] then i+1 is added
    // to the results as a 1-based position
    int index = 1;
    for (char ch : x[0]) {
        if (ch == c[0][0]) {
            out.push_back(index);
        }
        index++;
    }

    return Rcpp::wrap(out);
}

// tokenise a CIGAR string into CIGAR state and counts
// [[Rcpp::export]]
DataFrame
cigar_tokeniser_cpp(CharacterVector x) {
    std::vector<int> count;
    std::vector<char> state;

    int num;
    char ch;

    // tokenise and parse the string
    std::istringstream sstream(Rcpp::as<std::string>(x[0]));
    while (sstream >> num && sstream >> ch) {
        count.push_back(num);
        state.push_back(ch);
    }

    return DataFrame::create(_["state"] = state, _["count"] = count);
}

// get a mapping from read coordinates to reference coordinates
// [[Rcpp::export]]
IntegerVector
get_coord_map_cpp(std::string cigar) {
    // Tokenize cigar string
    DataFrame tokens = cigar_tokeniser_cpp(cigar);
    CharacterVector state_vec = tokens["state"];
    IntegerVector count_vec = tokens["count"];

    // Initialize seq and ref coordinates to 0
    int seq_coord = 0;
    int ref_coord = 0;

    size_t token_count = std::accumulate(count_vec.begin(), count_vec.end(), 0);
    std::vector<int> seq_map;
    seq_map.reserve(token_count);
    std::vector<int> ref_map;
    ref_map.reserve(token_count);

    // Loop over tokens and update coordinates and maps
    int num_tokens = state_vec.size();
    for (int i = 0; i < num_tokens; i++) {
        String tok = state_vec[i];
        int count = count_vec[i];
        if (tok == "M") {
            for (int j = 0; j < count; j++) {
                seq_coord++;
                ref_coord++;
                seq_map.push_back(seq_coord);
                ref_map.push_back(ref_coord);
            }
        } else if (tok == "I" || tok == "S") {
            for (int j = 0; j < count; j++) {
                seq_coord++;
                seq_map.push_back(seq_coord);
                ref_map.push_back(NA_REAL);
            }
        } else if (tok == "D" || tok == "N") {
            for (int j = 0; j < count; j++) {
                ref_coord++;
            }
        }
    }

    IntegerVector out(ref_map.begin(), ref_map.end());
    out.names() = seq_map;

    // Set names and return result
    return out;
}

// [[Rcpp::export]]
IntegerVector
get_coord_map_cpp2(std::string cigar) {
    // Tokenize cigar string
    std::vector<int> counts;
    std::vector<char> states;

    int num;
    char ch;

    // tokenise and parse the string
    std::istringstream sstream(cigar);
    while (sstream >> num && sstream >> ch) {
        counts.push_back(num);
        states.push_back(ch);
    }

    // Initialize seq and ref coordinates to 0
    int seq_coord = 0;
    int ref_coord = 0;

    size_t token_count = std::accumulate(counts.begin(), counts.end(), 0);
    std::vector<int> seq_map;
    seq_map.reserve(token_count);
    std::vector<int> ref_map;
    ref_map.reserve(token_count);

    // Loop over tokens and update coordinates and maps
    int num_tokens = states.size();
    for (int i = 0; i < num_tokens; i++) {
        int current_count = counts[i];
        switch (states[i]) {
            case 'M':
                for (int j = 0; j < current_count; ++j) {
                    ++seq_coord;
                    ++ref_coord;
                    seq_map.push_back(seq_coord);
                    ref_map.push_back(ref_coord);
                }
                break;
            case 'I':
            case 'S':
                for (int j = 0; j < current_count; ++j) {
                    ++seq_coord;
                    seq_map.push_back(seq_coord);
                    ref_map.push_back(NA_REAL);
                }
                break;
            case 'D':
            case 'N':
                ref_coord += current_count;
                break;
        }

    }

    IntegerVector out(ref_map.begin(), ref_map.end());
    out.names() = seq_map;

    // Set names and return result
    return out;
}

// [[Rcpp::export]]
DataFrame
mod_tokeniser_cpp(std::string string, std::string scores) {
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
    return DataFrame::create(
        _["mod_pos"] = mod_pos,
        _["mod_prob"] = NumericVector(mod_prob.begin(), mod_prob.end())
    );
}

// helper function for splitting string
std::vector<std::string_view>
split_string_view(std::string_view sv) {
    std::vector<std::string_view> tokens;
    tokens.reserve(sv.size() * 2/3);
    size_t start = 0;
    size_t end;

    // insert tokens into vector when split by comma or semicolon
    while ((end = sv.find_first_of(",;", start)) != std::string_view::npos) {
        tokens.emplace_back(sv.substr(start, end - start));
        start = end + 1;
    }

    // shrink to fit
    tokens.shrink_to_fit();

    return tokens;
}

// struct for storing genomic positions and mods
struct GenomicModPos {
    // constructor with no arguments
    GenomicModPos() = default;
    // constructor to reserve space for vectors
    GenomicModPos(int n_mods) {
        seq_pos.reserve(n_mods);
        pos.reserve(n_mods);
        mod_score.reserve(n_mods);
        base.reserve(n_mods);
        mod.reserve(n_mods);
    }

    std::vector<int> seq_pos;
    std::vector<int> pos;
    std::vector<double> mod_score;
    std::vector<char> base;
    std::vector<std::string> mod;
    size_t size() const {
        return pos.size();
    }
};

// complement IUPAC base
char comp_base(char b) {
    switch (b) {
        case 'A':
            return 'T';
        case 'C':
            return 'G';
        case 'G':
            return 'C';
        case 'T':
            return 'A';
        case 'U':
            return 'A';
        case 'W':
            return 'W';
        case 'S':
            return 'S';
        case 'M':
            return 'K';
        case 'K':
            return 'M';
        case 'R':
            return 'Y';
        case 'Y':
            return 'R';
        case 'B':
            return 'V';
        case 'D':
            return 'H';
        case 'H':
            return 'D';
        case 'V':
            return 'B';
        default:
            return b;
    }
}

// map positions in query to positions in genome
std::vector<int>
query_pos_to_genome_pos(const std::string& seq, int map_pos, const std::string& cigar) {
    try {
        std::vector<int> genomic_positions(seq.size(), -1);

        std::istringstream ss(cigar);
        int genomic_pos = map_pos;
        int seq_pos = 0;
        int len;
        char op;
        while (ss >> len >> op) {
            switch (op) {
                case 'M': case '=': case 'X':
                    for (int i = 0; i < len; ++i) {
                        genomic_positions.at(seq_pos++) = genomic_pos++;
                    }
                    break;
                case 'I': case 'S': case 'H':
                    seq_pos += len;
                    break;
                case 'D': case 'N':
                    genomic_pos += len;
                    break;
                default:
                    break;
            }
        }

        return genomic_positions;
    } catch (const std::exception& e) {
        throw std::runtime_error("CIGAR string does not match sequence length");
    }
}

int
get_n_mods(const std::string& str) {
    return std::count(str.begin(), str.end(), ',');
}

// enum for mode parsing mode
enum class ParseMode {
    skip_low_prob,
    skip_unknown
};

// function for calculating genomic positions
GenomicModPos
parse_bam(
    std::string const &mm_string,
    std::string const &ml_string,
    std::string const &seq,
    std::string const &cigar,
    std::string const &strand,
    int const map_pos
) {
    try {
        // get genomic positions
        std::vector<int> gpos_map = query_pos_to_genome_pos(seq, map_pos, cigar);

        // set up mod offsets
        std::istringstream ss_mm(mm_string);
        std::istringstream ss_ml(ml_string);
        std::string ml_token;

        // split mm_string into tokens for mod positions
        std::vector<std::string_view> mm_tokens = split_string_view(mm_string);

        // create output struct
        GenomicModPos output(mm_tokens.size());

        // declare variables
        size_t seq_ind = 0;
        char current_base;
        char target_base;
        // char current_strand;
        std::string current_mod;
        ParseMode parse_mode;
        int base_offset;
        int mod_prob;

        // iterate through mod positions
        for (std::string_view mm_tok : mm_tokens) {

            if (std::isdigit(mm_tok[0])) {
                // if token is a number, it is a base offset
                // store base offset and mod probability
                base_offset = std::stoi(std::string(mm_tok));
                std::getline(ss_ml, ml_token, ',');
                mod_prob = std::stoi(std::string(ml_token));
            } else {
                if (!mm_tok.empty()) {
                    // if it is mod base declaration
                    // store base, strand, and mod type

                    current_base = mm_tok[0];
                    // current_strand = mm_tok[1]; // currently unused
                    int mod_string_end = 3;
                    
                    // parse until non-alphanumeric character
                    while ( mod_string_end < mm_tok.size() && std::isalnum(mm_tok[mod_string_end]) ) {
                        mod_string_end++;
                    }
                    current_mod = std::string(mm_tok.substr(2, mod_string_end - 2));

                    // set parse mode
                    if (mod_string_end == mm_tok.size() || mm_tok[mod_string_end] == '.') {
                        parse_mode = ParseMode::skip_low_prob;
                    } else if (mm_tok[mod_string_end] == '?') {
                        parse_mode = ParseMode::skip_unknown;
                    } else {
                        throw std::runtime_error("Invalid mod parsing mode");
                    }

                    if (strand == "-") {
                        // if strand is negative, complement target base
                        target_base = comp_base(current_base);
                        seq_ind = seq.size() - 1;
                    } else {
                        target_base = current_base;
                        seq_ind = 0;
                    }
                    continue;
                } else {
                    // token is empty
                    continue;
                }
            }

            while (base_offset >= 0) {
                // if every base needs to be parsed
                if (parse_mode == ParseMode::skip_low_prob) {
                    if (seq.at(seq_ind) == target_base) {
                        --base_offset;
                        output.seq_pos.push_back(seq_ind);
                        output.pos.push_back(gpos_map[seq_ind-1]);
                        output.base.push_back(current_base);
                        output.mod.push_back(current_mod);
                        output.mod_score.push_back(0);
                    }

                // if only tagged bases need to be parsed
                } else if (parse_mode == ParseMode::skip_unknown) {
                    if (seq.at(seq_ind) == target_base) {
                        --base_offset;
                    }
                }
                if (strand == "-") {
                    seq_ind--;
                } else {
                    seq_ind++;
                }
            }

            output.seq_pos.push_back(seq_ind);
            if (strand == "-") {
                output.pos.push_back(gpos_map[seq_ind + 1]);
            } else {
                output.pos.push_back(gpos_map[seq_ind - 1]);
            }
            output.base.push_back(current_base);
            output.mod.push_back(current_mod);
            output.mod_score.push_back(mod_prob);
        }

        // filter out positions that are not in the genome
        for (size_t i = 0; i < output.pos.size(); ++i) {
            if (output.pos[i] == -1) {
                output.pos.erase(output.pos.begin() + i);
                output.seq_pos.erase(output.seq_pos.begin() + i);
                output.base.erase(output.base.begin() + i);
                output.mod.erase(output.mod.begin() + i);
                output.mod_score.erase(output.mod_score.begin() + i);
                --i;
            }
        }

        return output;
    } catch (const std::exception& e) {
        throw std::runtime_error("CIGAR string does not match sequence length");
    }

}

double
logit(double p) {
    return std::log(p / (1 - p));
}

// function for calculating genomic positions
// [[Rcpp::export]]
DataFrame
parse_bam_cpp(
    std::string const &seq,
    std::string const &cigar,
    std::string const &mm_string,
    std::string const &ml_string,
    int const map_pos,
    std::string const &strand,
    std::string mod_code
) {
    try {
        GenomicModPos output = parse_bam(
            mm_string,
            ml_string,
            seq,
            cigar,
            strand,
            map_pos
        );
        // filter out mods that don't match mod_code
        for (size_t i = 0; i < output.pos.size(); ++i) {
            if (output.mod[i] != mod_code) {
                output.pos.erase(output.pos.begin() + i);
                output.seq_pos.erase(output.seq_pos.begin() + i);
                output.base.erase(output.base.begin() + i);
                output.mod.erase(output.mod.begin() + i);
                output.mod_score.erase(output.mod_score.begin() + i);
                --i;
            }
        }

        std::vector<double> statistic(output.pos.size(), 0.0);
        for (size_t i = 0; i < output.pos.size(); ++i) {
            statistic[i] = logit(output.mod_score[i] / 255.0);
        }

        if (strand == "-") {
            if (mod_code == "m" || mod_code == "h") {
                // assume CG symmetry and shift pos by -1
                for (size_t i = 0; i < output.pos.size(); ++i) {
                    output.pos[i] -= 1;
                }
            }
        }

        return DataFrame::create(
            _["pos"] = output.pos,
            _["statistic"] = statistic
        );
    } catch (const std::exception& e) {
        // if an error occurs, return empty data frame
        Rcout << e.what() << std::endl;
        return DataFrame::create(
            _["pos"] = IntegerVector::create(),
            _["statistic"] = NumericVector::create()
        );
    }
}

// [[Rcpp::export]]
List
parse_bam_list_cpp(
    std::vector<std::string> const &seq,
    std::vector<std::string> const &cigar,
    std::vector<std::string> const &mm_string,
    std::vector<std::string> const &ml_string,
    std::vector<int> const &map_pos,
    std::vector<std::string> const &strand,
    std::string mod_code
) {
    std::vector<DataFrame> output_list;

    for (size_t i = 0; i < seq.size(); ++i) {
        // check for R interrupt every 100 reads
        if (i % 100 == 0) {
            Rcpp::checkUserInterrupt();
        }

        DataFrame output = parse_bam_cpp(
            seq[i],
            cigar[i],
            mm_string[i],
            ml_string[i],
            map_pos[i],
            strand[i],
            mod_code
        );

        output_list.push_back(output);
    }

    return Rcpp::wrap(output_list);
}
