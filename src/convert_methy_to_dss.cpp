#include "convert_methy_to_dss.hpp"

using namespace std;
using namespace Rcpp;

vector<string> samples_found;
map<string, bool> file_known;

pair<int, int>
find_nth_column(
        const string &str,
        const string &delim,
        int n
) {
    size_t start = 0;
    size_t end = str.find(delim);

    int n_found;
    for (n_found = 1; n_found < n && end != string::npos; n_found++) {
        start = end + 1;
        end = str.find(delim, start);
    }

    if (n_found < n) {
        throw runtime_error("end of string reached without finding");
    } else {
        return make_pair(start, end);
    }
}

entry
parse_line(const string &line) {
    pair<int, int> sample_pos = find_nth_column(line, "\t", 1);
    pair<int, int> chr_pos = find_nth_column(line, "\t", 2);
    pair<int, int> pos_pos = find_nth_column(line, "\t", 3);
    pair<int, int> stat_pos = find_nth_column(line, "\t", 5);
    entry e{
        line.substr(sample_pos.first, sample_pos.second - sample_pos.first),
        line.substr(chr_pos.first, chr_pos.second - chr_pos.first),
        stoi(line.substr(pos_pos.first, pos_pos.second - pos_pos.first)),
        stod(line.substr(stat_pos.first, stat_pos.second - stat_pos.first))
    };

    return(e);
}

void
flush_data(unordered_map<string, MethyData> const &sample_data, string const &prefix) {
    for (auto const &s_data : sample_data) {
        string const &sample_name = s_data.first;

        string out_path = make_path(prefix, sample_name + ".txt");

        ofstream out_file;

        if (file_known[sample_name]) {
            out_file.open(out_path, ios_base::out | ios_base::app);
        } else {
            out_file.open(out_path, ios_base::out | ios_base::trunc);
            out_file << "chr\tpos\ttotal\tmethylated\n";
            file_known[sample_name] = true;
            samples_found.push_back(sample_name);
        }

        for (auto const &m_count : s_data.second) {
            out_file
                << m_count.first.chr << "\t"
                << m_count.first.pos << "\t"
                << m_count.second.total << "\t"
                << m_count.second.methylated << "\n";
        }

        out_file.close();
    }
}


vector<string>
convert_methy_to_dss_cpp(
    string input,
    string output_dir
) {
    zstr::ifstream file(input, ios_base::in | ios_base::binary);

    unordered_map<string, MethyData> sample_data;
    std::string line;
    string current_chr = "";
    while (getline(file, line)) {
        entry e = parse_line(line);

        if (e.chr != current_chr) {
            current_chr = e.chr;
            stringstream ss;
            ss << "parsing " << current_chr << "...";
            LOG(ss.str());
            flush_data(sample_data, output_dir);
            sample_data.clear();
        }

        GenomicPos gpos = {e.chr, e.pos};
        sample_data[e.sample][gpos].total++;
        if (e.stat > 0) {
            sample_data[e.sample][gpos].methylated++;
        }
    }

    stringstream ss;
    ss << "samples found: ";
    for (auto const &samp : samples_found) {
        ss << samp << " ";
    }
    LOG(ss.str());

    vector<string> files_created;
    for (auto const &x : file_known) {
        files_created.push_back(x.first);
    }

    return files_created;
}
