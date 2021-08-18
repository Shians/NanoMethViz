#include "utils.h"

using namespace std;
using namespace Rcpp;

bool r_message(string txt) {
    Function msg("message");
    msg(txt);
    return true;
}

string timestamp() {
    time_t t = time(nullptr);

    char time_str[100];
    strftime(time_str, sizeof(time_str), "%Y-%m-%d %H:%M:%S", localtime(&t));
    return string(time_str);
}

string make_path(string const &prefix, string const &basename) {
    if (prefix != "") {
        return prefix + "/" + basename;
    } else {
        return basename;
    }
}
