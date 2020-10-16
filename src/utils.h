#pragma once
#include <Rcpp.h>

#include <chrono>
#include <ctime>

bool
r_message(std::string txt);

std::string
timestamp(void);

inline void
LOG(std::string const &txt) {
    r_message("[" + timestamp() + "]" + " " + txt);
}

std::string 
make_path(std::string const &prefix, std::string const &basename);
