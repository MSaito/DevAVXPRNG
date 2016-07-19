#pragma once
#ifndef SFMTAVX512FDC_H
#define SFMTAVX512FDC_H

#include "devavxprng.h"

class options {
public:
    int mexp;
    bool verbose;
    bool fixed;
    uint64_t seed;
    std::string filename;
    long count;
};

bool parse_opt(options& opt, int argc, char **argv);

int search(options& opt, int count);

#endif // SFMTAVX512FDC_H
