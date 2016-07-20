#pragma once
#ifndef DSFMTAVX512FDC_H
#define DSFMTAVX512FDC_H

#include "devavxprng.h"

class options {
public:
    int mexp;
    bool verbose;
    bool fixed;
    int fixedSL1;
    uint64_t seed;
    std::string filename;
    long count;
};

bool parse_opt(options& opt, int argc, char **argv);


int search(options& opt, int count);
#endif // DSFMTAVX512FDC_H
