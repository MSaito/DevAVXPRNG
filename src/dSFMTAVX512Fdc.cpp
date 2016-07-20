/**
 * @file dSFMTAVX512Fdc.cpp
 */

#include "devavxprng.h"
#include "dSFMTAVX512Fsearch.hpp"
#include "DCOptions.hpp"
#include "dSFMTAVXdc.hpp"

using namespace MTToolBox;

int main(int argc, char** argv) {
    DCOptions opt(1279);
    bool parse = opt.parse(argc, argv);
    if (!parse) {
        return -1;
    }
    return search<w512_t, dSFMTAVX512F, 512>(opt, opt.count);
}
