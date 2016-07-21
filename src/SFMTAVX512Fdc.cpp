/**
 * @file SFMTAVX512Fdc.cpp
 */
#include "devavxprng.h"
#include "DCOptions.hpp"
#include "SFMTAVX512Fsearch.hpp"
#include "SFMTAVXdc.hpp"

using namespace MTToolBox;

int main(int argc, char** argv) {
    DCOptions opt(1279);
    opt.useSR1 = true;
    opt.fixedSL1 = 19;
    opt.fixedSR1 = 4;
    opt.fixedPerm = 1;
    bool parse = opt.parse(argc, argv);
    if (!parse) {
        return -1;
    }
    return sfmtavx_search<w512_t, SFMTAVX512F, 512>(opt, opt.count);
}
