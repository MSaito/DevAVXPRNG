/**
 * @file dSFMTAVX2dc.cpp
 *
 * @brief The main function of parameter generator of dSFMTAVX2.
 *
 */
#include "devavxprng.h"
#include "dSFMTAVX2search.hpp"
#include "DCOptions.hpp"
#include "dSFMTAVXdc.hpp"

int main(int argc, char** argv) {
    using namespace MTToolBox;
    DCOptions opt(607);
    opt.useSR1 = false;
    opt.fixedSL1 = 19;
    opt.fixedPerm = 1;
    bool parse = opt.parse(argc, argv);
    if (!parse) {
        return -1;
    }
    return dsfmtavx_search<w256_t, dSFMTAVX2, 256>(opt, opt.count);
}
