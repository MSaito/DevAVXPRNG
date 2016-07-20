/**
 * @file SFMTAVX2dc.cpp
 *
 */

#include "devavxprng.h"
#include "SFMTAVX2search.hpp"
#include "SFMTAVXdc.hpp"

int main(int argc, char** argv) {
    using namespace MTToolBox;
    DCOptions opt(607);
    bool parse = opt.parse(argc, argv);
    if (!parse) {
        return -1;
    }
#if 0
    opt.d_p();
#endif
    return MTToolBox::search<w256_t, SFMTAVX2, 256>(opt, opt.count);
}
