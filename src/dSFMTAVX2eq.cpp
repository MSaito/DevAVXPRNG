#include "devavxprng.h"
#include "dSFMTAVX2search.hpp"
#include "EQOptions.hpp"
#include "dSFMTAVXeq.hpp"

using namespace MTToolBox;

int main(int argc, char * argv[])
{
    EQOptions<dSFMTAVX2_param> opt;
    if (!opt.parse(argc, argv)) {
        return -1;
    }
    int r
        = dsfmtavx_equidistribution<w256_t, dSFMTAVX2,
                                    dSFMTAVX2_param, 256>(opt);
    return r;
}
