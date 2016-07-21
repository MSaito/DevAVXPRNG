#include "devavxprng.h"
#include "dSFMTAVX512Fsearch.hpp"
#include "EQOptions.hpp"
#include "dSFMTAVXeq.hpp"

using namespace MTToolBox;

int main(int argc, char * argv[])
{
    EQOptions<dSFMTAVX512F_param> opt;
    if (!opt.parse(argc, argv)) {
        return -1;
    }
    int r
        = dsfmtavx_equidistribution<w512_t, dSFMTAVX512F,
                                   dSFMTAVX512F_param, 512>(opt);
    return r;
}
