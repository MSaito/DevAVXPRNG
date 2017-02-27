#pragma once
#ifndef DSFMTAVXEQ_HPP
#define DSFMTAVXEQ_HPP
/**
 * @file dSFMTAVXeq.hpp
 */

#include "devavxprng.h"
#include "AlgorithmDSFMTEquidistribution.hpp"
//#include <MTToolBox/AlgorithmEquidistribution.hpp>
//#include <MTToolBox/AlgorithmReducibleRecursionSearch.hpp>
//#include <MTToolBox/period.hpp>
//#include <NTL/GF2X.h>
#include "Annihilate.hpp"
#include "EQOptions.hpp"

namespace MTToolBox {

    template<typename U, typename G, typename P, int bitWidth>
    int dsfmtavx_equidistribution(EQOptions<P> opt)
    {
        using namespace std;
        G sf(opt.params);
        U wseed;
        wseed.u64[0] = opt.seed;
        sf.seed(wseed);
        Annihilate<G, U> annihilate;
        if (!annihilate.anni(sf)) {
            return -1;
        }
        int delta52 = 0;
        int veq52[52];
        DSFMTInfo info;
        info.bitSize = bitWidth; // IMPORTANT
        info.elementNo = bitWidth / 64;
        delta52 = calc_dSFMT_equidistribution<U, G>
            (sf, veq52, 52, info, opt.params.mexp);
        cout << sf.getParamString();
        cout << dec << delta52 << endl;
        if (opt.verbose) {
            cout << "52bit dimension of equidistribution at v-bit accuracy k(v)"
                 << endl;
            for (int j = 0; j < 52; j++) {
                cout << "k(" << dec << (j + 1) << ") = " << dec << veq52[j];
                cout << "\td(" << dec << (j + 1) << ") = " << dec
                     << (opt.params.mexp / (j + 1) - veq52[j]) << endl;
            }
        }
        return 0;
    }
}

#endif // DSFMTAVXEQ_HPP
