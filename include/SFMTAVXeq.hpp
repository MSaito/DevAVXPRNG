#pragma once
#ifndef SFMTAVXEQ_HPP
#define SFMTAVXEQ_HPP
/**
 * @file SFMTAVXeq.hpp
 */
#include "devavxprng.h"
#include "AlgorithmSIMDEquidistribution.hpp"
//#include <MTToolBox/AlgorithmEquidistribution.hpp>
//#include <MTToolBox/AlgorithmReducibleRecursionSearch.hpp>
//#include <MTToolBox/period.hpp>
//#include <NTL/GF2X.h>
#include "Annihilate.hpp"
#include "EQOptions.hpp"

namespace MTToolBox {

    template<typename U, typename G, typename P, int bitWidth>
    int sfmtavx_equidistribution(EQOptions<P> opt)
    {
        using namespace std;
        G sf(opt.params);
        cout << sf.getParamString() << endl;
        U wseed;
        wseed.u64[0] = opt.seed;
        sf.seed(wseed);
        Annihilate<G, U, uint32_t> annihilate;
        if (!annihilate.anni(sf)) {
            return -1;
        }
        //cout << "annihilate end" << endl;
        int delta32 = 0;
        int delta64 = 0;
        bool lsb = false;
        const char * lsb_str = "";
#if 0
        if (opt.reverse) {
            sf.set_reverse_bit();
            lsb = true;
            cout << "Equidistribution from LSB" << endl;
            lsb_str = " from LSB";
        }
#endif
        SIMDInfo info;
        info.bitSize = bitWidth;
        info.fastMode = false;
        sf.reset_reverse_bit();
        info.bitMode = 64;
        info.elementNo = bitWidth / 64;
        if (opt.verbose) {
            cout << "64bit dimension of equidistribution at v-bit accuracy k(v)"
                 << lsb_str << endl;
        }
        for (int v = 1; v <= 64; v++) {
            int veq = calc_SIMD_equidist<U, G>(v, sf, info,
                                               opt.params.mexp, lsb);
            int d = opt.params.mexp / v - veq;
            if (opt.verbose) {
                cout << "k(" << dec << v << ") = " << dec << veq;
                cout << "\td(" << dec << v << ") = " << dec << d << endl;
            }
            delta64 += d;
        }
        info.bitMode = 32;
        info.elementNo = bitWidth / 32;
        if (opt.verbose) {
            cout << "32bit dimension of equidistribution at v-bit accuracy k(v)"
                 << lsb_str << endl;
        }
        for (int v = 1; v <= 32; v++) {
            int veq = calc_SIMD_equidist<U, G>
                (v, sf, info, opt.params.mexp, lsb);
            int d = opt.params.mexp / v - veq;
            if (opt.verbose) {
                cout << "k(" << dec << v << ") = " << dec << veq;
                cout << "\td(" << dec << v << ") = " << dec << d << endl;
            }
            delta32 += d;
        }
        cout << sf.getParamString();
        cout << dec << delta32 << "," << delta64 << endl;
        return 0;
    }
}
#endif // SFMTAVXEQ_HPP
