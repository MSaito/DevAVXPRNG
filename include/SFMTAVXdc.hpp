#pragma once
#ifndef SFMTAVXDC_HPP
#define SFMTAVXDC_HPP
/**
 * @file SFMTAVXdc.hpp
 */
#include "devavxprng.h"
#include <MTToolBox/AlgorithmReducibleRecursionSearch.hpp>
#include <MTToolBox/AlgorithmCalculateParity.hpp>
#include <MTToolBox/AlgorithmEquidistribution.hpp>
#include <MTToolBox/MersenneTwister.hpp>
#include <NTL/GF2X.h>
#include "AlgorithmSIMDEquidistribution.hpp"
#include "Annihilate.hpp"
#include "DCOptions.hpp"

namespace MTToolBox {
    /**
     * search parameters using all_in_one function in the file search_all.hpp
     * @param opt command line options
     * @param count number of parameters user requested
     * @return 0 if this ends normally
     */
    template<typename U, typename G, int bitWidth>
    int sfmtavx_search(DCOptions& opt, int count) {
        using namespace std;
        using namespace NTL;
        MersenneTwister mt(opt.seed);
        G g(opt.mexp);

        if (opt.fixedL) {
            g.setFixedSL1(opt.fixedSL1);
        }
        if (opt.fixedR) {
            g.setFixedSR1(opt.fixedSR1);
        }
        if (opt.fixedP) {
            g.setFixedPerm(opt.fixedPerm);
        }
        cout << "seed = " << dec << opt.seed << endl;
        if (opt.verbose) {
            time_t t = time(NULL);
            cout << "search start at " << ctime(&t);
        }
        AlgorithmReducibleRecursionSearch<U, uint32_t> ars(g, mt);
        int i = 0;
        AlgorithmCalculateParity<U, G> cp;
        Annihilate<G, U, uint32_t> annihilate;
        cout << "# " << g.getHeaderString() << ", delta32, delta64"
             << endl;
        while (i < count) {
            if (ars.start(opt.mexp * 100)) {
                GF2X irreducible = ars.getIrreducibleFactor();
                GF2X characteristic = ars.getCharacteristicPolynomial();
                GF2X quotient = characteristic / irreducible;
                if (deg(irreducible) != opt.mexp) {
                    cout << "error not erreducible" << endl;
                    return -1;
                }
#if 0
                cout << "before parity" << endl;
#endif
                cp.searchParity(g, irreducible);
                U seed;
                setBitOfPos(&seed, 0, 1);
                g.seed(seed);
#if 0
                cout << "before annihilate" << endl;
#endif
                if (!annihilate.anni(g)) {
                    cout << "error can't annihilate" << endl;
                    return -1;
                }
                int delta32 = 0;
                int delta64 = 0;
                int veq32[32];
                SIMDInfo info;
                info.bitMode = 32;
                info.elementNo = bitWidth / 32;
                info.bitSize = bitWidth;
                info.fastMode = true;
#if 0
                cout << "bitMode:" << dec << info.bitMode << endl;
                cout << "elementNo:" << dec << info.elementNo << endl;
                cout << "bitSize:" << dec << info.bitSize << endl;
                cout << "fastmode:" << info.fastMode << endl;
#endif
                delta32
                    = calc_SIMD_equidistribution<U, G>
                    (g, veq32, 32, info, opt.mexp);
                int veq64[64];
                info.bitMode = 64;
                info.elementNo = bitWidth / 64;
#if 0
                cout << "bitMode:" << dec << info.bitMode << endl;
                cout << "elementNo:" << dec << info.elementNo << endl;
                cout << "bitSize:" << dec << info.bitSize << endl;
                cout << "fastmode:" << info.fastMode << endl;
#endif
                delta64
                    = calc_SIMD_equidistribution<U, G>
                    (g, veq64, 64, info, opt.mexp);
                cout << g.getParamString();
                cout << dec << delta32 << "," << delta64;
                cout << "," << dec << veq64[63] << endl;
                i++;
            } else {
                cout << "search failed" << endl;
                break;
            }
        }
        if (opt.verbose) {
            time_t t = time(NULL);
            cout << "search end at " << ctime(&t) << endl;
        }
        return 0;
    }

}
#endif // SFMTAVXDC_HPP
