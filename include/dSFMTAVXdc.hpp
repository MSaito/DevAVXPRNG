#pragma once
#ifndef DSFMTAVXDC_HPP
#define DSFMTAVXDC_HPP
/**
 * @file dSFMTAVXdc.hpp
 */

#include "devavxprng.h"
#include <MTToolBox/AlgorithmReducibleRecursionSearch.hpp>
#include <MTToolBox/AlgorithmCalculateParity.hpp>
#include <MTToolBox/AlgorithmEquidistribution.hpp>
#include <MTToolBox/MersenneTwister64.hpp>
#include <NTL/GF2X.h>
#include "AlgorithmDSFMTEquidistribution.hpp"
#include "Annihilate.hpp"
#include "AlgorithmCalcFixPoint.hpp"
#include "DCOptions.hpp"

namespace MTToolBox {
    /**
     * search parameters using all_in_one function in the file search_all.hpp
     * @param opt command line options
     * @param count number of parameters user requested
     * @return 0 if this ends normally
     */
    template<typename U, typename G, int bitWidth>
    int dsfmtavx_search(DCOptions& opt, int count) {
        using namespace std;
        using namespace NTL;

        MersenneTwister64 mt(opt.seed);
        G g(opt.mexp);

        cout << "#seed = " << dec << opt.seed << endl;
        if (opt.verbose) {
            time_t t = time(NULL);
            cout << "search start at " << ctime(&t);
        }
        if (opt.fixedL) {
            g.setFixedSL1(opt.fixedSL1);
        }
        if (opt.fixedP) {
            g.setFixedPerm(opt.fixedPerm);
        }
        AlgorithmReducibleRecursionSearch<U> ars(g, mt);
        int i = 0;
        AlgorithmCalculateParity<U, G> cp;
        Annihilate<G, U> annihilate;
        cout << "# " << g.getHeaderString() << ", delta52"
             << endl;
        while (i < count) {
            if (ars.start(opt.mexp * 1000)) {
                GF2X irreducible = ars.getIrreducibleFactor();
                GF2X characteristic = ars.getCharacteristicPolynomial();
                if (deg(irreducible) != opt.mexp) {
                    cout << "error mexp = " << dec << opt.mexp << " deg = "
                         << dec << deg(irreducible) << endl;
                    return -1;
                }
                annihilate.getLCMPoly(characteristic, g);
                GF2X quotient = characteristic / irreducible;
                U fixpoint
                    = calc_fixpoint<G, U>(g, irreducible, quotient);
                g.setFixPoint(fixpoint);
                cp.searchParity(g, irreducible);
                U seed;
                setBitOfPos(&seed, 0, 1);
                g.seed(seed);
                if (!annihilate.anni(g)) {
                    return -1;
                }
                int veq52[52];
                DSFMTInfo info;
                info.bitSize = bitWidth;
                info.elementNo = bitWidth / 64;
                int delta52
                    = calc_dSFMT_equidistribution<U, G>
                    (g, veq52, 52, info, opt.mexp);
                cout << g.getParamString();
                cout << dec << delta52;
                cout << "," << veq52[51] << endl;
                //cout << endl;
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

#endif // DSFMTAVX512FDC_HPP
