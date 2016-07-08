/**
 * @file dSFMTAVX2dc.cpp
 *
 * @brief The main function of parameter generator of dSFMTAVX2.
 *
 */
#include "devavxprng.h"
#include <string>
#include <sstream>
#include <fstream>
#include <MTToolBox/AlgorithmReducibleRecursionSearch.hpp>
#include <MTToolBox/AlgorithmCalculateParity.hpp>
#include <MTToolBox/AlgorithmEquidistribution.hpp>
#include <MTToolBox/MersenneTwister64.hpp>
#include <NTL/GF2X.h>
#include "dSFMTAVX2search.hpp"
#include "AlgorithmDSFMTEquidistribution.hpp"
#include "Annihilate.hpp"
#include "AlgorithmCalcFixPoint.hpp"
#include "dSFMTAVX2dc.h"

using namespace std;
using namespace MTToolBox;
using namespace NTL;

/**
 * search parameters using all_in_one function in the file search_all.hpp
 * @param opt command line options
 * @param count number of parameters user requested
 * @return 0 if this ends normally
 */
int search(options& opt, int count) {
    MersenneTwister64 mt(opt.seed);
    dSFMTAVX2 g(opt.mexp);

    g.setFixedSL1(19);
    g.setFixed(true);

    cout << "seed = " << dec << opt.seed << endl;
    if (opt.verbose) {
        time_t t = time(NULL);
        cout << "search start at " << ctime(&t);
    }
    if (opt.fixed) {
        g.setFixed(true);
        g.setFixedSL1(opt.fixedSL1);
    }
    AlgorithmReducibleRecursionSearch<w256_t, uint64_t> ars(g, mt);
    int i = 0;
    AlgorithmCalculateParity<w256_t, dSFMTAVX2> cp;
    Annihilate<dSFMTAVX2, w256_t, uint64_t> annihilate;
    cout << "# " << g.getHeaderString() << ", delta52"
         << endl;
    while (i < count) {
        if (ars.start(opt.mexp * 1000)) {
            GF2X irreducible = ars.getIrreducibleFactor();
            GF2X characteristic = ars.getCharacteristicPolynomial();
            //cout << "deg irreducible = " << dec << deg(irreducible) << endl;
            //cout << "deg characteristic = " << dec << deg(characteristic)
            //     << endl;
            if (deg(irreducible) != opt.mexp) {
                cout << "error mexp = " << dec << opt.mexp << " deg = "
                     << dec << deg(irreducible) << endl;
                return -1;
            }
            annihilate.getLCMPoly(characteristic, g);
            GF2X quotient = characteristic / irreducible;
            //cout << "deg quotient = " << dec << deg(quotient) << endl;
            w256_t fixpoint
                = calc_fixpoint<dSFMTAVX2, w256_t>(g, irreducible, quotient);
            g.setFixPoint(fixpoint);
            cp.searchParity(g, irreducible);
            w256_t seed = {{1, 0, 0, 0, 0, 0, 0, 0}};
            g.seed(seed);
            if (!annihilate.anni(g)) {
                return -1;
            }
            int veq52[52];
            DSFMTInfo info;
            info.bitSize = 256;
            info.elementNo = 4;
            int delta52
                = calc_dSFMT_equidistribution<w256_t, dSFMTAVX2>
                (g, veq52, 52, info, opt.mexp);
            cout << g.getParamString();
            cout << dec << delta52;
            cout << "," << dec << veq52[51] << endl;
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
    return i;
}

static void output_help(string& pgm);

/**
 * command line option parser
 * @param opt a structure to keep the result of parsing
 * @param argc number of command line arguments
 * @param argv command line arguments
 * @param start default start value
 * @return command line options have error, or not
 */
bool parse_opt(options& opt, int argc, char **argv) {
    opt.verbose = false;
    opt.mexp = 0;
    opt.count = 1;
    opt.seed = (uint64_t)clock();
    opt.filename = "";
    opt.fixed = false;
    opt.fixedSL1 = 19;
    int c;
    bool error = false;
    string pgm = argv[0];
    static struct option longopts[] = {
        {"verbose", no_argument, NULL, 'v'},
        {"file", required_argument, NULL, 'f'},
        {"count", required_argument, NULL, 'c'},
        {"seed", required_argument, NULL, 's'},
        {"fixed", optional_argument, NULL, 'x'},
        {NULL, 0, NULL, 0}};
    errno = 0;
    for (;;) {
        c = getopt_long(argc, argv, "vs:f:c:x::", longopts, NULL);
        if (error) {
            break;
        }
        if (c == -1) {
            break;
        }
        switch (c) {
        case 's':
            opt.seed = strtoull(optarg, NULL, 0);
            if (errno) {
                error = true;
                cerr << "seed must be a number" << endl;
            }
            break;
        case 'x':
            opt.fixed = true;
            //if ((optarg != NULL) && (strlen(optarg) > 0)) {
            if (optarg != NULL) {
                opt.fixedSL1 = strtoull(optarg, NULL, 0);
                if (errno) {
                    error = true;
                    cerr << "fixed sl1 must be a number" << endl;
                }
            } else {
                opt.fixedSL1 = 19;
            }
            break;
        case 'v':
            opt.verbose = true;
            break;
        case 'f':
            opt.filename = optarg;
            break;
        case 'c':
            opt.count = strtoll(optarg, NULL, 10);
            if (errno) {
                error = true;
                cerr << "count must be a number" << endl;
            }
            break;
        case '?':
        default:
            error = true;
            break;
        }
    }
    argc -= optind;
    argv += optind;
    if (argc < 1) {
        error = true;
    } else {
        long mexp = strtol(argv[0], NULL, 10);
        static const int allowed_mexp[] = {607, 1279, 2281, 3217, 4253, 4423,
                                           9689, 9941, 11213, 19937,
                                           21701, 44497, -1};
        if (! errno) {
            bool found = false;
            for (int i = 0; allowed_mexp[i] > 0; i++) {
                if (mexp == allowed_mexp[i]) {
                    found = true;
                    break;
                }
            }
            if (! found) {
                error = true;
            }
        }
        if (errno || error){
            error = true;
            cerr << "mexp must be one of ";
            for (int i = 0; allowed_mexp[i] > 0; i++) {
                cerr << dec << allowed_mexp[i] << " ";
            }
            cerr << endl;
        }
        opt.mexp = mexp;
    }
    if (!opt.filename.empty()) {
        ofstream ofs(opt.filename.c_str());
        if (ofs) {
            ofs.close();
        } else {
            error = true;
            cerr << "can't open file:" << opt.filename << endl;
        }
    }
    if (error) {
        output_help(pgm);
        return false;
    }
    return true;
}

/**
 * showing help message
 * @param pgm program name
 */
static void output_help(string& pgm) {
    cerr << "usage:" << endl;
    cerr << pgm
         << " [-s seed] [-v] [-c count]"
         << " [-f outputfile]"
         << " mexp"
         << endl;
    static string help_string1 = "\n"
"--verbose, -v        Verbose mode. Output parameters, calculation time, etc.\n"
"--file, -f filename  Parameters are outputted to this file. without this\n"
"                     option, parameters are outputted to standard output.\n"
"--count, -c count    Output count. The number of parameters to be outputted.\n"
"--seed, -s seed      seed of randomness.\n"
"--fixed, -x fixedSL  fix the parameter sl1 to given value.\n"
"mexp                 mersenne exponent, 1279 or more.\n"
        ;
    cerr << help_string1 << endl;
}

#if !defined(NO_MAIN)
int main(int argc, char** argv) {
    options opt;
    bool parse = parse_opt(opt, argc, argv);
    if (!parse) {
        return -1;
    }
    return search(opt, opt.count);
}
#endif
