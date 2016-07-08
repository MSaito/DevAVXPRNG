#include "devavxprng.h"
#include "dSFMTAVX2search.hpp"
#include "AlgorithmDSFMTEquidistribution.hpp"
#include <MTToolBox/AlgorithmEquidistribution.hpp>
#include <MTToolBox/AlgorithmReducibleRecursionSearch.hpp>
#include <MTToolBox/period.hpp>
#include <NTL/GF2X.h>
#include <errno.h>
#include <stdlib.h>
#include <getopt.h>
#include "Annihilate.hpp"

using namespace MTToolBox;
using namespace std;

class options {
public:
    bool verbose;
    uint64_t seed;
    dSFMTAVX2_param params;
};

static bool parse_opt(options& opt, int argc, char **argv);
static void output_help(string& pgm);

int main(int argc, char * argv[])
{
    options opt;
    if (!parse_opt(opt, argc, argv)) {
        return -1;
    }
    dSFMTAVX2 sf(opt.params);
    w256_t wseed;
    wseed.u64[0] = opt.seed;
    sf.seed(wseed);
    Annihilate<dSFMTAVX2, w256_t, uint64_t> annihilate;
    if (!annihilate.anni(sf)) {
        return -1;
    }
    int delta52 = 0;
    int veq52[52];
    DSFMTInfo info;
    info.bitSize = 256; // IMPORTANT
    info.elementNo = 4;
    delta52 = calc_dSFMT_equidistribution<w256_t, dSFMTAVX2>
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

/**
 * command line option parser
 * @param opt a structure to keep the result of parsing
 * @param argc number of command line arguments
 * @param argv command line arguments
 * @param start default start value
 * @return command line options have error, or not
 */
static bool parse_opt(options& opt, int argc, char **argv) {
    opt.verbose = false;
    opt.seed = (uint64_t)clock();
    int c;
    bool error = false;
    string pgm = argv[0];
    static struct option longopts[] = {
        {"verbose", no_argument, NULL, 'v'},
        {"seed", required_argument, NULL, 's'},
        {NULL, 0, NULL, 0}};
    errno = 0;
    for (;;) {
        c = getopt_long(argc, argv, "vs:", longopts, NULL);
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
        case 'v':
            opt.verbose = true;
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
        opt.params.readFromString(argv[0]);
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
static void output_help(string& pgm)
{
    cerr << "usage:" << endl;
    cerr << pgm
         << " [-v]"
         << " mexp,pos1,sl1,perm,msk1,fix1,"
         << "parity1"
         << endl;
    static string help_string1 = "\n"
"--verbose, -v        Verbose mode. Output detailed information.\n"
        ;
    cerr << help_string1 << endl;
}
