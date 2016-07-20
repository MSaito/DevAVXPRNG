#include "devavxprng.h"
#include <mpi.h>
#include "SFMTAVX512Fsearch.hpp"
#include "AlgorithmSIMDEquidistribution.hpp"
#include <MTToolBox/AlgorithmEquidistribution.hpp>
#include <MTToolBox/AlgorithmReducibleRecursionSearch.hpp>
#include <MTToolBox/period.hpp>
#include <NTL/GF2X.h>
#include <errno.h>
#include "Annihilate.hpp"

using namespace MTToolBox;
using namespace std;

class options {
public:
    bool verbose;
    bool reverse;
    uint64_t seed;
    SFMTAVX512F_param params;
};

static bool parse_opt(options& opt, int argc, char **argv);
static void output_help(string& pgm);

int main(int argc, char * argv[])
{
    int rank;
    int num_process;
    // MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_process);

    options opt;
    if (!parse_opt(opt, argc, argv)) {
        MPI_Finalize();
        return -1;
    }
    SFMTAVX512F sf(opt.params);
    w512_t wseed;
    wseed.u64[0] = opt.seed;
    sf.seed(wseed);
    Annihilate<SFMTAVX512F, w512_t, uint32_t> annihilate;
    if (!annihilate.anni(sf)) {
        MPI_Finalize();
        return -1;
    }

    // file manipulation
    char * pgm = argv[0];
    char fname[500];
    sprintf(fname, "%s-%d-%04d.txt", pgm, opt.mexp, seed);
    int fd;
    errno = 0;
    fd = open(fname, O_WRONLY | O_CREAT, 0660);
    if (fd < 0) {
        printf("%s: open file error\n", pgm);
        MPI_Finalize();
        return 1;
    }
    dup2(fd, 1);
    if (errno) {
        printf("%s: dup error.\n", argv[0]);
        close(fd);
        MPI_Finalize();
        return 1;
    }

    bool lsb = false;
    const char * lsb_str = "";
    if (opt.reverse) {
        sf.set_reverse_bit();
        lsb = true;
        cout << "Equidistribution from LSB" << endl;
        lsb_str = " from LSB";
    }
    SIMDInfo info;
    info.fastMode = false;
    info.bitSize = 512;

    info.bitMode = 64;
    info.elementNo = 8;
    if (rank == 0) {
        cout << sf.getParamString() << endl;
        cout << "64bit dimension of equidistribution at v-bit accuracy k(v)"
             << lsb_str << endl;
    }
    for (int v = 1; v <= 64; v++) {
        if (v % num_process == rank) {
            int veq = calc_SIMD_equidistribution<w512_t, SFMTAVX512F>
                (v, sf, info, opt.params.mexp, lsb);
            int d = opt.params.mexp / v - veq;
            cout << "k(" << dec << v << ") = " << dec << veq;
            cout << "\td(" << dec << v << ") = " << dec << d << endl;
        }
    }
    info.bitMode = 32;
    info.elementNo = 16;
    sf.reset_reverse_bit();
    if (rank == 0) {
        cout << "32bit dimension of equidistribution at v-bit accuracy k(v)"
             << lsb_str << endl;
    }
    for (int v = 1; v <= 32; v++) {
        if (v % num_process == rank) {
            int veq = calc_SIMD_equidistribution<w512_t, SFMTAVX512F>
                (v, sf, info, opt.params.mexp, lsb);
            int d = opt.params.mexp / v - veq;
            cout << "k(" << dec << v << ") = " << dec << veq;
            cout << "\td(" << dec << v << ") = " << dec << d << endl;
        }
    }
    MPI_Finalize();
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
    opt.reverse = false;
    int c;
    bool error = false;
    string pgm = argv[0];
    static struct option longopts[] = {
        {"verbose", no_argument, NULL, 'v'},
        {"seed", required_argument, NULL, 's'},
        {"reverse", no_argument, NULL, 'r'},
        {NULL, 0, NULL, 0}};
    errno = 0;
    for (;;) {
        c = getopt_long(argc, argv, "vrs:", longopts, NULL);
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
        case 'r':
            opt.reverse = true;
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
         << " [-v] [-r]"
         << " mexp,pos1,sl1,sr1,msk1,parity1"
         << endl;
    static string help_string1 = "\n"
"--verbose, -v        Verbose mode. Output detailed information.\n"
"--reverse, -r        Reverse mode. Calculate equidistribution from LSB.\n"
        ;
    cerr << help_string1 << endl;
}
