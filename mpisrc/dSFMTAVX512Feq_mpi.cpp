#include "devavxprng.h"
#include <mpi.h>
#include "dSFMTAVX512Fsearch.hpp"
#include "AlgorithmDSFMTEquidistribution.hpp"
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
    dSFMTAVX512F_param params;
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
    dSFMTAVX512F sf(opt.params);
    w512_t wseed;
    wseed.u64[0] = opt.seed;
    sf.seed(wseed);
    Annihilate<dSFMTAVX512F, w512_t, uint64_t> annihilate;
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

    DSFMTInfo info;
    info.bitSize = 512; // IMPORTANT
    info.elementNo = 8;
    if (rank == 0) {
        cout << sf.getParamString() << endl;
        cout << "52bit dimension of equidistribution at v-bit accuracy k(v)"
             << endl;
    }
    for (int v = 1; v <= 52; v++) {
        if (v % num_process == rank) {
            int veq = calc_dSFMT_equidistribution<w512_t, dSFMTAVX512F>
                (v, sf, info, opt.params.mexp);
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
         << " [-v] [-r]"
         << " mexp,pos1,sl1,msk1,msk2,fix1,fix2,"
         << "parity1,parity2"
         << endl;
    static string help_string1 = "\n"
"--verbose, -v        Verbose mode. Output detailed information.\n"
        ;
    cerr << help_string1 << endl;
}
