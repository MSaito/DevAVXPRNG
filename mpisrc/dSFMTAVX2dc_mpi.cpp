/**
 * dSFMTAVX2dc-mpi.cpp
 */
#include "devavxprng.h"
#include <mpi.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include "dSFMTAVX2search.hpp"
#include "DCOptions.hpp"
#include "dSFMTAVXdc.hpp"

using namespace MTToolBox;
using namespace std;

int main(int argc, char *argv[]) {
    int rank;
    int num_process;
    // MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_process);
    DCOptions opt(607);
    opt.useSR1 = false;
    opt.fixedSL1 = 19;
    opt.fixedPerm = 1;
    bool parse = opt.parse(argc, argv);
    if (!parse) {
        MPI_Finalize();
        return -1;
    }
    char * pgm = argv[0];
    char fname[500];
    sprintf(fname, "%s-%d-%04d.txt", pgm, opt.mexp, rank);
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
    opt.seed = opt.seed + rank * 127;
    dsfmtavx_search<w256_t, dSFMTAVX2, 256>(opt, opt.count);
    close(fd);
    MPI_Abort(MPI_COMM_WORLD, 0);
    MPI_Finalize();
    return 0;
}
