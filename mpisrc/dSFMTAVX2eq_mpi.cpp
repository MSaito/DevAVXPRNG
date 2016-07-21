#include "devavxprng.h"
#include <mpi.h>
#include <fcntl.h>
#include "dSFMTAVX2search.hpp"
#include "dSFMTAVXeqmpi.hpp"
#include "EQOptions.hpp"

using namespace MTToolBox;
using namespace std;

int main(int argc, char * argv[])
{
    int rank;
    int num_process;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_process);

    EQOptions<dSFMTAVX2_param> opt;
    if (!opt.parse(argc, argv)) {
        MPI_Finalize();
        return -1;
    }

    // file manipulation
    char * pgm = argv[0];
    char fname[500];
    sprintf(fname, "%s-%d-%04d.txt", pgm, opt.params.mexp, (uint32_t)opt.seed);
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

    int r = dsfmtavxmpi_equidist<w256_t, dSFMTAVX2,
                                dSFMTAVX2_param, 256>(opt, rank, num_process);
    MPI_Finalize();
    return r;
}
