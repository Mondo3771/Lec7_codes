#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

int main(int argc, char *argv[])
{
    int num_procs, myrank;
    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    char g[12] = "Hello World!";
    MPI_Bcast(g, 12, MPI_CHAR, 0, MPI_COMM_WORLD);
    printf("Process %d received token %s from process %d.\n", myrank, g, 0);

    MPI_Finalize();
    return 0;
}