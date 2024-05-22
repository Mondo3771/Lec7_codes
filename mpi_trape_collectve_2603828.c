#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

int main(int argc, char *argv[])
{
    int num_procs, myrank, data[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    printf("FROM process %d out of %d , Hello World!\n", myrank, num_procs);
    if (myrank == 0)
    {
        MPI_Send(data, 10, MPI_INT, 1, 0, MPI_COMM_WORLD);
    }
    else
        f(myrank == 1)
        {
            MPI_Recv(data, 10, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    MPI_Finalize();
    return 0;
}