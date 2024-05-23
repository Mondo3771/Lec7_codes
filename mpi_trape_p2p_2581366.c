#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

double f(double x)
{
    double result = 4.0 / (1 + x * x);
    return result;
}

double Trapezoid(double x_1, double x_2, double h, int n)
{
    double result = (f(x_1) + f(x_2)) * h / 2;
    for (int i = 1; i <= n - 1; i++)
    {
        double x_i = x_1 + i * h;
        result += f(x_i);
    }
    result = h * result;
    return result;
}

int main(int argc, char *argv[])
{
    double a = atoi(argv[1]), b = atoi(argv[2]), number = atoi(argv[3]);
    double h = (b - a) / number;

    int rank, num_procs;
    double start = MPI_Wtime(), end;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    // n is divided between all processes
    double local_n = number / num_procs;
    double local_a = a + rank * local_n * h;
    double local_b = local_a + (local_n * h);
    double total;

    double local_result = Trapezoid(local_a, local_b, h, local_n);

    if (rank != 0)
    {
        MPI_Send(&local_result, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    else
    {
        total = local_result;
        for (int i = 1; i < num_procs; i++)
        {
            MPI_Recv(&local_result, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            total += local_result;
        }
        end = MPI_Wtime();

        printf("With n = %.0f trapezoids, our estimate of the integral from %.6f to %.6f = %.15f in %.6fs\n", number, a, b, total, end - start);
    }
    MPI_Finalize();
    return 0;
}