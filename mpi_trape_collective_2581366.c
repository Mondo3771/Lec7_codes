#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

double f(double x)
{
    return (double)(4.0 / (1 + x * x));
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

    int rank, num_procs;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    double start = MPI_Wtime(), end;
    double a = atof(argv[1]);
    double b = atof(argv[2]);
    double number = atof(argv[3]);
    double h = (b - a) / number;
    double local_n = number / num_procs;
    double local_a = a + rank * local_n * h;
    double local_b = local_a + (local_n * h);
    double total;

    double local_result = Trapezoid(local_a, local_b, h, local_n);
    MPI_Reduce(&local_result, &total, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        end = MPI_Wtime();
        printf("With n = %.0f trapezoids, our estimate of the integral from %.6f to %.6f = %.15f in %.6fs\n", number, a, b, total, end - start);
    }
    MPI_Finalize();

    return 0;
}