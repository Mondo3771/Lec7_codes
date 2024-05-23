#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

double f(double x)
{
    double result = 4 / (1 + x * x);
    return result;
}

double Trap(double x_1, double x_2, double h, int n)
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
    double start, end;
    double cpu_time_used;

    start = MPI_Wtime();
    MPI_Init(&argc, &argv);
    if (argc != 4)
    {
        printf("Usage: %s <a> <b> <n>\n", argv[0]);
        return 1;
    }
    int num_procs, myrank;
    double a = atof(argv[1]);
    double b = atof(argv[2]);
    double n = atof(argv[3]);
    double h;
    h = (b - a) / n;

    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    double local_approx = 0;
    double local_n = n / num_procs;
    double local_a = a + myrank * local_n * h;
    double local_b = local_a + local_n * h;
    local_approx = Trap(local_a, local_b, h, local_n);

    if (myrank != 0)
    {
        MPI_Send(&local_approx, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    else
    {
        double total_approx = local_approx;
        for (int i = 1; i < num_procs; i++)
        {
            MPI_Recv(&local_approx, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            total_approx += local_approx;
        }
        end = MPI_Wtime();
        cpu_time_used = ((double)(end - start));
        printf("The approximated value of the integral is: %f \n", total_approx);
        printf("With n = %.0f trapezoids, our estimate of the integral\n", n);
        printf("from %.6f to %.6f = %.15e in %.6fs\n", a, b, total_approx, cpu_time_used);
    }
    MPI_Finalize();

    return 0;
}