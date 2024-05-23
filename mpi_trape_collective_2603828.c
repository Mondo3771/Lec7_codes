#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <time.h>

double f(double x)
{
    double result = 4.0 / (1 + x * x);
    return result;
}

double Trap(double x_1, double x_2, double h, int n)
{
    double result = (f(x_1) + f(x_2)) * h / 2;
    for (int i = 1; i < n; i++)
    {
        double x_i = x_1 + i * h;
        result += f(x_i);
    }
    result = h * result;
    return result;
}

int main(int argc, char *argv[])
{
    clock_t start, end;
    start = clock();
    printf("Starting\n");
    MPI_Init(&argc, &argv);
    int num_procs, myrank;

    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    double cpu_time_used;

    if (argc != 4)
    {
        printf("Usage: %s <a> <b> <n>\n", argv[0]);
        return 1;
    }
    double a, b, n, h;
    if (myrank == 0)
    {
        a = atof(argv[1]);
        b = atof(argv[2]);
        n = atof(argv[3]);
        h = (b - a) / n;
    }
    // approx = (Trap(a, b, h));
    // for (int i = 1; i < n - 1; i++)
    // {
    //     double x_i = a + i * h;
    //     approx += f(x_i);
    // }
    // approx = h * approx;
    // printf("The approximated value of the integral is: %f", approx);

    MPI_Bcast(&n, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&b, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&a, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&h, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    double total_approx = 0;
    double local_approx = 0;
    int local_n = n / num_procs;
    double local_a = a + myrank * local_n * h;
    double local_b = local_a + local_n * h;
    local_approx = Trap(local_a, local_b, h, local_n);
    MPI_Reduce(&local_approx, &total_approx, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (myrank == 0)
    {
        printf("The approximated value of the integral is: %f \n", total_approx);
    }
    MPI_Finalize();
    end = clock();
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("With n = %.0f trapezoids, our estimate of the integral\n", n);
    printf("from %.6f to %.6f = %.15e in %.6fs\n", a, b, total_approx, cpu_time_used);
    return 0;
}