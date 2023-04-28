#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"
#include "func.c"
#pragma comment(lib, "C:/Program Files (x86)/Microsoft SDKs/MPI/Lib/x64/msmpi.lib")
#define NT (int)(t_max / t_step)
#define NH (int)(x_max / x_step)

double u_cross(double ** u, int k, int m) {
    if (k == 0){ return 2 * t_step * func(k * t_step, m * x_step) - t_step / x_step * (u[1][m + 1] - u[1][m - 1]); }
    if (m == 0) { return t_step * func(k * t_step, m * x_step) - t_step / x_step * (u[1][m + 1] - u[1][m]) + u[1][m]; }
    else { if (m == NH - 1) { return t_step * func(k * t_step, m * x_step) - t_step / x_step * (u[1][m] - u[1][m-1]) + u[1][m]; } }
    if (m >= NH) {
        printf("\n\nINDICES PROBLEM\n\n");
        return 0;
    }
    return 2 * t_step * func(k * t_step, m * x_step) - t_step / x_step * (u[1][m + 1] - u[1][m - 1]) + u[0][m];
};

void laba(int rank, int nproc, int argc, char** argv) {
    printf("IT'S WORKING, proc %i\n", rank);
    int len  = 2 + (NH - 1) / (nproc - 1); // 2 - для координат, -1 - первый элемент знаем
    int lenn = 2 + (NH - 1) / (nproc - 1) + (NH - 1) % (nproc - 1);
    if (rank == 0) {
        double** u = malloc(NT * sizeof(double*));
        for (int i = 0; i < NT; i++) {
            u[i] = calloc(NH, sizeof(double));
        }
        double* kmu  = malloc(len  * sizeof(double));
        double* kmun = malloc(lenn * sizeof(double));
        for (int i = 0; i < NT; i++) {
            u[i][0] = ksi(i * t_step);
        }
        for (int i = 0; i < NH; i++) {
            u[0][i] = fi(i * x_step);
        }
        for (int i = 1; i < nproc; i++) {
            MPI_Send(&u[0], NH, MPI_DOUBLE, i, 5, MPI_COMM_WORLD);
            MPI_Send(&u[0], NH, MPI_DOUBLE, i, 5, MPI_COMM_WORLD);
            if (i < nproc - 1) {
                MPI_Recv(kmu, len, MPI_DOUBLE, i, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for (int j = 2; j < len; j++) {
                    int m = (int)kmu[1] + j - 2;
                    u[1][m] = kmu[j];
                }
            }
            else {
                MPI_Recv(kmun, lenn, MPI_DOUBLE, i, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for (int j = 2; j < lenn; j++) {
                    int m = (int)kmun[1] + j - 2;
                    u[1][m] = kmun[j];
                }
            }
        }
        for (int iter = 2; iter < NT; iter++) {
            for (int i = 1; i < nproc; i++) {
                MPI_Send(u[iter - 1], NH, MPI_DOUBLE, i, 5, MPI_COMM_WORLD);
                MPI_Send(u[iter], NH, MPI_DOUBLE, i, 5, MPI_COMM_WORLD);
                if (i < nproc - 1) {
                    MPI_Recv(kmu, len, MPI_DOUBLE, i, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    int k = (int)kmu[0];
                    if (k == 0) { printf("NOOOOO\n"); }
                    for (int j = 2; j < len; j++) {
                        int m = (int)kmu[1] + j - 2;
                        u[k][m] = kmu[j];
                    }
                }
                else {
                    MPI_Recv(kmun, lenn, MPI_DOUBLE, i, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    int k = (int)kmun[0];
                    if (k == 0) { printf("NOOOOO\n"); }
                    for (int j = 2; j < lenn; j++) {
                        int m = (int)kmun[1] + j - 2;
                        u[k][m] = kmun[j];
                    }
                }
            }
        }
        //printf("\nkmu:\n");
        //for (int i = 0; i < lenn; i++) {
        //    printf("%f ", kmu[i]);
        //}
        //printf("\n");
        // 
        //for (int i = 0; i < NT; i++) {
        //    for (int j = 0; j < NH; j++) {
        //        printf("%lf ", u[i][j]);
        //    }
        //    printf("\n");
        //}

        FILE * file;
        file = fopen("C:\\parproga\\Parproga\\output.csv", "w+");
        if (file) {
            printf("writing to file");
            fprintf(file, "x,t,u\n"); // записываем заголовок
            for (int k = 0; k < NT; k++) {
                double t = k *  t_step ;
                for (int m = 0; m < NH; m++) {
                    double x = m * x_step;
                    fprintf(file, "%.6f,%.6f,%.6f\n", x, t, u[k][m]);
                }
            }
            fclose(file);
        }
        else {
            printf("Failed to open the file\n");
        }
    }
    else {
        double** u = malloc(2 * sizeof(double*));
        for (int i = 0; i < 2; i++) {
            u[i] = malloc(NH * sizeof(double));
        }
        int l = len;
        if (rank == nproc - 1) {
            l = lenn;
        }

        double* kmu = calloc(l, sizeof(double));

        for (int iter = 1; iter < NT; iter++) {
            MPI_Recv(u[0], NH, MPI_DOUBLE, 0, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(u[1], NH, MPI_DOUBLE, 0, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            if (kmu != NULL) {

                kmu[0] = iter;
                kmu[1] = (double)((rank - 1.0) * (len - 2.0) + 1.0); // +1, тк первый элемент знаем
                if (kmu[1] >= NH) {
                    printf("AAAAAAAAA\n");
                }
                for (int i = 2; i < l; i++) { kmu[i] = u_cross(u, iter, (i - 2 + (int)kmu[1])); }
                //for (int i = 2; i < l; i++) { kmu[i] = u_cross(u, (i - 2 + (int)kmu[1]), iter); }
                MPI_Send(kmu, l, MPI_DOUBLE, 0, 5, MPI_COMM_WORLD);
            }
        }
        free(kmu);
    }
};

int main(int argc, char** argv) {
    int rank;
    int nproc;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    if (nproc == 1) { 
        MPI_Finalize();
        return 0; 
    }

    laba(rank, nproc, argc, argv);

    MPI_Finalize();
    return 0;
}