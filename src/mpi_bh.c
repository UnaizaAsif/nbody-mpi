#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "common.h"
#include "bh.h"

static void fail(const char *msg) {
    fprintf(stderr, "ERROR: %s\n", msg);
    MPI_Abort(MPI_COMM_WORLD, 1);
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank, P;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &P);

    if (argc < 4) fail("Usage: mpirun -np P ./bin/mpi_bh N STEPS THETA [outfile]");
    int N = atoi(argv[1]), T = atoi(argv[2]);
    double theta = atof(argv[3]);
    const char *out = (argc > 4) ? argv[4] : NULL;
    if (N <= 0 || T <= 0 || theta <= 0) fail("N, STEPS and THETA must be > 0");
    if (P > N) fail("More processes than bodies");

    int *counts = malloc(P * sizeof(int));
    int *displs = malloc(P * sizeof(int));
    if (!counts || !displs) fail("malloc failed");
    int base = N / P, rem = N % P, start = 0;
    for (int r = 0; r < P; r++) {
        int nb = base + (r < rem ? 1 : 0);
        counts[r] = 3 * nb;
        displs[r] = 3 * start;
        start += nb;
    }
    int lo = displs[rank] / 3;
    int hi = lo + counts[rank] / 3;

    double *pos = malloc(3*N*sizeof(double));
    double *vel = malloc(3*N*sizeof(double));
    double *acc = malloc(3*N*sizeof(double));
    double *m   = malloc(N*sizeof(double));
    if (!pos || !vel || !acc || !m) fail("malloc failed");

    init_bodies(N, pos, vel, m);
    double e0 = 0;
    if (rank == 0) e0 = total_energy(N, pos, vel, m);

    double comp = 0, comm = 0;
    MPI_Barrier(MPI_COMM_WORLD);
    double t_start = MPI_Wtime();

    double a = MPI_Wtime();
    bh_accel(lo, hi, N, pos, m, acc, theta);
    comp += MPI_Wtime() - a;

    for (int s = 0; s < T; s++) {
        double t0 = MPI_Wtime();
        for (int k = 3*lo; k < 3*hi; k++) {
            vel[k] += 0.5 * DT * acc[k];
            pos[k] += DT * vel[k];
        }
        double t1 = MPI_Wtime();

        MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DOUBLE,
                       pos, counts, displs, MPI_DOUBLE, MPI_COMM_WORLD);
        double t2 = MPI_Wtime();

        bh_accel(lo, hi, N, pos, m, acc, theta);   /* tree build + walk */
        for (int k = 3*lo; k < 3*hi; k++)
            vel[k] += 0.5 * DT * acc[k];
        double t3 = MPI_Wtime();

        comp += (t1 - t0) + (t3 - t2);
        comm += (t2 - t1);
    }
    double total = MPI_Wtime() - t_start;

    double tot_max, comp_max, comp_sum, comm_max;
    MPI_Reduce(&total, &tot_max,  1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&comp,  &comp_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&comp,  &comp_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&comm,  &comm_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DOUBLE,
                   vel, counts, displs, MPI_DOUBLE, MPI_COMM_WORLD);

    if (rank == 0) {
        double e1 = total_energy(N, pos, vel, m);
        double imb = comp_max / (comp_sum / P);
        printf("%d,%d,%d,%.6f,%.6f,%.6f,%.4f,%.10e,%.10e\n",
               N, P, T, tot_max, comp_max, comm_max, imb, e0, e1);
        if (out) write_positions(out, N, pos);
    }

    bh_free();
    free(pos); free(vel); free(acc); free(m); free(counts); free(displs);
    MPI_Finalize();
    return 0;
}
