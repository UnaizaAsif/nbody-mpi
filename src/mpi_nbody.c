#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "common.h"

static void fail(const char *msg) {
    fprintf(stderr, "ERROR: %s\n", msg);
    MPI_Abort(MPI_COMM_WORLD, 1);
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank, P;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &P);

    /* ---- error handling ---- */
    if (argc < 3) fail("Usage: mpirun -np P ./bin/mpi_nbody N STEPS [outfile]");
    int N = atoi(argv[1]), T = atoi(argv[2]);
    const char *out = (argc > 3) ? argv[3] : NULL;
    if (N <= 0 || T <= 0) fail("N and STEPS must be > 0");
    if (P > N) fail("More processes than bodies");

    /* ---- data decomposition (handles N % P != 0) ----
       counts/displs are in DOUBLES (3 per body) because pos is interleaved */
    int *counts = malloc(P * sizeof(int));
    int *displs = malloc(P * sizeof(int));
    if (!counts || !displs) fail("malloc failed");
    int base = N / P, rem = N % P, start = 0;
    for (int r = 0; r < P; r++) {
        int nb = base + (r < rem ? 1 : 0);   /* first 'rem' ranks get +1 body */
        counts[r] = 3 * nb;
        displs[r] = 3 * start;
        start += nb;
    }
    int lo = displs[rank] / 3;
    int hi = lo + counts[rank] / 3;

    if (rank == 0)
        for (int r = 0; r < P; r++)
            fprintf(stderr, "rank %d owns bodies [%d, %d) = %d bodies\n",
                    r, displs[r]/3, displs[r]/3 + counts[r]/3, counts[r]/3);

    double *pos = malloc(3*N*sizeof(double));
    double *vel = malloc(3*N*sizeof(double));
    double *acc = malloc(3*N*sizeof(double));
    double *m   = malloc(N*sizeof(double));
    if (!pos || !vel || !acc || !m) fail("malloc failed");

    init_bodies(N, pos, vel, m);   /* same seed on every rank */
    double e0 = 0;
    if (rank == 0) e0 = total_energy(N, pos, vel, m);

    /* ---- timed simulation ---- */
    double comp = 0, comm = 0;
    MPI_Barrier(MPI_COMM_WORLD);
    double t_start = MPI_Wtime();

    double a = MPI_Wtime();
    compute_accel(lo, hi, N, pos, m, acc);
    comp += MPI_Wtime() - a;

    for (int s = 0; s < T; s++) {
        double t0 = MPI_Wtime();
        for (int k = 3*lo; k < 3*hi; k++) {       /* half kick + drift (own block) */
            vel[k] += 0.5 * DT * acc[k];
            pos[k] += DT * vel[k];
        }
        double t1 = MPI_Wtime();

        /* share new positions: everyone gets everyone's block */
        MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DOUBLE,
                       pos, counts, displs, MPI_DOUBLE, MPI_COMM_WORLD);
        double t2 = MPI_Wtime();

        compute_accel(lo, hi, N, pos, m, acc);    /* expensive part */
        for (int k = 3*lo; k < 3*hi; k++)
            vel[k] += 0.5 * DT * acc[k];
        double t3 = MPI_Wtime();

        comp += (t1 - t0) + (t3 - t2);
        comm += (t2 - t1);
    }
    double total = MPI_Wtime() - t_start;

    /* ---- statistics across ranks ---- */
    double tot_max, comp_max, comp_sum, comm_max;
    MPI_Reduce(&total, &tot_max,  1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&comp,  &comp_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&comp,  &comp_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&comm,  &comm_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    /* gather velocities (untimed) so rank 0 can compute final energy */
    MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DOUBLE,
                   vel, counts, displs, MPI_DOUBLE, MPI_COMM_WORLD);

    if (rank == 0) {
        double e1 = total_energy(N, pos, vel, m);
        double imb = comp_max / (comp_sum / P);
        printf("%d,%d,%d,%.6f,%.6f,%.6f,%.4f,%.10e,%.10e\n",
               N, P, T, tot_max, comp_max, comm_max, imb, e0, e1);
        if (out) write_positions(out, N, pos);
    }

    free(pos); free(vel); free(acc); free(m); free(counts); free(displs);
    MPI_Finalize();
    return 0;
}
