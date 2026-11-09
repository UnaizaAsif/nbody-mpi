/* mpi_live.c - MPI N-body that streams positions for real-time visualization.
 *
 * Same decomposition, Allgatherv and velocity-Verlet steps as mpi_nbody.c.
 * After each Allgatherv, rank 0 already holds every position, so every EVERY
 * steps it writes one frame to OUTFILE (normally a named pipe read by
 * scripts/live_view.py). No extra MPI communication is needed for the display.
 *
 * Usage: mpirun -np P ./bin/mpi_live N STEPS EVERY OUTFILE
 * Easier: ./scripts/live.sh N STEPS EVERY P
 *
 * Stream format (binary, native byte order):
 *   header: int N, int P, int bodies_per_rank[P]
 *   frame:  int step, double sim_time_so_far, double comm_time_so_far, double pos[3N]
 *   end:    int -1
 */
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
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
    signal(SIGPIPE, SIG_IGN);   /* closing the viewer must not kill rank 0 */

    if (argc < 5) fail("Usage: mpirun -np P ./bin/mpi_live N STEPS EVERY OUTFILE");
    int N = atoi(argv[1]), T = atoi(argv[2]), every = atoi(argv[3]);
    const char *out = argv[4];
    if (N <= 0 || T <= 0 || every <= 0) fail("N, STEPS and EVERY must be > 0");
    if (P > N) fail("More processes than bodies");

    /* block decomposition, counts/displs in doubles (3 per body) */
    int *counts = malloc(P * sizeof(int));
    int *displs = malloc(P * sizeof(int));
    int *nbody  = malloc(P * sizeof(int));
    if (!counts || !displs || !nbody) fail("malloc failed");
    int base = N / P, rem = N % P, start = 0;
    for (int r = 0; r < P; r++) {
        nbody[r]  = base + (r < rem ? 1 : 0);
        counts[r] = 3 * nbody[r];
        displs[r] = 3 * start;
        start += nbody[r];
    }
    int lo = displs[rank] / 3, hi = lo + nbody[rank];

    double *pos = malloc(3*N*sizeof(double));
    double *vel = malloc(3*N*sizeof(double));
    double *acc = malloc(3*N*sizeof(double));
    double *m   = malloc(N*sizeof(double));
    if (!pos || !vel || !acc || !m) fail("malloc failed");
    init_bodies(N, pos, vel, m);

    FILE *f = NULL;
    int viewer_ok = 1;
    double t_start = MPI_Wtime(), comm = 0;
    if (rank == 0) {
        f = fopen(out, "wb");                     /* blocks until the viewer opens the pipe */
        if (!f) fail("cannot open OUTFILE");
        fwrite(&N, sizeof(int), 1, f);
        fwrite(&P, sizeof(int), 1, f);
        fwrite(nbody, sizeof(int), P, f);
        int s0 = 0; double z = 0;
        fwrite(&s0, sizeof(int), 1, f); fwrite(&z, sizeof(double), 1, f);
        fwrite(&z, sizeof(double), 1, f); fwrite(pos, sizeof(double), 3*N, f);
        fflush(f);
        fprintf(stderr, "streaming N=%d with P=%d processes to %s\n", N, P, out);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    t_start = MPI_Wtime();

    compute_accel(lo, hi, N, pos, m, acc);
    for (int s = 1; s <= T; s++) {
        for (int k = 3*lo; k < 3*hi; k++) {       /* half kick + drift (own block) */
            vel[k] += 0.5 * DT * acc[k];
            pos[k] += DT * vel[k];
        }
        double t1 = MPI_Wtime();
        MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DOUBLE,
                       pos, counts, displs, MPI_DOUBLE, MPI_COMM_WORLD);
        comm += MPI_Wtime() - t1;
        compute_accel(lo, hi, N, pos, m, acc);
        for (int k = 3*lo; k < 3*hi; k++)
            vel[k] += 0.5 * DT * acc[k];

        if (rank == 0 && viewer_ok && s % every == 0) {   /* rank 0 already has all positions */
            double el = MPI_Wtime() - t_start;
            fwrite(&s, sizeof(int), 1, f);
            fwrite(&el, sizeof(double), 1, f);
            fwrite(&comm, sizeof(double), 1, f);
            fwrite(pos, sizeof(double), 3*N, f);
            if (fflush(f) != 0 || ferror(f)) {
                viewer_ok = 0;
                fprintf(stderr, "viewer closed at step %d; finishing the run without display\n", s);
            }
        }
    }

    if (rank == 0) {
        int end = -1;
        if (viewer_ok) fwrite(&end, sizeof(int), 1, f);
        fclose(f);
        fprintf(stderr, "done: %d steps in %.2f s\n", T, MPI_Wtime() - t_start);
    }
    free(pos); free(vel); free(acc); free(m); free(counts); free(displs); free(nbody);
    MPI_Finalize();
    return 0;
}
