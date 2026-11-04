#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "common.h"

static void fail(const char *msg) {
    fprintf(stderr, "ERROR: %s\n", msg);
    MPI_Abort(MPI_COMM_WORLD, 1);
}

/* Same maths as compute_accel in common.c, but x, y, z are separate arrays */
static void accel3(int lo, int hi, int N, const double *x, const double *y,
                   const double *z, const double *m,
                   double *ax, double *ay, double *az) {
    const double eps2 = EPS * EPS;
    for (int i = lo; i < hi; i++) {
        double xi = x[i], yi = y[i], zi = z[i];
        double fx = 0, fy = 0, fz = 0;
        for (int j = 0; j < N; j++) {
            double dx = x[j] - xi;
            double dy = y[j] - yi;
            double dz = z[j] - zi;
            double r2 = dx*dx + dy*dy + dz*dz + eps2;
            double inv = m[j] / (r2 * sqrt(r2));
            fx += dx * inv; fy += dy * inv; fz += dz * inv;
        }
        ax[i] = G * fx; ay[i] = G * fy; az[i] = G * fz;
    }
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank, P;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &P);

    if (argc < 3) fail("Usage: mpirun -np P ./bin/mpi_nbody_3calls N STEPS [outfile]");
    int N = atoi(argv[1]), T = atoi(argv[2]);
    const char *out = (argc > 3) ? argv[3] : NULL;
    if (N <= 0 || T <= 0) fail("N and STEPS must be > 0");
    if (P > N) fail("More processes than bodies");

    /* counts/displs are in BODIES here (each array holds 1 double per body) */
    int *counts = malloc(P * sizeof(int));
    int *displs = malloc(P * sizeof(int));
    if (!counts || !displs) fail("malloc failed");
    int base = N / P, rem = N % P, start = 0;
    for (int r = 0; r < P; r++) {
        counts[r] = base + (r < rem ? 1 : 0);
        displs[r] = start;
        start += counts[r];
    }
    int lo = displs[rank];
    int hi = lo + counts[rank];

    double *x  = malloc(N*sizeof(double)), *y  = malloc(N*sizeof(double)), *z  = malloc(N*sizeof(double));
    double *vx = malloc(N*sizeof(double)), *vy = malloc(N*sizeof(double)), *vz = malloc(N*sizeof(double));
    double *ax = malloc(N*sizeof(double)), *ay = malloc(N*sizeof(double)), *az = malloc(N*sizeof(double));
    double *m  = malloc(N*sizeof(double));
    double *pos = malloc(3*N*sizeof(double)), *vel = malloc(3*N*sizeof(double));
    if (!x||!y||!z||!vx||!vy||!vz||!ax||!ay||!az||!m||!pos||!vel) fail("malloc failed");

    /* same initial data as every other version, split into separate arrays */
    init_bodies(N, pos, vel, m);
    for (int i = 0; i < N; i++) {
        x[i]  = pos[3*i];  y[i]  = pos[3*i+1];  z[i]  = pos[3*i+2];
        vx[i] = vel[3*i];  vy[i] = vel[3*i+1];  vz[i] = vel[3*i+2];
    }
    double e0 = 0;
    if (rank == 0) e0 = total_energy(N, pos, vel, m);

    double comp = 0, comm = 0;
    MPI_Barrier(MPI_COMM_WORLD);
    double t_start = MPI_Wtime();

    double a = MPI_Wtime();
    accel3(lo, hi, N, x, y, z, m, ax, ay, az);
    comp += MPI_Wtime() - a;

    for (int s = 0; s < T; s++) {
        double t0 = MPI_Wtime();
        for (int i = lo; i < hi; i++) {
            vx[i] += 0.5*DT*ax[i]; vy[i] += 0.5*DT*ay[i]; vz[i] += 0.5*DT*az[i];
            x[i] += DT*vx[i];      y[i] += DT*vy[i];      z[i] += DT*vz[i];
        }
        double t1 = MPI_Wtime();

        /* THREE collective calls per step (the thing we are measuring) */
        MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DOUBLE, x, counts, displs, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DOUBLE, y, counts, displs, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DOUBLE, z, counts, displs, MPI_DOUBLE, MPI_COMM_WORLD);
        double t2 = MPI_Wtime();

        accel3(lo, hi, N, x, y, z, m, ax, ay, az);
        for (int i = lo; i < hi; i++) {
            vx[i] += 0.5*DT*ax[i]; vy[i] += 0.5*DT*ay[i]; vz[i] += 0.5*DT*az[i];
        }
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

    /* untimed: gather velocities so rank 0 can compute final energy */
    MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DOUBLE, vx, counts, displs, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DOUBLE, vy, counts, displs, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DOUBLE, vz, counts, displs, MPI_DOUBLE, MPI_COMM_WORLD);

    if (rank == 0) {
        for (int i = 0; i < N; i++) {
            pos[3*i] = x[i];  pos[3*i+1] = y[i];  pos[3*i+2] = z[i];
            vel[3*i] = vx[i]; vel[3*i+1] = vy[i]; vel[3*i+2] = vz[i];
        }
        double e1 = total_energy(N, pos, vel, m);
        double imb = comp_max / (comp_sum / P);
        printf("%d,%d,%d,%.6f,%.6f,%.6f,%.4f,%.10e,%.10e\n",
               N, P, T, tot_max, comp_max, comm_max, imb, e0, e1);
        if (out) write_positions(out, N, pos);
    }

    free(x); free(y); free(z); free(vx); free(vy); free(vz);
    free(ax); free(ay); free(az); free(m); free(pos); free(vel);
    free(counts); free(displs);
    MPI_Finalize();
    return 0;
}
