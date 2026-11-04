#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "common.h"

/* Same seed => every run and every MPI rank gets identical bodies */
void init_bodies(int N, double *pos, double *vel, double *m) {
    srand48(SEED);
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < 3; k++) pos[3*i+k] = 2.0 * drand48() - 1.0;
        for (int k = 0; k < 3; k++) vel[3*i+k] = 0.2 * drand48() - 0.1;
        m[i] = 1.0 / N;
    }
}

/* Acceleration for bodies lo..hi-1, against ALL N bodies.
   Results stored at global index (acc[3*i]), not i-lo. */
void compute_accel(int lo, int hi, int N, const double *pos,
                   const double *m, double *acc) {
    const double eps2 = EPS * EPS;
    for (int i = lo; i < hi; i++) {
        double xi = pos[3*i], yi = pos[3*i+1], zi = pos[3*i+2];
        double fx = 0, fy = 0, fz = 0;
        for (int j = 0; j < N; j++) {
            double dx = pos[3*j]   - xi;
            double dy = pos[3*j+1] - yi;
            double dz = pos[3*j+2] - zi;
            double r2 = dx*dx + dy*dy + dz*dz + eps2;
            double inv = m[j] / (r2 * sqrt(r2));
            fx += dx * inv; fy += dy * inv; fz += dz * inv;
        }
        acc[3*i] = G * fx; acc[3*i+1] = G * fy; acc[3*i+2] = G * fz;
    }
}

double total_energy(int N, const double *pos, const double *vel,
                    const double *m) {
    double ke = 0, pe = 0;
    const double eps2 = EPS * EPS;
    for (int i = 0; i < N; i++) {
        double v2 = vel[3*i]*vel[3*i] + vel[3*i+1]*vel[3*i+1] + vel[3*i+2]*vel[3*i+2];
        ke += 0.5 * m[i] * v2;
    }
    for (int i = 0; i < N; i++)
        for (int j = i + 1; j < N; j++) {
            double dx = pos[3*j]-pos[3*i], dy = pos[3*j+1]-pos[3*i+1], dz = pos[3*j+2]-pos[3*i+2];
            pe -= G * m[i] * m[j] / sqrt(dx*dx + dy*dy + dz*dz + eps2);
        }
    return ke + pe;
}

void write_positions(const char *fname, int N, const double *pos) {
    FILE *f = fopen(fname, "w");
    if (!f) { fprintf(stderr, "Cannot open %s\n", fname); exit(1); }
    for (int i = 0; i < N; i++)
        fprintf(f, "%.17g %.17g %.17g\n", pos[3*i], pos[3*i+1], pos[3*i+2]);
    fclose(f);
}

