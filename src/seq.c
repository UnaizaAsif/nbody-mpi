#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "common.h"

static double now(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s N STEPS [outfile]\n", argv[0]);
        return 1;
    }
    int N = atoi(argv[1]), T = atoi(argv[2]);
    const char *out = (argc > 3) ? argv[3] : NULL;
    if (N <= 0 || T <= 0) { fprintf(stderr, "N and STEPS must be > 0\n"); return 1; }

    double *pos = malloc(3*N*sizeof(double));
    double *vel = malloc(3*N*sizeof(double));
    double *acc = malloc(3*N*sizeof(double));
    double *m   = malloc(N*sizeof(double));
    if (!pos || !vel || !acc || !m) { fprintf(stderr, "malloc failed\n"); return 1; }

    init_bodies(N, pos, vel, m);
    double e0 = total_energy(N, pos, vel, m);

    double t0 = now();
    compute_accel(0, N, N, pos, m, acc);
    for (int s = 0; s < T; s++) {
        for (int k = 0; k < 3*N; k++) {           /* half kick + drift */
            vel[k] += 0.5 * DT * acc[k];
            pos[k] += DT * vel[k];
        }
        compute_accel(0, N, N, pos, m, acc);      /* new forces */
        for (int k = 0; k < 3*N; k++)             /* half kick */
            vel[k] += 0.5 * DT * acc[k];
    }
    double total = now() - t0;

    double e1 = total_energy(N, pos, vel, m);
    printf("%d,%d,%.6f,%.10e,%.10e\n", N, T, total, e0, e1);
    if (out) write_positions(out, N, pos);

    free(pos); free(vel); free(acc); free(m);
    return 0;
}

