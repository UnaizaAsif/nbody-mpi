#include <stdio.h>
#include <stdlib.h>
#include "common.h"

static void write_frame(FILE *f, int step, int N, const double *pos) {
    for (int i = 0; i < N; i++)
        fprintf(f, "%d,%d,%.10g,%.10g,%.10g\n",
                step, i, pos[3*i], pos[3*i+1], pos[3*i+2]);
}

int main(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr, "Usage: %s N STEPS EVERY [outfile]\n", argv[0]);
        return 1;
    }
    int N = atoi(argv[1]), T = atoi(argv[2]), every = atoi(argv[3]);
    const char *out = (argc > 4) ? argv[4] : "results/trajectory.csv";
    if (N <= 0 || T <= 0 || every <= 0) {
        fprintf(stderr, "N, STEPS and EVERY must be > 0\n");
        return 1;
    }

    double *pos = malloc(3*N*sizeof(double));
    double *vel = malloc(3*N*sizeof(double));
    double *acc = malloc(3*N*sizeof(double));
    double *m   = malloc(N*sizeof(double));
    if (!pos || !vel || !acc || !m) { fprintf(stderr, "malloc failed\n"); return 1; }

    FILE *f = fopen(out, "w");
    if (!f) { fprintf(stderr, "Cannot open %s\n", out); return 1; }
    fprintf(f, "step,body,x,y,z\n");

    init_bodies(N, pos, vel, m);
    compute_accel(0, N, N, pos, m, acc);
    write_frame(f, 0, N, pos);

    for (int s = 1; s <= T; s++) {
        for (int k = 0; k < 3*N; k++) {          /* half kick + drift */
            vel[k] += 0.5 * DT * acc[k];
            pos[k] += DT * vel[k];
        }
        compute_accel(0, N, N, pos, m, acc);     /* new forces */
        for (int k = 0; k < 3*N; k++)            /* half kick */
            vel[k] += 0.5 * DT * acc[k];
        if (s % every == 0) write_frame(f, s, N, pos);
    }

    fclose(f);
    printf("Saved %s\n", out);
    free(pos); free(vel); free(acc); free(m);
    return 0;
}
