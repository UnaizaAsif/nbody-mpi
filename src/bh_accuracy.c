#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "common.h"
#include "bh.h"

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "Usage: %s N\n", argv[0]); return 1; }
    int N = atoi(argv[1]);
    if (N <= 0) { fprintf(stderr, "N must be > 0\n"); return 1; }

    double *pos = malloc(3*N*sizeof(double)), *vel = malloc(3*N*sizeof(double));
    double *ad  = malloc(3*N*sizeof(double)), *ab  = malloc(3*N*sizeof(double));
    double *m   = malloc(N*sizeof(double));
    if (!pos || !vel || !ad || !ab || !m) { fprintf(stderr, "malloc failed\n"); return 1; }

    init_bodies(N, pos, vel, m);
    compute_accel(0, N, N, pos, m, ad);               /* exact reference */

    double thetas[] = {0.2, 0.3, 0.5, 0.7};
    printf("N,theta,mean_rel_error,max_rel_error\n");
    for (int t = 0; t < 4; t++) {
        bh_accel(0, N, N, pos, m, ab, thetas[t]);
        double sum = 0, mx = 0;
        for (int i = 0; i < N; i++) {
            double dx = ab[3*i]-ad[3*i], dy = ab[3*i+1]-ad[3*i+1], dz = ab[3*i+2]-ad[3*i+2];
            double e = sqrt(dx*dx+dy*dy+dz*dz);
            double a = sqrt(ad[3*i]*ad[3*i]+ad[3*i+1]*ad[3*i+1]+ad[3*i+2]*ad[3*i+2]);
            double rel = e / (a + 1e-300);
            sum += rel; if (rel > mx) mx = rel;
        }
        printf("%d,%.1f,%.6e,%.6e\n", N, thetas[t], sum / N, mx);
    }
    bh_free();
    return 0;
}
