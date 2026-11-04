#ifndef COMMON_H
#define COMMON_H

#define G    1.0     /* gravitational constant */
#define DT   0.01    /* time step */
#define EPS  0.05    /* softening */
#define SEED 42      /* random seed */

/* pos, vel, acc are interleaved arrays: x0,y0,z0,x1,y1,z1,... (size 3N) */
void   init_bodies(int N, double *pos, double *vel, double *m);
void   compute_accel(int lo, int hi, int N, const double *pos,
                     const double *m, double *acc);
double total_energy(int N, const double *pos, const double *vel,
                    const double *m);
void   write_positions(const char *fname, int N, const double *pos);

#endif
