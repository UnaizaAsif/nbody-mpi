#ifndef BH_H
#define BH_H

/* Barnes-Hut acceleration for bodies lo..hi-1, using ALL N bodies.
   Builds the octree from pos each call. theta = opening angle (0.5 typical). */
void bh_accel(int lo, int hi, int N, const double *pos, const double *m,
              double *acc, double theta);
void bh_free(void);

#endif
