#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "common.h"
#include "bh.h"

#define MAXDEPTH 40

typedef struct {
    double cx, cy, cz, half;      /* cell centre and half-width            */
    double mass, mx, my, mz;      /* total mass and mass-weighted position */
    int    child[8];              /* child node indices, -1 = none         */
    int    body;                  /* body index if leaf with one body      */
    int    count;                 /* number of bodies inside this cell     */
} Node;

static Node *nodes = NULL;
static int n_nodes = 0, cap = 0;
static const double *g_pos, *g_m;

static int new_node(double cx, double cy, double cz, double half) {
    if (n_nodes == cap) {
        cap = cap ? cap * 2 : 1024;
        nodes = realloc(nodes, cap * sizeof(Node));
        if (!nodes) { fprintf(stderr, "bh: out of memory\n"); exit(1); }
    }
    Node *n = &nodes[n_nodes];
    n->cx = cx; n->cy = cy; n->cz = cz; n->half = half;
    n->mass = n->mx = n->my = n->mz = 0.0;
    for (int k = 0; k < 8; k++) n->child[k] = -1;
    n->body = -1;
    n->count = 0;
    return n_nodes++;
}

static void insert(int ni, int b, int depth);

/* put body b into the correct child of node ni (creating it if needed) */
static void insert_child(int ni, int b, int depth) {
    double x = g_pos[3*b], y = g_pos[3*b+1], z = g_pos[3*b+2];
    int o = (x > nodes[ni].cx) | ((y > nodes[ni].cy) << 1) | ((z > nodes[ni].cz) << 2);
    if (nodes[ni].child[o] < 0) {
        double h  = nodes[ni].half * 0.5;
        double cx = nodes[ni].cx + ((o & 1) ? h : -h);
        double cy = nodes[ni].cy + ((o & 2) ? h : -h);
        double cz = nodes[ni].cz + ((o & 4) ? h : -h);
        int c = new_node(cx, cy, cz, h);      /* may realloc: use indices only */
        nodes[ni].child[o] = c;
    }
    insert(nodes[ni].child[o], b, depth + 1);
}

static void insert(int ni, int b, int depth) {
    double mb = g_m[b];
    nodes[ni].mass += mb;
    nodes[ni].mx += mb * g_pos[3*b];
    nodes[ni].my += mb * g_pos[3*b+1];
    nodes[ni].mz += mb * g_pos[3*b+2];

    if (nodes[ni].count == 0) {               /* empty cell: becomes a leaf */
        nodes[ni].body = b;
        nodes[ni].count = 1;
        return;
    }
    if (depth >= MAXDEPTH) {                  /* safety: stop subdividing   */
        nodes[ni].count++;
        nodes[ni].body = -1;
        return;
    }
    if (nodes[ni].count == 1 && nodes[ni].body >= 0) {
        int old = nodes[ni].body;             /* push the old body down     */
        nodes[ni].body = -1;
        insert_child(ni, old, depth);
    }
    nodes[ni].count++;
    insert_child(ni, b, depth);
}

static void walk(int ni, int i, double xi, double yi, double zi,
                 double theta2, double eps2,
                 double *fx, double *fy, double *fz) {
    const Node *n = &nodes[ni];
    if (n->count == 0) return;
    if (n->body == i) return;                 /* never interact with itself */

    double dx = n->mx / n->mass - xi;
    double dy = n->my / n->mass - yi;
    double dz = n->mz / n->mass - zi;
    double r2 = dx*dx + dy*dy + dz*dz;
    double size = 2.0 * n->half;

    /* leaf, or cell far enough away: treat as ONE body at its centre of mass */
    int direct = (n->body >= 0) || (size * size < theta2 * r2);
    if (!direct) {
        int any = 0;
        for (int k = 0; k < 8; k++)
            if (n->child[k] >= 0) {
                any = 1;
                walk(n->child[k], i, xi, yi, zi, theta2, eps2, fx, fy, fz);
            }
        if (any) return;
    }
    r2 += eps2;
    double inv = n->mass / (r2 * sqrt(r2));
    *fx += dx * inv; *fy += dy * inv; *fz += dz * inv;
}

void bh_accel(int lo, int hi, int N, const double *pos, const double *m,
              double *acc, double theta) {
    g_pos = pos; g_m = m;

    /* bounding cube around all bodies */
    double mn[3] = {pos[0], pos[1], pos[2]};
    double mx[3] = {pos[0], pos[1], pos[2]};
    for (int i = 1; i < N; i++)
        for (int k = 0; k < 3; k++) {
            double v = pos[3*i+k];
            if (v < mn[k]) mn[k] = v;
            if (v > mx[k]) mx[k] = v;
        }
    double ext = 0;
    for (int k = 0; k < 3; k++) if (mx[k] - mn[k] > ext) ext = mx[k] - mn[k];
    double half = 0.5 * ext * 1.000001 + 1e-9;

    n_nodes = 0;
    int root = new_node(0.5*(mn[0]+mx[0]), 0.5*(mn[1]+mx[1]), 0.5*(mn[2]+mx[2]), half);
    for (int b = 0; b < N; b++) insert(root, b, 0);   /* same order on every rank */

    const double eps2 = EPS * EPS, theta2 = theta * theta;
    for (int i = lo; i < hi; i++) {
        double fx = 0, fy = 0, fz = 0;
        walk(root, i, pos[3*i], pos[3*i+1], pos[3*i+2], theta2, eps2, &fx, &fy, &fz);
        acc[3*i] = G * fx; acc[3*i+1] = G * fy; acc[3*i+2] = G * fz;
    }
}

void bh_free(void) { free(nodes); nodes = NULL; n_nodes = cap = 0; }
