#pragma once

/* Conjugate gradients solver for symmetric positive-definite systems.
 *
 * All storage is plain double arrays; no external library dependencies.
 * Designed for repeated solves against the same or slowly-changing operator,
 * with pre-allocated workspace to avoid malloc overhead at n ~ 500k. */

/* Callback: compute out = A * in.  n is the vector length; data is opaque. */
typedef void (*cg_matvec_fn)(double *out, const double *in, int n, void *data);

typedef struct {
    int    max_iter;      /* 0 → use n as the limit */
    double tol;           /* convergence: ||r||/||b|| < tol.  Default 1e-6 */
    double abs_tol;       /* also stop if ||r|| < abs_tol.  0 → disabled */
    double *diag_precond; /* optional length-n diagonal preconditioner (entries = 1/d_i);
                             NULL → no preconditioning */
} cg_opts_t;

typedef struct {
    int     n;
    double *r;            /* residual */
    double *p;            /* search direction */
    double *q;            /* A * p */
    double *z;            /* preconditioned residual (only used when diag_precond != NULL) */
    cg_opts_t opts;
    /* diagnostics from most recent solve */
    int    last_iters;
    double last_resid;    /* final ||r|| / ||b|| */
} cg_state_t;

/* Allocate state and workspace for vectors of length n.
 * opts may be NULL for defaults (tol=1e-6, max_iter=n, no preconditioner).
 * Returns NULL on allocation failure. */
cg_state_t *cg_state_new(int n, const cg_opts_t *opts);

void cg_state_free(cg_state_t *st);

/* Update options on an existing state (e.g. to change tol between solves).
 * Does NOT reallocate workspace; diag_precond pointer is copied by reference. */
void cg_state_set_opts(cg_state_t *st, const cg_opts_t *opts);

/* Solve A x = b via (preconditioned) conjugate gradients.
 *
 * x[n]: initial guess on entry, approximate solution on exit.
 *   - Zero x before calling for a cold start.
 *   - Leave x as the previous solution for a warm start.
 *
 * Returns the number of iterations used, or -1 on error (e.g. b = 0). */
int cg_solve(cg_state_t *st, cg_matvec_fn Av, void *data,
             const double *b, double *x);
