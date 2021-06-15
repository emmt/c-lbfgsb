// clbfgsb_test1.c -
//
// This simple example demonstrates how to call the L-BFGS-B code to solve a
// simple problem (the extended Rosenbrock function subject to bounds on the
// variables).
//
// The dimension `N` of this problem and/or the maximum number `M` of steps to
// memorize can be set by compiling with `-DN=...` and/or `-DM=...`.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lbfgsb.h>

// Number of variables.
#ifndef N
# define N 25
#endif

// Number of steps to memorize.
#ifndef M
# define M 5
#endif

static inline double pow2(double x) { return x*x; }

static double compute_fg(
    const double x[],
    double       g[],
    long         n)
{
    // Compute function value f for the sample problem.
    double f = pow2(x[0] - 1);
    for (long i = 1; i < n; ++i) {
        f += 4*pow2(x[i] - pow2(x[i-1]));
    }

    // Compute gradient g for the sample problem.
    double t1 = x[1] - pow2(x[0]);
    g[0] = 2*(x[0] - 1) - 16*x[0]*t1;
    for (long i = 1; i < n-1; ++i) {
        double t2 = t1;
        t1 = x[i+1] - pow2(x[i]);
        g[i] = 8*t2 - 16*x[i]*t1;
    }
    g[n-1] = 8*t1;

    return f;
}

int main(int argc, char* argv[])
{
    // Problem size and maximum number of memorized steps.
    long n = N, m = M;

     // Create new context.
    lbfgsb_context* ctx = lbfgsb_create(n, m);
    if (ctx == NULL) {
        fprintf(stderr, "failed to allocate context\n");
        return EXIT_FAILURE;
    }

    // Initialize bounds.
    double* lower = ctx->lower;
    double* upper = ctx->upper;
    for (long i = 0; i < n; ++i) {
        lower[i] = (i&1) == 0 ? 1.0 : -1.0e2;
        upper[i] = 1.0e2;
    }

    // We wish to have output at every iteration.
    ctx->print = 1;

    // We specify the tolerances in the stopping criteria.
    ctx->factr = 1.0e+7;
    ctx->pgtol = 1.0e-5;

    // Allocate and initialize variables.
    double x[n];
    for (long i = 0; i < n; ++i) {
        x[i] = 3.0;
    }

    // Variables to store function value and its gradient.
    double f = LBFGSB_NAN; // initial value is irrelevant
    double g[n];

    // Run algorithm.
    printf("\n     %s\n      %s\n\n",
           "Solving sample problem.",
           "(f = 0.0 at the optimal solution.)");
    while (1) {
        // Iterate algorithm.
        int task = lbfgsb_iterate(ctx, x, &f, g);

        if (task == LBFGSB_FG) {
            // The minimization routine has requested the function f and
            // gradient g values at the current x.
            f = compute_fg(x, g, n);
            continue;
        }

        if (task == LBFGSB_NEW_X) {
            // A new iterate is available for inspection.
            continue;
        }

        // Convergence or error.
        break;
    }

    // Release resources.
    lbfgsb_destroy(ctx);
    return EXIT_SUCCESS;
}
