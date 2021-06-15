// clbfgsb_test2.c -
//
// This example shows how to replace the default stopping test by other
// termination criteria. It also illustrates how to print the values of several
// parameters during the course of the iteration.  The test problem used here
// is the same as in `clbfgsb_test1.c` (the extended Rosenbrock function with
// bounds on the variables).
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

static void print_variables(
    const double x[],
    long         n)
{
    for (long i = 0; i < n; ++i) {
        printf("%s%12.4E%s", ((i%6) == 0 ? " " : ""), x[i],
               (i == n-1 || (i%6) == 5 ? "\n" : ""));
    }
}

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

    // We suppress the default output.
    ctx->print = -1;

    // We suppress both code-supplied stopping tests because the user is
    // providing his own stopping criteria.
    ctx->factr = 0.0;
    ctx->pgtol = 0.0;

    // Initialize bounds.
    double* lower = ctx->lower;
    double* upper = ctx->upper;
    for (long i = 0; i < n; ++i) {
        lower[i] = (i&1) == 0 ? 1.0 : -1.0e2;
        upper[i] = 1.0e2;
    }

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
            // The minimization routine has returned with a new iterate.  At
            // this point we have the opportunity of stopping the iteration or
            // observing the values of certain parameters.
            //
            // First are two examples of stopping tests.

            // Note: task must be assigned the value "STOP..." to terminate the
            // iteration and ensure that the final results are printed in the
            // default format. The rest of the character string task may be
            // used to store other information.
            //
            // 1) Terminate if the total number of f and g evaluations
            //    exceeds 99.
            if (LBFGSB_NTOT_FG(ctx) >= 99) {
                task = lbfgsb_set_task(
                    ctx,
                    "STOP: TOTAL NO. of f AND g EVALUATIONS EXCEEDS LIMIT");
            }

            // 2) Terminate if  |proj g|/(1+|f|) < 1.0e-10, where
            // "proj g" denoted the projected gradient
            if (LBFGSB_PG_NORMINF(ctx) <= 1.0e-10*(1.0 + fabs(f))) {
                task = lbfgsb_set_task(
                    ctx,
                    "STOP: THE PROJECTED GRADIENT IS SUFFICIENTLY SMALL");
            }

            // We now wish to print the following information at each
            // iteration:
            //
            // 1) the current iteration number,
            // 2) the total number of f and g evaluations,
            // 3) the value of the objective function f,
            // 4) the norm of the projected gradient,
            printf("Iterate %4d    nfg = %4d    "
                   "f =%12.5E    |proj g| =%12.5E\n",
                   LBFGSB_NUM_ITER(ctx), LBFGSB_NTOT_FG(ctx),
                   f, LBFGSB_PG_NORMINF(ctx));

            // If the run is to be terminated, we print also the information
            // contained in task as well as the final value of x.
            if (task == LBFGSB_STOP) {
                char buf[LBFGSB_TASK_LENGTH+1];
                printf(" %s\n Final X=\n",
                       lbfgsb_get_task_string(ctx, buf, sizeof(buf)));
                print_variables(x, n);
                break;
            }
            continue;
        }

        // Convergence or error.
        break;
    }

    // Release resources.
    lbfgsb_destroy(ctx);
    return EXIT_SUCCESS;
}
