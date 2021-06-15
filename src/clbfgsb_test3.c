// clbfgsb_test3.c -
//
// This time-controlled example shows that it is possible to terminate a run by
// elapsed CPU time, and yet be able to print all desired information. This
// driver also illustrates the use of two stopping criteria that may be used in
// conjunction with a limit on execution time. The test problem used here is
// the same as in `clbfgsb_test1.c` and `clbfgsb_test2.c` (the extended
// Rosenbrock function with bounds on the variables).
//
// The dimension `n` of this problem and/or the maximum number `m` of steps to
// memorize can be set by compiling with `-DN=...` and/or `-DM=...`.
//
// The time limit (in seconds) can be set by compiling with `-DTLIMIT=...`.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lbfgsb.h>

// Number of variables.
#ifndef N
# define N 1000
#endif

// Number of steps to memorize.
#ifndef M
# define M 10
#endif

// Limit on the CPU time (in seconds).
#ifndef TLIMIT
# define TLIMIT 0.2
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

    // We suppress both code-supplied stopping tests because we will provide
    // our own termination conditions.
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

    // We begin counting the CPU time.
    double t0 = lbfgsb_timer();

    // Run algorithm.
    printf("\n     %s\n      %s\n\n",
           "Solving sample problem.",
           "(f = 0.0 at the optimal solution.)");
    while (1) {
        int task = lbfgsb_iterate(ctx, x, &f, g);
        if (task == LBFGSB_FG) {
            // The minimization routine has requested the function f and
            // gradient g values at the current x.

            // Before evaluating f and g we check the CPU time spent.
            double t = lbfgsb_timer() - t0;
            if (t > TLIMIT) {
                // Note: Assigning `task[0..3] = "STOP"` will terminate the
                // run; setting `task[6..9] = "CPU"` will restore the
                // information at the latest iterate generated by the code so
                // that it can be correctly printed by this example.
                task = lbfgsb_set_task(
                    ctx,
                    "STOP: CPU EXCEEDING THE TIME LIMIT.");

                // In this example we have chosen to disable the printing
                // options of the code (we set `print = -1`); instead we are
                // using customized output: we print the latest value of `x`,
                // the corresponding function value `f` and the norm of the
                // projected gradient `|proj g|`.

                // We print out the information contained in task.
                char buf[LBFGSB_TASK_LENGTH+1];
                printf(" %s\n Latest X=\n",
                       lbfgsb_get_task_string(ctx, buf, sizeof(buf)));
                print_variables(lbfgsb_get_latest_x(ctx), n);

                // We print the function value `f` and the norm of the
                // projected gradient `|proj g|` at the last iterate.
                printf(" At latest iterate   f =%12.5E    |proj g| =%12.5E\n",
                       LBFGSB_PREV_F(ctx), LBFGSB_PG_NORMINF(ctx));
            } else {
                // The time limit has not been reached and we compute
                // the function value f for the sample problem.
                f = compute_fg(x, g, n);
                continue;
            }
        }

        if (task == LBFGSB_NEW_X) {
            // The minimization routine has returned with a new iterate.  The
            // time limit has not been reached, and we test whether the
            // following two stopping tests are satisfied:
            //
            // 1) Terminate if the total number of f and g evaluations
            //    exceeds 900.
            if (LBFGSB_NTOT_FG(ctx) >= 900) {
                task = lbfgsb_set_task(
                    ctx,
                    "STOP: TOTAL NO. of f AND g EVALUATIONS EXCEEDS LIMIT");
            }

            // 2) Terminate if `|proj g|/(1+|f|) < 1.0e-10`.
            if (LBFGSB_PG_NORMINF(ctx) <= 1.0e-10*(1.0 + fabs(f))) {
                task = lbfgsb_set_task(
                    ctx,
                    "STOP: THE PROJECTED GRADIENT IS SUFFICIENTLY SMALL");
            }

            // We wish to print the following information at each iteration:
            //
            // 1) the current iteration number,
            // 2) the total number of f and g evaluations,
            // 3) the value of the objective function f,
            // 4) the norm of the projected gradient,
            printf("Iterate %4d    nfg = %4d    "
                   "f =%12.5E    |proj g| =%12.5E\n",
                   LBFGSB_NUM_ITER(ctx), LBFGSB_NTOT_FG(ctx),
                   f, LBFGSB_PG_NORMINF(ctx));
            continue;
        }

        if (task == LBFGSB_STOP) {
            // If the run is to be terminated, we print also the information
            // contained in task as well as the final value of x.
            char buf[LBFGSB_TASK_LENGTH+1];
            printf(" %s\n Final X=\n",
                   lbfgsb_get_task_string(ctx, buf, sizeof(buf)));
            print_variables(lbfgsb_get_latest_x(ctx), n);
            break;
        }

        // Convergence or error.
        break;
    }

    // Release resources.
    lbfgsb_destroy(ctx);
    return EXIT_SUCCESS;
}