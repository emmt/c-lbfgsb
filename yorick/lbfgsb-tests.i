/*
 * lbfgsb-tests.i --
 *
 * Tests for the Yorick interface to L-BFGS-B algorithm.
 *
 *-----------------------------------------------------------------------------
 *
 * This file is part of software licensed under the MIT license
 * (https://github.com/emmt/clbfgsb).
 *
 * Copyright (C) 2021, Éric Thiébaut.
 */

func lbfgsb_test1(nil, n=, m=)
/* DOCUMENT lbfgsb_test1, n=..., m=..., fg=...;

   This simple example demonstrates how to call the L-BFGS-B code to solve a
   simple problem (the extended Rosenbrock function subject to bounds on the
   variables).

   Keywords `n` and `m` are to specify the dimension `n` of this problem and
   the maximum number `m` of memorized steps.  By default, `n=25` and `m=5`.

   Keywords `fg` can be used to specify a different objective function
   than the extended Rosenbrock function `lbfgsb_test_fg`.

   SEE ALSO: lbfgsb_test_fg, lbfgsb_create, lbfgsb_config, lbfgsb_iterate,
             lbfgsb_reset, lbfgsb_stop,
 */
{
    // Default problem, problem size and maximum number of memorized steps.
    if (is_void(fg)) fg = lbfgsb_test_fg;
    if (is_void(n)) n = 25;
    if (is_void(m)) m =  5;

    // Create context.
    ctx = lbfgsb_create(n, m);

    // Choose bounds.
    lower = array(double, n);
    lower(1:n:2) =    1.0;
    lower(2:n:2) = -100.0;
    upper = 100.0;

    // Configure context.
    lbfgsb_config, ctx, lower = lower, upper = upper,
        // We wish to have output at every iteration.
        print = 1,
        // We specify the tolerances in the stopping criteria.
        factr = 1.0e+7,
        pgtol = 1.0e-5;

    // Allocate and initialize variables.
    x = array(3.0, n);

    // Variables to store function value and its gradient.
    f = LBFGSB_NAN; // initial value is irrelevant
    g = array(double, n);

    // Run algorithm.
    printf, "\n     %s\n      %s\n\n", "Solving sample problem.",
        "(f = 0.0 at the optimal solution.)";
    while (1) {
        // Iterate algorithm.
        task = lbfgsb_iterate(ctx, x, f, g);

        if (task == LBFGSB_FG) {
            // The minimization routine has requested the function f and
            // gradient g values at the current x.
            f = fg(x, g);
            continue;
        }

        if (task == LBFGSB_NEW_X) {
            // A new iterate is available for inspection.
            continue;
        }

        // Convergence or error.
        break;
    }
}

func lbfgsb_test2(nil, n=, m=, fg=)
/* DOCUMENT lbfgsb_test2, n=..., m=..., fg=...;

   This example shows how to replace the default stopping test by other
   termination criteria. It also illustrates how to print the values of several
   parameters during the course of the iteration.  The test problem used here
   is the same as in `clbfgsb_test1` (the extended Rosenbrock function with
   bounds on the variables).

   Keywords `n` and `m` are to specify the dimension `n` of this problem and
   the maximum number `m` of memorized steps.  By default, `n=25` and `m=5`.

   Keywords `fg` can be used to specify a different objective function
   than the extended Rosenbrock function `lbfgsb_test_fg`.

   SEE ALSO: lbfgsb_test_fg, lbfgsb_create, lbfgsb_config, lbfgsb_iterate,
             lbfgsb_reset, lbfgsb_stop,
 */
{
    // Default problem, problem size and maximum number of memorized steps.
    if (is_void(fg)) fg = lbfgsb_test_fg;
    if (is_void(n)) n = 25;
    if (is_void(m)) m =  5;

    // Choose bounds.
    lower = array(double, n);
    lower(1:n:2) =    1.0;
    lower(2:n:2) = -100.0;
    upper = 100.0;

    // Create and configure context.
    ctx = lbfgsb_config(
        lbfgsb_create(n, m),
        // Specifiy the bounds.
        lower = lower,
        upper = upper,
        // We suppress the default output.
        print = -1,
        // We suppress both code-supplied stopping tests because the user is
        // providing his own stopping criteria.
        factr = 0.0,
        pgtol = 0.0);

    // Allocate and initialize variables.
    x = array(3.0, n);

    // Variables to store function value and its gradient.
    f = 0.0; // initial value is irrelevant
    g = array(double, n);

    // Run algorithm.
    printf, "\n     %s\n      %s\n\n", "Solving sample problem.",
        "(f = 0.0 at the optimal solution.)";
    while (1) {
        // Iterate algorithm.
        lbfgsb_iterate, ctx, x, f, g;

        if (ctx.task == LBFGSB_FG) {
            // The minimization routine has requested the function f and
            // gradient g values at the current x.
            f = fg(x, g);
            continue;
        }

        if (ctx.task == LBFGSB_NEW_X) {
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
            if (ctx.nevals >= 99) {
                lbfgsb_stop, ctx,
                    "STOP: Total no. of f and g evaluations exceeds limit";
            }

            // 2) Terminate if  |proj g|/(1+|f|) < 1.0e-10, where
            // "proj g" denoted the projected gradient
            if (ctx.pgnorm <= 1.0e-10*(1.0 + abs(f))) {
                lbfgsb_stop, ctx,
                    "STOP: the projected gradient is sufficiently small";
            }

            // We now wish to print the following information at each
            // iteration:
            //
            // 1) the current iteration number,
            // 2) the total number of f and g evaluations,
            // 3) the value of the objective function f,
            // 4) the norm of the projected gradient,
            printf,
                "Iterate%5d    nfg =%5d    f =%12.5E    |proj g| =%12.5E\n",
                ctx.niters, ctx.nevals, f, ctx.pgnorm;

            if (ctx.task == LBFGSB_NEW_X) {
                continue;
            }
        }

        // The run is to be terminated, we print also the information
        // contained in task as well as the final value of x.
        printf, " %s\n Final X=\n", ctx.reason;
        lbfgsb_print_variables, x;
        break;
    }
}

func lbfgsb_test3(nil, n=, m=, fg=, maxtime=, maxeval=, maxiter=)
/* DOCUMENT lbfgsb_test3, n=..., m=..., fg=...,
                maxtime=..., maxeval=..., maxiter=...;

   This time-controlled example shows that it is possible to terminate a run by
   elapsed CPU time, and yet be able to print all desired information. This
   driver also illustrates the use of two stopping criteria that may be used in
   conjunction with a limit on execution time. The test problem used here is
   the same as in `clbfgsb_test1.c` and `clbfgsb_test2.c` (the extended
   Rosenbrock function with bounds on the variables).

   Keywords `n` and `m` are to specify the dimension `n` of this problem and
   the maximum number `m` of memorized steps.  By default, `n=1000` and `m=10`.

   Keywords `fg` can be used to specify a different objective function
   than the extended Rosenbrock function `lbfgsb_test_fg`.

   Keyword `maxtime` is to specify the time limit in seconds (default 0.2
   seconds).

   Keyword `maxiter` is to specify the maximum number of iterations (no limits
   by default).

   Keyword `maxeval` is to specify the maximum number of objective function
   evaluations (no limits by default).

   SEE ALSO: lbfgsb_test_fg, lbfgsb_create, lbfgsb_config, lbfgsb_iterate,
             lbfgsb_reset, lbfgsb_stop,
 */
{
    // Default problem, problem size and maximum number of memorized steps.
    if (is_void(fg)) fg = lbfgsb_test_fg;
    if (is_void(n)) n = 1000;
    if (is_void(m)) m =   10;

    // Other settings.
    Inf = LBFGSB_INFINITY;
    if (is_void(maxtime)) maxtime = 0.2;
    if (is_void(maxiter)) maxiter = Inf;
    if (is_void(maxeval)) maxeval = Inf;

    // Choose bounds.
    lower = array(double, n);
    lower(1:n:2) =    1.0;
    lower(2:n:2) = -100.0;
    upper = 100.0;

    // Create and configure context.
    ctx = lbfgsb_config(
        lbfgsb_create(n, m),
        // Specifiy the bounds.
        lower = lower,
        upper = upper,
        // We suppress the default output.
        print = -1,
        // We suppress both code-supplied stopping tests because the user is
        // providing his own stopping criteria.
        factr = 0.0,
        pgtol = 0.0);

    // Allocate and initialize variables.
    x = array(3.0, n);

    // Variables to store function value and its gradient.
    f = 0.0; // initial value is irrelevant
    g = array(double, n);

    // Initial time.
    ts = array(double, 3);
    timer, ts;
    t0 = ts(1); // CPU time

    // Run algorithm.
    printf, "\n     %s\n      %s\n\n", "Solving sample problem.",
        "(f = 0.0 at the optimal solution.)";
    while (1) {
        // Iterate algorithm.
        task = lbfgsb_iterate(ctx, x, f, g);

        if (task == LBFGSB_FG) {
            // The minimization routine has requested the function f and
            // gradient g values at the current x.

            // Before evaluating f and g we check the CPU time spent.
            timer, ts;
            t = ts(1) - t0; // CPU time
            if (t > maxtime) {
                task = lbfgsb_stop(
                    ctx, "STOP: CPU exceeding the time limit");
            }

            // Before evaluating f and g we check the number of evaluations
            // limit.
            if (task == LBFGSB_FG && ctx.nevals >= maxeval) {
                task = lbfgsb_stop(
                    ctx, "STOP: Too many function evaluations");
            }

            // If no limits have been overreached, we compute
            // the function value f for problem.
            if (task == LBFGSB_FG) {
                f = fg(x, g);
                continue;
            }
        }

        if (task == LBFGSB_NEW_X) {
            // The minimization routine has returned with a new iterate.  The
            // time limit has not been reached, and we test whether the
            // following two stopping tests are satisfied:
            //
            // 1) Terminate if the total number of f and g iterations
            //    exceeds `maxiter`.
            if (ctx.niters >= maxiter) {
                task = lbfgsb_stop(
                    ctx,
                    "STOP: Too many iterations");
            }

            // 2) Terminate if `|proj g|/(1+|f|) < 1.0e-10`.
            if (ctx.pgnorm <= 1.0e-10*(1.0 + abs(f))) {
                task = lbfgsb_stop(
                    ctx,
                    "STOP: The projected gradient is sufficiently small");
            }

            // We now wish to print the following information at each
            // iteration:
            //
            // 1) the current iteration number,
            // 2) the total number of f and g evaluations,
            // 3) the value of the objective function f,
            // 4) the norm of the projected gradient,
            printf,
                "Iterate%5d    nfg =%5d    f =%12.5E    |proj g| =%12.5E\n",
                ctx.niters, ctx.nevals, f, ctx.pgnorm;

            if (ctx.task == LBFGSB_NEW_X) {
                continue;
            }
        }

        // The run is to be terminated, we print also the information
        // contained in task as well as the final value of x.
        printf, " %s\n Final X=\n", ctx.reason;
        lbfgsb_print_variables, x;
        break;
    }
}

func lbfgsb_test_fg(x, &g)
/* DOCUMENT f = lbfgsb_test_fg(x, g);

     Compute the extended Rosenbrock function at `x` and stores its gradient
     in `g`.

   SEE ALSO:
 */
{
    n = numberof(x);
    T = structof(x);
    if (T != float && T != double) {
        T = double;
        x = T(x);
    }
    if (numberof(g) != n || structof(g) != T) {
        g = array(T, dimsof(x));
    }
    s = x(1:-1);
    t = x(2:0) - s*s;

    // Compute gradient g for the sample problem.
    g(1) = (2 - 16*t(1))*x(1) - 2;
    if (n > 2) {
        g(2:-1) = 8*t(1:-1) - 16*x(2:-1)*t(2:0);
    }
    g(0) = 8*t(0);

    // Compute function value f for the sample problem.
    r = x(1) - 1;
    f = r*r + 4*sum(t*t);
    return f;
}

func lbfgsb_print_variables(x)
{
    n = numberof(x);
    for (i = 1; i <= n; ++i) {
        printf, "%s%12.4E%s",
            ((i%6) == 1 ? " " : ""), x(i),
            (i == n || (i%6) == 0 ? "\n" : "");
    }
}

local sprintf, printf, fprintf;
/* DOCUMENT printf, fmt, a1, a2,  ...;
         or fprintf, io, fmt, a1, a2,  ...;
         or str = sprintf(fmt, a1, a2, ...);

     Use format string `fmt` to print arguments `a1`, `a2`, etc. to standard
     or `io` output stream or to build a string `str`.

   SEE ALSO: write.
 */
func printf(fmt, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10)
{
    if (!(is_string(fmt) && is_scalar(fmt))) {
        error, "usage: printf, fmt, a1, a2, ...;";
    }
    if (a1 == []) {
        write, format="%s", fmt;
    } else if (a2 == []) {
        write, format=fmt, a1;
    } else if (a3 == []) {
        write, format=fmt, a1, a2;
    } else if (a4 == []) {
        write, format=fmt, a1, a2, a3;
    } else if (a5 == []) {
        write, format=fmt, a1, a2, a3, a4;
    } else if (a6 == []) {
        write, format=fmt, a1, a2, a3, a4, a5;
    } else if (a7 == []) {
        write, format=fmt, a1, a2, a3, a4, a5, a6;
    } else if (a8 == []) {
        write, format=fmt, a1, a2, a3, a4, a5, a6, a7;
    } else if (a9 == []) {
        write, format=fmt, a1, a2, a3, a4, a5, a6, a7, a8;
    } else if (a10 == []) {
        write, format=fmt, a1, a2, a3, a4, a5, a6, a7, a8, a9;
    } else {
        write, format=fmt, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10;
    }
}

func fprintf(io, fmt, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10)
{
    if (!(is_string(fmt) && is_scalar(fmt))) {
        error, "usage: fprintf, io, fmt, a1, a2, ...;";
    }
    if (a1 == []) {
        write, io, format="%s", fmt;
    } else if (a2 == []) {
        write, io, format=fmt, a1;
    } else if (a3 == []) {
        write, io, format=fmt, a1, a2;
    } else if (a4 == []) {
        write, io, format=fmt, a1, a2, a3;
    } else if (a5 == []) {
        write, io, format=fmt, a1, a2, a3, a4;
    } else if (a6 == []) {
        write, io, format=fmt, a1, a2, a3, a4, a5;
    } else if (a7 == []) {
        write, io, format=fmt, a1, a2, a3, a4, a5, a6;
    } else if (a8 == []) {
        write, io, format=fmt, a1, a2, a3, a4, a5, a6, a7;
    } else if (a9 == []) {
        write, io, format=fmt, a1, a2, a3, a4, a5, a6, a7, a8;
    } else if (a10 == []) {
        write, io, format=fmt, a1, a2, a3, a4, a5, a6, a7, a8, a9;
    } else {
        write, io, format=fmt, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10;
    }
}

func sprintf(fmt, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10)
{
    if (!(is_string(fmt) && is_scalar(fmt))) {
        error, "usage: sprintf(fmt, a1, a2, ...);";
    }
    if (a1 == []) {
        return fmt;
    } else if (a2 == []) {
        return swrite(format=fmt, a1);
    } else if (a3 == []) {
        return swrite(format=fmt, a1, a2);
    } else if (a4 == []) {
        return swrite(format=fmt, a1, a2, a3);
    } else if (a5 == []) {
        return swrite(format=fmt, a1, a2, a3, a4);
    } else if (a6 == []) {
        return swrite(format=fmt, a1, a2, a3, a4, a5);
    } else if (a7 == []) {
        return swrite(format=fmt, a1, a2, a3, a4, a5, a6);
    } else if (a8 == []) {
        return swrite(format=fmt, a1, a2, a3, a4, a5, a6, a7);
    } else if (a9 == []) {
        return swrite(format=fmt, a1, a2, a3, a4, a5, a6, a7, a8);
    } else if (a10 == []) {
        return swrite(format=fmt, a1, a2, a3, a4, a5, a6, a7, a8, a9);
    } else {
        return swrite(format=fmt, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10);
    }
}
