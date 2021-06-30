/*
 * lbfgsb.i --
 *
 * Yorick interface to L-BFGS-B algorithm.
 *
 *-----------------------------------------------------------------------------
 *
 * This file is part of software licensed under the MIT license
 * (https://github.com/emmt/clbfgsb).
 *
 * Copyright (C) 2021, Éric Thiébaut.
 */

if (is_func(plug_in)) plug_in, "ylbfgsb";

func lbfgsb(fg, x0, &f, &g, &status, lower=, upper=, mem=, cputime=,
            ftol=, gtol=, xtol=, maxiter=, maxeval=, verb=, output=)
/* DOCUMENT x = lbfgsb(fg, x0, [f, g, status,] lower=, upper=, mem=);

     Apply L-BFGS-B algorithm to minimize a multi-variate differentiable
     objective function possibly under separable bound constraints.

     The method has two required arguments: `fg`, the function to call to
     compute the objective function and its gradient, and `x0`, the initial
     variables (VMLMB is an iterative method).  The initial variables may be an
     array of any dimensions.

     The method returns `x` the best solution found during iterations.
     Arguments `f`, `g` and `status` are optional output variables to store the
     value and the gradient of the objective at `x` and an integer code
     indicating the reason of the termination of the algorithm.

     The function `fg` shall be implemented as follows:

         func fg(x, &gx)
         {
             fx = ...; // value of the objective function at `x`
             gx = ...; // gradient of the objective function at `x`
             return fx;
         }

     All other settings are specified by keywords:

     - Keywords `upper` and `lower` are to specify a lower and/or an upper
       bounds for the variables.  If unspecified or set to an empty array, a
       given bound is considered as unlimited.  Bounds must be conformable with
       the variables.

     - Keyword `mem` specifies the memory used by the algorithm, that is the
       number of previous steps memorized to approximate the Hessian of the
       objective function.  With `mem=0`, the algorithm behaves as a steepest
       descent method.  The default is `mem=5`.

     - Keywords `ftol`, `gtol` and `xtol` specify tolerances for deciding the
       convergence of the algorithm.

       Convergence in the function occurs if one of the following conditions
       hold:

           f ≤ fatol
           |f - fp| ≤ frtol⋅max(|f|, |fp|)

       where `f` and `fp` are the values of the objective function at the
       current and previous iterates.  In these conditions, `fatol` and `frtol`
       are absolute and relative tolerances specified by `ftol` which can be
       `ftol=[fatol,frtol]` or `ftol=frtol` and assume that `fatol=-Inf`.  The
       default is `ftol=1e-8`.

       Convergence in the gradient occurs if the following condition holds:

           ‖g‖ ≤ max(0, gatol, grtol⋅‖g0‖)

       where `‖g‖` is the Euclidean norm of the projected gradient, `g0` is the
       projected gradient at the initial solution.  In this condition, `gatol`
       and `grtol` are absolute and relative gradient tolerances specified by
       `gtol` which can be `gtol=[gatol,grtol]` or `gtol=grtol` and assume that
       `gatol=0`.  The default is `gtol=1e-5`.

       Convergence in the variables occurs if the following condition holds:

           ‖x - xp‖ ≤ max(0, xatol, xrtol*‖x‖)

       where `x` and `xp` are the current and previous variables.  In this
       condition, `xatol` and `xrtol` are absolute and relative tolerances
       specified by `xtol` which can be `xtol=[fatol,frtol]` or `xtol=xrtol`
       and assume that `xatol=0`.  The default is `xtol=1e-6`.

     - Keywords `maxiter` and `maxeval` are to specify a maximum number of
       algorithm iterations or or evaluations of the objective function
       implemented by `fg`.  By default, these are unlimited.

     - Keyword `verb`, if positive, specifies to print information every `verb`
       iterations.  Nothing is printed if `verb ≤ 0`.  By default, `verb = 0`.

     - Keyword `output` specifies the file stream to print information.

     - If keyword `cputime` is true, print the CPU time instead of the WALL
       time.

   SEE ALSO: lbfgsb_config, lbfgsb_create, lbfgsb_iterate, lbfgsb_stop.
 */
{
    // Constants.
    INF = LBFGSB_INFINITY;
    NAN = LBFGSB_NAN;
    TRUE = 1n;
    FALSE = 0n;

    // Parse settings.
    if (is_void(mem)) mem = 5;
    if (is_void(maxiter)) maxiter = INF;
    if (is_void(maxeval)) maxeval = INF;
    if (is_void(ftol)) ftol = 1.0E-8;
    if (is_void(gtol)) gtol = 1.0E-5;
    if (is_void(xtol)) xtol = 1.0E-6;
    if (is_void(verb)) verb = FALSE;

    // Tolerances.  Most of these are forced to be nonnegative to simplify
    // tests.
    if (is_scalar(ftol)) {
        fatol = -INF;
        frtol = max(0.0, ftol);
    } else {
        fatol = max(0.0, ftol(1));
        frtol = max(0.0, ftol(2));
    }
    if (is_scalar(gtol)) {
        gatol = 0.0;
        grtol = max(0.0, gtol);
    } else {
        gatol = max(0.0, gtol(1));
        grtol = max(0.0, gtol(2));
    }
    if (is_scalar(xtol)) {
        xatol = 0.0;
        xrtol = max(0.0, xtol);
    } else {
        xatol = max(0.0, xtol(1));
        xrtol = max(0.0, xtol(2));
    }

    // Bound constraints.  For faster code, unlimited bounds are preferentially
    // represented by empty arrays.
    if (is_array(lower) && allof(lower == -INF)) lower = [];
    if (is_array(upper) && allof(upper == +INF)) upper = [];

    // Create and configure context.
    ctx = lbfgsb_create(dimsof(x0), mem);
    lbfgsb_config,
        ctx,
        // Specify the bounds.
        lower = lower,
        upper = upper,
        // Suppress the default output and code-supplied stopping tests.
        print = -1,
        factr = 0.0,
        pgtol = 0.0;

    // Other initialization.
    x = double(x0);// initial iterate (forcing copy)
    f = 0.0;       // objective function
    g = array(double, dimsof(x)); // gradient
    best_f = INF;  // best function value so far
    best_g = [];   // corresponding gradient
    best_x = [];   // corresponding variables
    evals = 0;     // number of calls to fg
    iters = 0;     // number of iterations
    if (verb > 0) {
        time_index = (cputime ? 1 : 3);
        elapsed = array(double, 3);
        timer, elapsed;
        t0 = elapsed(time_index);
    }
    f0 = f;
    gtest = [];
    xtest = (xatol > 0 || xrtol > 0);
    while (TRUE) {
        task = lbfgsb_iterate(ctx, x, f, g);
        if (task == LBFGSB_FG) {
            if (evals >= maxeval) {
                task = lbfgsb_stop(
                    ctx, "WARNING: Too many function evaluations");
            } else {
                f = fg(x, g);
                ++evals;
                if (f < best_f) {
                    best_f = f;
                    best_g = g;
                    best_x = x;
                }
                if (evals == 1 && verb > 0) {
                    write, output, format="%s%s\n%s%s\n",
                        "# Iter.   Time (ms)   Eval.   Skips ",
                        "       Obj. Func.           Grad.       Step",
                        "# ----------------------------------",
                        "-----------------------------------------------";
                    gnorm = lbfgsb_pgnorm2(ctx, x, g);
                    _lbfgsb_print;
                }
                continue;
            }
        }
        if (task == LBFGSB_NEW_X) {
            gnorm = lbfgsb_pgnorm2(ctx, x, g);
            if (gtest == []) {
                gtest =  max(gatol, grtol*gnorm);
            }
            if (gnorm <= gtest) {
                task = lbfgsb_stop(
                    ctx, "STOP: ‖∇f(x)‖ ≤ max(gatol, grtol⋅‖∇f(x0)‖)");
            } else if (f <= fatol) {
                task = lbfgsb_stop(
                    ctx, "STOP: f(x) ≤ fatol");
            } else if (iters > 0 && abs(f - f0) <= frtol*max(abs(f), abs(f0))) {
                task = lbfgsb_stop(
                    ctx, "STOP: |Δf(x)| ≤ frtol⋅|f(x)|");
            } else if (xtest) {
                s = x - x0;
                snorm = sqrt(sum(s*s));
                s = [];
                if (snorm <= xatol ||
                    xrtol > 0 && snorm <= xrtol*sqrt(sum(x*x))) {
                    task = lbfgsb_stop(
                        ctx, "STOP: ‖Δx| ≤ max(xatol, xrtol⋅‖x‖)");
                }
            }
            ++iters;
            if (iters >= maxiter) {
                task = lbfgsb_stop(
                    ctx, "WARNING: Too many algorithm iterations");
            }
        }
        if (verb > 0 && ((iters % verb) == 0 || task != LBFGSB_NEW_X)) {
            _lbfgsb_print;
        }
        if (task != LBFGSB_NEW_X) {
            break;
        }
        if (xtest && evals > 1) {
            x0 = x;
        }
    }

    // Restore best solution so far and return solution (and status).
    if (best_f < f) {
        f = best_f;
        eq_nocopy, g, best_g;
        eq_nocopy, x, best_x;
    }
    if (verb > 0) {
        write, output, format="# Termination: %s\n", ctx.reason;
    }
    return x;
}

// Print LBFGSB iteration, all parameters are external.
func _lbfgsb_print
{
    timer, elapsed;
    t = (elapsed(time_index) - t0)*1E3; // elapsed milliseconds
    skips = ctx.nskips;
    alpha = (iters < 1 ? 0.0 : ctx.step);
    write, output,
        format="%7d %11.3f %7d %7d %23.15e %11.3e %11.3e\n",
        iters, t, evals, skips, f, gnorm, alpha;
}

extern lbfgsb_create;
/* DOCUMENT ctx = lbfgsb_create(dims, mem);

     Creates a new context for solving a box constrained minimization problem
     with L-BFGS-B algorithm.  Argument `dims` specifies the dimension list of
     the variables and argument `mem` specifies the maximum number of steps
     memorized by the algorithm.

     The L-BFGS-B context has a number of members:

     - `ctx.siz`: the number of variables;

     - `ctx.dims`: the dimensions of the variables;

     - `ctx.mem`: the maximum number of memorized steps;

     - `ctx.epsmch`: the computed machine precision (set after the first call
       to `lbfgsb_iterate`);

     - `ctx.factr`: the relative function decrease threshold (see
       `lbfgsb_config`);

     - `ctx.pgnorm`: the infinite norm of the projected gradient;

     - `ctx.pgtol`: the threshold on the infinite norm of the projected
       gradient (see `lbfgsb_config`);

     - `ctx.niters`: the number of iterations of the algorithm;

     - `ctx.nevals`: the number of computations of the objective function and
       its gradient;

     - `ctx.nskips`: total number of skipped BFGS updates before the current
       iteration;

     - `ctx.lower`: the lower bounds, initially all set to `-Inf`;

     - `ctx.upper`: the upper bounds, initially all set to `+Inf`;

     - `ctx.print`: the level of verbosity (negative to suppress all printed
       information);

     - `ctx.task`: the current stage of the algorithm;

     - `ctx.reason`: a detailed textual information of the current stage of the
       algorithm.

     - `ctx.step`: the relative length of the step along the search direction.

     - `ctx.theta`: the scaling parameter of the BFGS matrix;

     The bounds of the problem `ctx.lower` and `ctx.upper` and parameters
     `ctx.factr`, `ctx.pgtol`, and `ctx.print` can be set with `lbfgsb_config`.

     The possible values for `ctx.task` are:

     - `LBFGSB_START`: The algorithm has not yet started or has been restarted
       (see `lbfgsb_reset`).  The user can set the constraints of the problem
       and the parameters of the algorithm with `lbfgsb_config`.  The user
       should then call `lbfgsb_iterate` with the initial variables to
       effectively start a new minimization by the L-BFGS-B algorithm.  For
       this initial call, the bounds will be checked and the variables will be
       made feasible, there is no needs to compute the objective function value
       and its gradient.

     - `LBFGSB_FG`: The caller is requested to compute the value of the
       objective function and its gradient at the current variables `x`.

     - `LBFGSB_NEW_X`: The current variables `x` are available for inspection.
       This occurs for the initial variables (after they have been made
       feasible and the corresponding objective function and gradient computed)
       and after each iteration of the algorithm (that is at the end of each
       line-search).

     - `LBFGSB_CONVERGENCE`: The algorithm has converged, `x`, `f` and `g`
       contain the solution, the corresponding objective function and its
       gradient.

     - `LBFGSB_STOP`: The algorithm has been stopped by the caller.

     - `LBFGSB_WARNING`: Algorithm terminated on a non-fatal error.

     - `LBFGSB_ERROR`: Abnormal termination or ba algorithm settings.

   SEE ALSO: lbfgsb_config, lbfgsb_iterate, lbfgsb_reset, lbfgsb_stop.
 */

extern lbfgsb_config;
/* DOCUMENT lbfgsb_config, ctx, factr=, pgtol=, print=, upper=, lower=;

     Configure parameters of L-BFGS-B algorithm in context `ctx`.  The current
     tast `ctx.task` must be `LBFGSB_START`, that is `lbfgsb_iterate` must not
     have been called since the creation of the context or the last call to
     `lbfgsb_reset`.

     Keywords `lower` and `upper` are to specify the bounds on the variables.
     Bounds may be specified as a scalar value (for uniform bounds) or as an
     array of dimensions `ctx.dims` (for element-wise bounds).

     Keyword `factr ≥ 0` is to specify a relative function decrease threshold.
     The iteration will stop when:

         (f_{k} - f_{k+1})/max{|f_{k}|,|f_{k+1}|,1} <= factr⋅epsmch

     where `f_{k}` denotes the objective function value at `k`-th iteration,
     `epsmch` is the machine precision, which is automatically generated by the
     code.  Typical values for `factr`: 1e+12 for low accuracy; 1e+7 for
     moderate accuracy; 1e+1 for extremely high accuracy.

     Keyword `pgtol ≥ 0` is to specify a threshold on the infinite norm of the
     projected gradient.  The iteration will stop when:

         max{|pg_i | i = 1, ..., n} ≤ pgtol

     where `pg_i` is the `i`-th component of the projected gradient.

     If called as a function, returns the context `ctx`.

   SEE ALSO: lbfgsb_create, lbfgsb_iterate, lbfgsb_reset, lbfgsb_stop.
 */

extern lbfgsb_reset;
/* DOCUMENT lbfgsb_reset, ctx, bounds=0/1;

     Restart L-BFGS-B algorithm for context `ctx`.  This is useful to re-use
     the same context (avoiding re-allocation) for the same variables
     dimensions and maximum number of memorized steps.  The problem bounds and
     oher algorithm settings can be configured with `lbfgsb_config`.

     If keyword `bounds` is set true, the lower and upper bounds are reset to
     `-Inf` and `+Inf` (where `Inf` donotes the numerical infinity as given by
     `LBFGSB_INFINITY`).

     If called as a function, returns the context `ctx`.

   SEE ALSO: lbfgsb_create, lbfgsb_config, lbfgsb_iterate, lbfgsb_stop.
 */

extern lbfgsb_stop;
/* DOCUMENT task = lbfgsb_stop(ctx, reason);

     Stop L-BFGS-B algorithm for context `ctx` with a string indicating the
     resaon of the termination.  Argument `reason` is used to set `ctx.reason`
     and `ctx.task` (which is returned).  Argument `reason` must starts with
     "STOP:", "WARNING:", or "ERROR:" to have `ctx.task` respectively set to
     `LBFGSB_STOP`, `LBFGSB_WARNING`, or `LBFGSB_ERROR`.

   SEE ALSO: lbfgsb_create, lbfgsb_config, lbfgsb_iterate, lbfgsb_stop.
 */

local LBFGSB_START, LBFGSB_FG, LBFGSB_NEW_X, LBFGSB_CONVERGENCE;
local LBFGSB_STOP, LBFGSB_WARNING, LBFGSB_ERROR;
extern lbfgsb_iterate;
/* DOCUMENT task = lbfgsb_iterate(ctx, x, f, g);

     Execute next step of L-BFGS-B algorithm for context `ctx` for the
     variables `x`, objective function value `f` and its gradient `g`.
     Arguments `x`, `f` and `g` are used for input and output, they must not be
     expressions.

     For the first algorithm step or after a restart (when `ctx.task` is
     `LBFGSB_START`), the bounds (`ctx.lower` and `ctx.upper`) will be checked
     and the variables `x` will be made feasible, there is no needs to compute
     the objective function value and its gradient (variables `f` and `g` must
     however be specified, their value is irrelevant).

   SEE ALSO: lbfgsb_create, lbfgsb_config, lbfgsb_reset, lbfgsb_stop.
 */

extern lbfgsb_pgnorm2;
/* DOCUMENT norm = lbfgsb_pgnorm2(ctx, x, g);

     Compte the Euclidean norm of the projected gradient for the bounds stored
     by L-BFGS-B context `ctx` given `g` the gradient of the objective function
     at variables `x`.  The variables must be feasible.

   SEE ALSO: lbfgsb_iterate.
 */

local LBFGSB_INFINITY, LBFGSB_NAN;
extern lbfgsb_init;
/* DOCUMENT lbfgsb_init;

     (Re)define constants of L-BFGS-B plugin.

     Among these constants are the possible task values `LBFGSB_START`,
     `LBFGSB_FG`, `LBFGSB_NEW_X`, `LBFGSB_CONVERGENCE`, `LBFGSB_STOP`,
     `LBFGSB_WARNING`, and `LBFGSB_ERROR`.

     Other constants:

     - `LBFGSB_INFINITY`: the floating-point representation of positive
       infinity;

     - `LBFGSB_NAN`: the floating-point representation of a quiet NaN (Not a
       Number).

   SEE ALSO: lbfgsb_create, lbfgsb_config, lbfgsb_iterate, lbfgsb_reset,
             lbfgsb_stop.
 */
lbfgsb_init;
