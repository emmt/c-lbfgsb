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

     - `ctx.lower`: the lower bounds, initially all set to `-Inf`;

     - `ctx.upper`: the upper bounds, initially all set to `+Inf`;

     - `ctx.print`: the level of verbosity (negative to suppress all printed
       information);

     - `ctx.task`: the current stage of the algorithm;

     - `ctx.reason`: a detailed textula infromation of the current stage of the
       algorithm.

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