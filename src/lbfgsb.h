#ifndef LBFGSB_H
#define LBFGSB_H 1

#include <math.h>

/*
 * Since C99, the macro `NAN`, defined in `<math.h>`, expands a to constant
 * expression of type float which evaluates to a quiet not-a-number (QNaN)
 * value.  If the implementation does not support QNaNs, this macro
 * constant is not defined.
 *
 * See <https://en.cppreference.com/w/c>.
 */
#define LBFGSB_NAN NAN

/*
 * Since C99, the macro `INFINITY`, defined in `<math.h>`, expands to a
 * constant expression of type float which evaluates to positive or
 * unsigned infinity.
 *
 * Since C99, the macros `HUGE_VALF`, `HUGE_VAL` and `HUGE_VALL`, defined
 * in `<math.h>`, expand to positive floating point constant expressions
 * which compare equal to the values returned by floating-point functions
 * and operators in case of overflow.
 *
 * See <https://en.cppreference.com/w/c>.
 */
#define LBFGSB_INF INFINITY

#ifdef __cplusplus
extern "C" {
#endif

// The following definitions must match FORTRAN compiler settings.
typedef int  logical;
typedef int  integer;
typedef char character;

// Link name of a few FORTRAN subroutines in L-BFGS-B code.
#define LBFGSB_SETULB_ setulb_
#define LBFGSB_TIMER_  timer_

extern void LBFGSB_TIMER_(
    double* t);

extern void LBFGSB_SETULB_(
    const integer* n,
    const integer* m,
    double         x[],
    const double   l[],
    const double   u[],
    const integer  nbd[],
    double*        f,
    double         g[],
    const double*  factr,
    const double*  pgtol,
    double         wa[],
    integer        iwa[],
    character      task[],
    integer*       iprint,
    character      csave[],
    logical        lsave[],
    integer        isave[],
    double         dsave[]);

/**
 * L-BFGS-B task codes
 *
 * @see lbfgsb_get_task().
 */
typedef enum {
    LBFGSB_START       =  0,
    LBFGSB_FG          =  1,
    LBFGSB_NEW_X       =  2,
    LBFGSB_CONVERGENCE =  3,
    LBFGSB_STOP        =  4,
    LBFGSB_WARNING     =  5,
    LBFGSB_ERROR       =  6,
} lbfgsb_task;

#define LBFGSB_TASK_LENGTH 60

typedef struct lbfgsb_context {
    long        siz;   ///> Size of the problem (number of variables).
    long        mem;   ///> Maximum number of memorized steps.
    double*     lower; ///> Array of lower bounds.
    double*     upper; ///> Array of upper bounds.
    double      factr; ///> Tolerance factor for convergence in function value.
    double      pgtol; ///> Tolerance for convergence in projected gradient.
    int         task;  ///> Task to execute.
    int         print; ///> Verbosity setting.
    // Private workspaces.
    struct {
        integer*   nbd;
        double*    wa;
        integer*   iwa;
        character  task[LBFGSB_TASK_LENGTH];
        character  csave[LBFGSB_TASK_LENGTH];
        logical    lsave[4];
        integer    isave[44];
        double     dsave[29];
    } wrks;
} lbfgsb_context;

/**
 * @brief Create a new L-BFGS-B context.
 *
 * This method allocates memory for solving a box constrained optimization
 * problem by the L-BFGS-B algorithm.  The lower and upper bounds stored by the
 * context are initialized to `-Inf` and `+Inf` as if the problemn is
 * unconstrained.
 *
 * The task assigned to the newly created context is `LBFGSB_START`,
 * see lbfgsb_get_task() for a description of how to initiate an optimization
 * with this context.
 *
 * The same context can be used to solve several problems with the same numebr
 * of variables and maximum number of memorized steps, see lbfgsb_reset().
 *
 * It is the caller's responsibility to release allocated resources by calling
 * lbfgsb_destroy().
 *
 * @param siz     The number of variables of the problem.
 * @param mem     The maximum number of memorized steps.
 *
 * @return The address of the new context or `NULL` in case of failure.
 */
extern lbfgsb_context* lbfgsb_create(
    long siz,
    long mem);

/**
 * @brief Destroy L-BFGS-B context.
 *
 * Release resources associated with L-BFGS-B context created by
 * lbfgsb_create().
 *
 * @param ctx   The L-BFGS-B context.
 *
 * @return The next task for the given context.
 */
extern void lbfgsb_destroy(
    lbfgsb_context* ctx);

/**
 * @brief Restart L-BFGS-B algorithm.
 *
 * This function can be called to restart L-BFGS-B algorithm with the same
 * context.  This may be to solve another problem (with the same number of
 * variables and maximum number of memorized steps) with the same bounds or
 * not.  Calling lbfgsb_reset() is needed after an error.
 *
 * If `full` is non-zero, the lower and upper bounds stored by the context are
 * initialized to `-Inf` and `+Inf` as if the problemn is unconstrained.
 *
 * @param ctx   The L-BFGS-B context.
 * @param full  Also reset bounds if non-zero.
 *
 * @see lbfgsb_get_task(), lbfgsb_create(), lbfgsb_iterate().
 */
extern void lbfgsb_reset(
    lbfgsb_context* ctx,
    int full);

/**
 * @brief Iterate L-BFGS-B algorithm.
 *
 * @param ctx   The L-BFGS-B context.
 * @param x     The variables of the problem.
 * @param f     A pointer to the objective function value.
 * @param g     The gradient of the objective function.
 *
 * @return The next task for the given context.
 */
extern lbfgsb_task lbfgsb_iterate(
    lbfgsb_context* ctx,
    double          x[],
    double*         f,
    double          g[]);

/**
 * @brief Get current task.
 *
 * This function yields `ctx->task` and is provided for bindings in other
 * language than C/C++.  The possible values are:
 *
 * - `LBFGSB_START`: The algorithm has not yet started or has been restarted
 *   (see lbfgsb_reset()).  The user can set the constraints of the problem
 *   (that is the bounds on the variables in `ctx->lower` and `ctx->upper`) and
 *   the parameters of the algorithm (`ctx->factr`, `ctx->pgtol` and
 *   `ctx->print`).  The user should then call lbfgsb_iterate() with the
 *   initial variables to effectively start a new minimization by the L-BFGS-B
 *   algorithm.  For this initial call, the bounds will be checked and the
 *   variables will be made feasible, there is no needs to compute the
 *   objective function value and its gradient.
 *
 * - `LBFGSB_FG`: The caller is requested to compute the value of the objective
 *   function and its gradient at the current variables `x`.
 *
 * - `LBFGSB_NEW_X`: The current variables `x` are available for inspection.
 *   This occurs for the initial variables (after they have been made feasible
 *   and the corresponding objective function and gradient computed) and after
 *   each iteration of the algorithm (that is at the end of each line-search).
 *
 * - `LBFGSB_CONVERGENCE`: The algorithm has converged, `x`, `*f` and `g`
 *   contain the solution, the corresponding objective function and its
 *   gradient.
 *
 * - `LBFGSB_STOP`: The algorithm has been stopped by the caller.
 *
 * - `LBFGSB_WARNING`: Algorithm terminated on a non-fatal error.
 *
 * - `LBFGSB_ERROR`: Abnormal termination or ba algorithm settings.
 *
 * @param ctx   The L-BFGS-B context.
 *
 * @return The current task for the given context.
 */
extern lbfgsb_task lbfgsb_get_task(
    const lbfgsb_context* ctx);

extern lbfgsb_task lbfgsb_set_task(
    lbfgsb_context* ctx,
    const char*     str);

extern char* lbfgsb_get_task_string(
    lbfgsb_context* ctx,
    char*           buf,
    long            siz);

/**
 * @brief Get maximum number of memorized steps.
 *
 * This function yields `ctx->mem` and is provided for bindings in other
 * language than C/C++.
 *
 * @param ctx   The L-BFGS-B context.
 *
 * @return The maximum number of previous steps memorized by the given context.
 */
extern long lbfgsb_get_mem(
    const lbfgsb_context* ctx);

/**
 * @brief Get number of variables.
 *
 * This function yields `ctx->siz` and is provided for bindings in other
 * language than C/C++.
 *
 * @param ctx   The L-BFGS-B context.
 *
 * @return The number of variables in the given context.
 */
extern long lbfgsb_get_siz(
    const lbfgsb_context* ctx);

/**
 * @brief Get array storing lower bounds.
 *
 * This function yields `ctx->lower` and is provided for bindings in other
 * language than C/C++.
 *
 * @param ctx   The L-BFGS-B context.
 *
 * @return The array of lower bound values in the given context.
 */
extern double *lbfgsb_get_lower(
    const lbfgsb_context* ctx);

/**
 * @brief Get array storing upper bounds.
 *
 * This function yields `ctx->upper` and is provided for bindings in other
 * language than C/C++.
 *
 * @param ctx   The L-BFGS-B context.
 *
 * @return The array of upper bound values in the given context.
 */
extern double *lbfgsb_get_upper(
    const lbfgsb_context* ctx);

extern double lbfgsb_get_factr(
    const lbfgsb_context* ctx);

extern void lbfgsb_set_factr(
    lbfgsb_context* ctx,
    double factr);

extern double lbfgsb_get_pgtol(
    const lbfgsb_context* ctx);

extern void lbfgsb_set_pgtol(
    lbfgsb_context* ctx,
    double pgtol);

extern long lbfgsb_get_print(
    const lbfgsb_context* ctx);

extern void lbfgsb_set_print(
    lbfgsb_context* ctx,
    long print);

extern double lbfgsb_timer(
    void);

extern const double* lbfgsb_get_latest_x(
    const lbfgsb_context* ctx);

#ifdef __cplusplus
}
#endif

#define LBFGSB_DSAVE_(ctx, i) ((ctx)->wrks.dsave[i])
#define LBFGSB_ISAVE_(ctx, i) ((ctx)->wrks.isave[i])
#define LBFGSB_LSAVE_(ctx, i) ((ctx)->wrks.lsave[i])

// On exit with `task == LBFGSB_NEW_X`, the following information is available:

// - `dsave[0]` = current `theta` in the BFGS matrix;
#define LBFGSB_THETA(ctx) LBFGSB_DSAVE_(ctx,0)

// - `dsave[1] = f(x)` in the previous iteration;
#define LBFGSB_PREV_F(ctx) LBFGSB_DSAVE_(ctx,1)

// - `dsave[2] = factr*epsmch`;
#define LBFGSB_F_TEST(ctx) LBFGSB_DSAVE_(ctx,2)

// - `dsave[3]` = 2-norm of the line search direction vector;
#define LBFGSB_D_NORM2(ctx) LBFGSB_DSAVE_(ctx,3)

// - `dsave[4]` = the machine precision epsmch generated by the code;
#define LBFGSB_EPSMCH(ctx) LBFGSB_DSAVE_(ctx,4)

// - `dsave[6]` = the accumulated time spent on searching for Cauchy points;
#define LBFGSB_CAUCHY_TIME(ctx) LBFGSB_DSAVE_(ctx,6)

// - `dsave[7]` = the accumulated time spent on subspace minimization;
#define LBFGSB_SUBSPACE_TIME(ctx) LBFGSB_DSAVE_(ctx,7)

// - `dsave[8]` = the accumulated time spent on line search;
#define LBFGSB_LNSRCH_TIME(ctx) LBFGSB_DSAVE_(ctx,8)

// - `dsave[10]` = the slope of the line search function at the current
//   point of line search;
#define LBFGSB_DF(ctx) LBFGSB_DSAVE_(ctx,10)

// - `dsave[11]` = the maximum relative step length imposed in line search;
#define LBFGSB_MAX_STEP(ctx) LBFGSB_DSAVE_(ctx,11)

// - `dsave[12]` = the infinity norm of the projected gradient;
#define LBFGSB_PG_NORMINF(ctx) LBFGSB_DSAVE_(ctx,12)

// - `dsave[13]` = the relative step length in the line search;
#define LBFGSB_STEP(ctx) LBFGSB_DSAVE_(ctx,13)

// - `dsave[14]` = the slope of the line search function at the starting
//   point of the line search;
#define LBFGSB_DF0(ctx) LBFGSB_DSAVE_(ctx,14)

// - `dsave[15]` = the square of the 2-norm of the line search direction
//   vector.
#define LBFGSB_D_NORM2_SQUARED(ctx) LBFGSB_DSAVE_(ctx,15)

// - If `lsave[0] == true` then the initial X has been replaced by its
//   projection in the feasible set;
#define LBFGSB_INITIAL_X_UNFEASIBLE(ctx) (LBFGSB_LSAVE_(ctx,0) != 0)

// - If `lsave[1] == true` then the problem is constrained;
#define LBFGSB_CONSTRAINED(ctx) (LBFGSB_LSAVE_(ctx,1) != 0)

// - If `lsave[2] == true` then each variable has upper and lower bounds;
#define LBFGSB_FULLY_CONSTRAINED(ctx) (DLBFGSB_LSAVE_(ctx,2) != 0)

// - `isave[21]` is the total number of intervals explored in the search of
//   Cauchy points;
#define LBFGSB_NTOT_CAUCHY(ctx) LBFGSB_ISAVE_(ctx,21)

// - `isave[25]` is the total number of skipped BFGS updates before the
//   current iteration;
#define LBFGSB_NTOT_SKIP(ctx) LBFGSB_ISAVE_(ctx,25)

// - `isave[29]` is the number of current iteration;
#define LBFGSB_NUM_ITER(ctx) LBFGSB_ISAVE_(ctx,29)

// - `isave[30]` is the total number of BFGS updates prior the current
//   iteration;
#define LBFGSB_NTOT_UPDT(ctx) LBFGSB_ISAVE_(ctx,30)

// - `isave[32]` is the number of intervals explored in the search of Cauchy
//   point in the current iteration;
#define LBFGSB_NUM_CAUCHY(ctx) LBFGSB_ISAVE_(ctx,32)

// - `isave[33]` is the total number of function and gradient evaluations;
#define LBFGSB_NTOT_FG(ctx) LBFGSB_ISAVE_(ctx,33)

// - `isave[35]` is the number of function value or gradient evaluations in
//   the current iteration;
#define LBFGSB_NUM_FG(ctx) LBFGSB_ISAVE_(ctx,35)

// - if `isave[36] == 0` then the subspace argmin is within the box;
//   if `isave[36] == 1` then the subspace argmin is beyond the box;
#define LBFGSB_WITHIN_BOX(ctx) (LBFGSB_ISAVE_(ctx,36) == 0)

// - `isave[37]` is the number of free variables in the current iteration;
#define LBFGSB_NUM_FREE(ctx) LBFGSB_ISAVE_(ctx,37)

// - `isave[38]` is the number of active constraints in the current
//   iteration;
#define LBFGSB_NUM_ACTIVE(ctx) LBFGSB_ISAVE_(ctx,38)

// - `n + 1 - isave[39]` is the number of variables leaving the set of
//   active constraints in the current iteration;
#define LBFGSB_NUM_LEAVING(ctx) ((ctx)->siz + 1 - LBFGSB_ISAVE_(ctx,39))

// - `isave[40]` is the number of variables entering the set of active
//   constraints in the current iteration.
#define LBFGSB_NUM_ENTERING(ctx) LBFGSB_ISAVE_(ctx,40)

#endif // LBFGS_H
