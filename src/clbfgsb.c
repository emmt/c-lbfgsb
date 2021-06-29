#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <math.h>
#include "lbfgsb.h"

// Allocate a dynamic array of `n` elements of type `T`.
#define NEW_ARRAY(n, T)  ((T*)malloc((n)*sizeof(T)))

// Allocate a dynamic array of `n` elements of type `T` filled with zeros.
#define ZEROS(n, T)  ((T*)malloc((n)*sizeof(T)))

double lbfgsb_timer(void)
{
    double t;
    LBFGSB_TIMER_(&t);
    return t;
}

static inline lbfgsb_task get_task(
    const character* buf)
{
    switch (buf[0]) {
    case 'F':
        if (buf[1] == 'G') {
            return LBFGSB_FG;
        }
        break;
    case 'N':
        if (buf[1] == 'E' && buf[2] == 'W' &&
            buf[3] == '_' && buf[4] == 'X') {
            return LBFGSB_NEW_X;
        }
        break;
    case 'C':
        if (buf[1] == 'O' && buf[2] == 'N' && buf[3] == 'V') {
            return LBFGSB_CONVERGENCE;
        }
        break;
    case 'S':
        if (buf[1] == 'T' && buf[2] == 'A' &&
            buf[3] == 'R' && buf[4] == 'T') {
            return LBFGSB_START;
        }
        if (buf[1] == 'T' && buf[2] == 'O' && buf[3] == 'P') {
            return LBFGSB_STOP;
        }
        break;
    case 'W':
        if (buf[1] == 'A' && buf[2] == 'R' && buf[3] == 'N') {
            return LBFGSB_WARNING;
        }
    }
    return LBFGSB_ERROR;
}

char* lbfgsb_get_task_string(
    lbfgsb_context* ctx,
    char*           buf,
    long            siz)
{
    if (buf == NULL || siz < 1) {
        return "";
    }
    long len = siz - 1;
    if (len > LBFGSB_TASK_LENGTH) {
        len = LBFGSB_TASK_LENGTH;
    }
    const character* task = ctx->wrks.task;
    long ilast = -1;
    for (long i = 0; i < len; ++i) {
        int c = task[i];
        if (c == 0) {
            break;
        }
        if (c != ' ') {
            ilast = i;
        }
        buf[i] = c;
    }
    buf[ilast+1] = '\0';
    return buf;
}

lbfgsb_task lbfgsb_set_task(
    lbfgsb_context* ctx,
    const char*     str)
{
    long len1 = (str == NULL ? 0 : strlen(str));
    long len2 = LBFGSB_TASK_LENGTH;
    if (len1 > len2) {
        len1 = len2;
    }
    character* task = ctx->wrks.task;
    for (long i = 0; i < len1; ++i) {
        task[i] = str[i];
    }
    for (long i = len1; i < len2; ++i) {
        task[i] = ' ';
    }
    ctx->task = get_task(ctx->wrks.task);
    return ctx->task;
}

void lbfgsb_reset(
    lbfgsb_context* ctx,
    int full)
{
    if (full != 0) {
        long         n = ctx->siz;
        double*  lower = ctx->lower;
        double*  upper = ctx->upper;
        integer* bound = ctx->wrks.nbd;
        for (long i = 0; i < n; ++i) {
            lower[i] = -INFINITY;
            upper[i] = +INFINITY;
            bound[i] = 0;
        }
    }
    lbfgsb_set_task(ctx, "START");
}

lbfgsb_context* lbfgsb_create(
    long n,
    long m)
{
    // FIXME: check for integer overflow (in particular in n_wa below)
    if (n < 1 || m < 1) {
        errno = EINVAL;
        return NULL;
    }
    lbfgsb_context* ctx = malloc(sizeof(lbfgsb_context));
    if (ctx == NULL) {
        return NULL;
    }
    memset(ctx, 0, sizeof(lbfgsb_context));
    ctx->mem = m;
    ctx->siz = n;
    ctx->factr = 1.0e+7;
    ctx->pgtol = 1.0e-6;
    ctx->print = -1; // No output.
#ifdef OLD_LBFGSB_VERSION
   long n_wa = (2*m + 4)*n + 12*m*(m + 1);
#else
    long n_wa = (2*m + 5)*n + (11*m + 8)*m;
#endif
    if ((ctx->lower    = ZEROS(n,    double))  == NULL ||
        (ctx->upper    = ZEROS(n,    double))  == NULL ||
        (ctx->wrks.nbd = ZEROS(n,    integer)) == NULL ||
        (ctx->wrks.wa  = ZEROS(n_wa, double))  == NULL ||
        (ctx->wrks.iwa = ZEROS(3*n,  integer)) == NULL) {
        lbfgsb_destroy(ctx);
        return NULL;
    }
    lbfgsb_reset(ctx, 1);
    return ctx;
}


#define free_memory(ptr)                        \
    do {                                        \
        void *ptr_ = (ptr);                     \
        if (ptr_ != NULL) {                     \
            (ptr) = NULL;                       \
            free(ptr_);                         \
        }                                       \
    } while (0)

void lbfgsb_destroy(lbfgsb_context* ctx)
{
    if (ctx != NULL) {
        free_memory(ctx->lower);
        free_memory(ctx->upper);
        free_memory(ctx->wrks.nbd);
        free_memory(ctx->wrks.wa);
        free_memory(ctx->wrks.iwa);
        free(ctx);
    }
}

static void check_bounds(
    lbfgsb_context* ctx,
    double          x[])
{
    long           n     = ctx->siz;
    const double*  lower = ctx->lower;
    const double*  upper = ctx->upper;
    integer*       bound = ctx->wrks.nbd;
    for (long i = 0; i < n; ++i) {
        double lo = lower[i];
        double hi = upper[i];
        if (isnan(lo)) {
            lbfgsb_set_task(ctx, "ERROR: Invalid lower bound value");
            break;
        }
        if (isnan(hi)) {
            lbfgsb_set_task(ctx, "ERROR: Invalid upper bound value");
            break;
        }
        if (lo > hi) {
            lbfgsb_set_task(ctx, "ERROR: Incompatible bounds");
            break;
        }
        if (lo > -INFINITY) {
            if (hi < +INFINITY) {
                bound[i] = 2;
            } else {
                bound[i] = 1;
            }
        } else {
            if (hi < +INFINITY) {
                bound[i] = 3;
            } else {
                bound[i] = 0;
            }
        }
#if 0 // Not needed, this is done by main L-BFGS-B subroutine.
        if (x[i] > hi) {
            x[i] = hi;
        }
        if (x[i] < lo) {
            x[i] = lo;
        }
#endif
    }
}

lbfgsb_task lbfgsb_iterate(
    lbfgsb_context* ctx,
    double          x[],
    double*         f,
    double          g[])
{
    if (ctx->task == LBFGSB_START) {
        check_bounds(ctx, x);
    }
    if (ctx->task != LBFGSB_ERROR) {
        integer m = ctx->mem;
        integer n = ctx->siz;
        integer print = ctx->print;
        LBFGSB_SETULB_(
            &n, &m, x, ctx->lower, ctx->upper, ctx->wrks.nbd, f, g,
            &ctx->factr, &ctx->pgtol, ctx->wrks.wa, ctx->wrks.iwa,
            ctx->wrks.task, &print, ctx->wrks.csave, ctx->wrks.lsave,
            ctx->wrks.isave, ctx->wrks.dsave);
        ctx->task = get_task(ctx->wrks.task);
    }
    return ctx->task;
}

lbfgsb_task lbfgsb_get_task(
    const lbfgsb_context* ctx)
{
    return ctx->task;
}

long lbfgsb_get_mem(
    const lbfgsb_context* ctx)
{
    return ctx->mem;
}

long lbfgsb_get_siz(
    const lbfgsb_context* ctx)
{
    return ctx->siz;
}

double *lbfgsb_get_lower(
    const lbfgsb_context* ctx)
{
    return ctx->lower;
}

double *lbfgsb_get_upper(
    const lbfgsb_context* ctx)
{
    return ctx->upper;
}

double lbfgsb_get_factr(
    const lbfgsb_context* ctx)
{
    return ctx->factr;
}

void lbfgsb_set_factr(
    lbfgsb_context* ctx,
    double factr)
{
    ctx->factr = factr;
}

double lbfgsb_get_pgtol(
    const lbfgsb_context* ctx)
{
    return ctx->pgtol;
}

void lbfgsb_set_pgtol(
    lbfgsb_context* ctx,
    double pgtol)
{
    ctx->pgtol = pgtol;
}

long lbfgsb_get_print(
    const lbfgsb_context* ctx)
{
    return ctx->print;
}

void lbfgsb_set_print(
    lbfgsb_context* ctx,
    long print)
{
    ctx->print = print;
}

const double* lbfgsb_get_latest_x(
    const lbfgsb_context* ctx)
{
    long n = ctx->siz;
    long m = ctx->mem;

    // The latest iterate is contained in `wa[j..j+n-1]` where:
    long j = 3*n + 2*m*n + 11*m*m;
    return ctx->wrks.wa + j;
}
