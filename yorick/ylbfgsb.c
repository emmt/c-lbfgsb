// yor-tao-rt.c --
//
// Implements Yorick interface to TAO real-time software.
//
//-----------------------------------------------------------------------------
//
// Copyright (C) 2018-2021, Éric Thiébaut <eric.thiebaut@univ-lyon1.fr>
//
// See LICENSE.md for details.

#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#include <lbfgsb.h>

#include <pstdlib.h>
#include <play.h>
#include <yapi.h>
#include <ydata.h>

// We assume that Yorick char type is unsigned.
typedef unsigned char byte;
#define Y_BYTE Y_CHAR

// Define some macros to get rid of some GNU extensions when not compiling
// with GCC.
#if ! (defined(__GNUC__) && __GNUC__ > 1)
#   define __attribute__(x)
#   define __inline__
#   define __FUNCTION__        ""
#   define __PRETTY_FUNCTION__ ""
#endif

PLUG_API void y_error(const char *) __attribute__ ((noreturn));

static void push_string(const char* str)
{
    ypush_q(NULL)[0] = p_strcpy(str);
}

static inline long numberof(const long dims[])
{
    long n = 1;
    if (dims != NULL) {
        long ndims = dims[0];
        for (long i = 1; i <= ndims; ++i) {
            n *= dims[i];
        }
    }
    return n;
}

static void grow_dims(
    long dims[], int iarg, long max_ndims)
{
    long ndims = dims[0];
    int type = yarg_typeid(iarg);
    if (max_ndims == -1) {
        max_ndims = Y_DIMSIZE - 1;
    }
    if (type <= Y_LONG) {
        int rank = yarg_rank(iarg);
        if (rank == 0) {
            long dim = ygets_l(iarg);
            if (dim < 1) {
                goto bad_dim;
            }
            if (ndims + 1 > max_ndims) {
                goto too_many_dims;
            }
            ++ndims;
            dims[0] = ndims;
            dims[ndims] = dim;
            return;
        } else if (rank == 1) {
            long i, ntot;
            long* vals = ygeta_l(iarg, &ntot, NULL);
            if (ntot < 1 || vals[0] != ntot - 1) {
                goto bad_dimlist;
            }
            if (ndims + ntot - 1 > max_ndims) {
                goto too_many_dims;
            }
            for (i = 1; i < ntot; ++i) {
                long dim = vals[i];
                if (dim < 1) {
                    goto bad_dim;
                }
                ++ndims;
                dims[0] = ndims;
                dims[ndims] = dim;
            }
            return;
        }
    } else if (type == Y_VOID) {
        return;
    }
 bad_dimlist:
    y_error("invalid dimension list");
 bad_dim:
    y_error("invalid dimension");
 too_many_dims:
    y_error("too many dimensions");
}

static inline int same_dims(
    const long* a,
    const long* b)
{
    long z[1] = {0};
    if (a == NULL) {
        a = z;
    }
    if (b == NULL) {
        b = z;
    }
    if (a == b) {
        return 1;
    }
    long ndims;
    if ((ndims = a[0]) != b[0]) {
        return 0;
    }
    for (long i = 1; i <= ndims; ++i) {
        if (a[i] != b[i]) {
            return 0;
        }
    }
    return 1;
}

static double* push_copy_d(
    const double src[],
    long dims[])
{
    double* dst = ypush_d(dims);
    long n = numberof(dims);
    for (long i = 0; i < n; ++i) {
        dst[i] = src[i];
    }
    return dst;
}

static const char* task_name(
    lbfgsb_task task)
{
    switch (task) {
    case LBFGSB_START:       return "START";
    case LBFGSB_FG:          return "FG";
    case LBFGSB_NEW_X:       return "NEW_X";
    case LBFGSB_CONVERGENCE: return "CONVERGENCE";
    case LBFGSB_STOP:        return "STOP";
    case LBFGSB_WARNING:     return "WARNING";
    case LBFGSB_ERROR:       return "ERROR";
    default:                 return "UNKNOWN";
    }
}

//-----------------------------------------------------------------------------
// OBJECTS

typedef struct context {
    lbfgsb_context* ctx;
    long dims[Y_DIMSIZE];
} context;

static void free_context(void* addr)
{
    context* obj = (context*)addr;
    lbfgsb_destroy(obj->ctx);
}

static void print_context(void* addr)
{
    char buf[LBFGSB_TASK_LENGTH+1];
    context* obj = (context*)addr;
    lbfgsb_context* ctx = obj->ctx;
    y_print("L-BFGS-B context (siz=", 0);
    sprintf(buf, "%ld", ctx->siz);
    y_print(buf, 0);
    y_print(", mem=", 0);
    sprintf(buf, "%ld", ctx->mem);
    y_print(buf, 0);
    y_print(", dims=[", 0);
    for (long i = 0; i <= obj->dims[0]; ++i) {
        buf[0] = ',';
        sprintf((i == 0 ? buf : buf+1), "%ld", obj->dims[i]);
        y_print(buf, 0);
    }
    y_print("], print=", 0);
    sprintf(buf, "%ld", (long)ctx->print);
    y_print(buf, 0);
    y_print(", factr=", 0);
    sprintf(buf, "%.2g", (double)ctx->factr);
    y_print(buf, 0);
    y_print(", pgtol=", 0);
    sprintf(buf, "%.2g", (double)ctx->pgtol);
    y_print(buf, 0);
    y_print(", task=", 0);
    y_print(task_name(ctx->task), 0);
    y_print(")", 1);
}

// FIXME: ctx(x, f, g) -> task (like lbfgsb_iterate).
static void eval_context(void* addr, int argc)
{
    ypush_nil();
}

static void extract_context(void* addr, char* name)
{
    char buf[LBFGSB_TASK_LENGTH+1];
    context* obj = (context*)addr;
    lbfgsb_context* ctx = obj->ctx;
    switch (name[0]) {
    case 'd':
        if (strcmp(name, "dims") == 0) {
            int ndims = obj->dims[0];
            long* dims = ypush_l((long[2]){1, ndims+1});
            for (int d = 0; d <= ndims; ++d) {
                dims[d] =  obj->dims[d];
            }
            return;
        }
        break;
    case 'e':
        if (strcmp(name, "epsmch") == 0) {
            ypush_double(LBFGSB_EPSMCH(ctx));
            return;
        }
        break;
    case 'f':
        if (strcmp(name, "factr") == 0) {
            ypush_double(ctx->factr);
            return;
        }
        break;
   case 'l':
        if (strcmp(name, "lower") == 0) {
            push_copy_d(ctx->lower, obj->dims);
            return;
        }
        break;
    case 'm':
        if (strcmp(name, "mem") == 0) {
            ypush_long(ctx->mem);
            return;
        }
        break;
    case 'n':
        if (strcmp(name, "niters") == 0) {
            ypush_long(LBFGSB_NUM_ITER(ctx));
            return;
        }
        if (strcmp(name, "nevals") == 0) {
            ypush_long(LBFGSB_NTOT_FG(ctx));
            return;
        }
        if (strcmp(name, "nskips") == 0) {
            ypush_long(LBFGSB_NTOT_SKIP(ctx));
            return;
        }
        break;
    case 'p':
        if (strcmp(name, "pgnorm") == 0) {
            ypush_double(LBFGSB_PG_NORMINF(ctx));
            return;
        }
        if (strcmp(name, "pgtol") == 0) {
            ypush_double(ctx->pgtol);
            return;
        }
        if (strcmp(name, "print") == 0) {
            ypush_long(ctx->print);
            return;
        }
        break;
    case 'r':
        if (strcmp(name, "reason") == 0) {
            push_string(lbfgsb_get_task_string(ctx, buf, sizeof(buf)));
            return;
        }
        break;
    case 's':
        if (strcmp(name, "siz") == 0) {
            ypush_long(ctx->siz);
            return;
        }
        if (strcmp(name, "step") == 0) {
            ypush_double(LBFGSB_STEP(ctx));
            return;
        }
        break;
    case 't':
        if (strcmp(name, "task") == 0) {
            ypush_long(ctx->task);
            return;
        }
        if (strcmp(name, "theta") == 0) {
            ypush_double(LBFGSB_THETA(ctx));
            return;
        }
        break;
    case 'u':
        if (strcmp(name, "upper") == 0) {
            push_copy_d(ctx->upper, obj->dims);
            return;
        }
        break;
    }
    y_error("bad member");
}

static y_userobj_t context_type = {
    "L-BFGS-B context",
    free_context,
    print_context,
    eval_context,
    extract_context,
    NULL
};

static context* get_context(
    int iarg)
{
    return (context*)yget_obj(iarg, &context_type);
}

//-----------------------------------------------------------------------------
// BUILT-IN FUNCTIONS

void Y_lbfgsb_create(
    int argc)
{
    if (argc != 2) {
        y_error("usage: lbfgsb_create(dims, mem)");
    }
    long dims[Y_DIMSIZE] = {0};
    grow_dims(dims, argc - 1, -1);
    long mem =  ygets_l(argc - 2);
    if (mem < 1) {
        y_error("argument `mem` must be at least 1");
    }
    long siz = numberof(dims);
    context* obj = ypush_obj(&context_type, sizeof(context));
    for (long i = 0; i <= dims[0]; ++i) {
        obj->dims[i] = dims[i];
    }
    obj->ctx = lbfgsb_create(siz, mem);
    if (obj->ctx == NULL) {
        int code = errno;
        y_error(
            (code == ENOMEM ?
             "insufficient memory for creating L-BFGS-B context" :
             (code == EINVAL ?
              "invalid argument for creating L-BFGS-B context" :
              "failed to create L-BFGS-B context")));
    }
}

static void set_bound(context* obj, double* dst, int iarg)
{
    int rank = yarg_rank(iarg);
    if (rank == 0) {
        // Scalar.
        double val = ygets_d(iarg);
        long n = obj->ctx->siz;
        for (long i = 0; i < n; ++i) {
            dst[i] =val;
        }
    } else if (rank > 0) {
        long dims[Y_DIMSIZE], n;
        const double* src = ygeta_d(iarg, &n, dims);
        if (!same_dims(dims, obj->dims)) {
            y_error((dst == obj->ctx->lower ?
                     "value of `lower` has incompatible dimensions" :
                     "value of `upper` has incompatible dimensions"));
        }
        for (long i = 0; i < n; ++i) {
            dst[i] = src[i];
        }
    } else {
        y_error((dst == obj->ctx->lower ?
                 "non-array value for `lower`" :
                 "non-array value for `upper`"));
    }
}

void Y_lbfgsb_config(
    int argc)
{
    // Keyword unique indices.
    static long factr_index = -1L;
    static long lower_index = -1L;
    static long pgtol_index = -1L;
    static long print_index = -1L;
    static long upper_index = -1L;
    if (upper_index == -1L) {
        factr_index = yget_global("factr", 0);
        lower_index = yget_global("lower", 0);
        pgtol_index = yget_global("pgtol", 0);
        print_index = yget_global("print", 0);
        upper_index = yget_global("upper", 0);
    }

    // First pass on positional arguments.
    int drop = 0;
    context* obj = NULL;
    for (int iarg = argc - 1; iarg >= 0; --iarg) {
        long index = yarg_key(iarg);
        if (index < 0) {
            // Positional argument.
            if (obj == NULL) {
                obj = get_context(iarg);
                drop = iarg;
            } else {
                y_error("too many arguments");
            }
        } else {
            // Keyword argument.
            --iarg;
        }
    }
    if (obj == NULL) {
        y_error("missing L-BFGS-B context");
    }
    lbfgsb_context* ctx = obj->ctx;
    if (ctx->task != LBFGSB_START) {
        y_error("call `lbfgsb_reset` before `lbfgsb_config`");
    }

    // Second pass on keywords.
    for (int iarg = argc - 1; iarg >= 0; --iarg) {
        long index = yarg_key(iarg);
        if (index < 0) {
            // Positional argument.
            continue;
        } else {
            // Keyword argument.
            --iarg;
            if (index == factr_index) {
                if (!yarg_nil(iarg)) {
                    double factr = ygets_d(iarg);
                    if (isnan(factr) || factr < 0) {
                        y_error("bad value for parameter `factr`");
                    }
                    ctx->factr = factr;
                }
            } else if (index == lower_index) {
                if (!yarg_nil(iarg)) {
                    set_bound(obj, ctx->lower, iarg);
                }
            } else if (index == pgtol_index) {
                if (!yarg_nil(iarg)) {
                    double pgtol = ygets_d(iarg);
                    if (isnan(pgtol) || pgtol < 0) {
                        y_error("bad value for parameter `pgtol`");
                    }
                    ctx->pgtol = pgtol;
                }
            } else if (index == print_index) {
                if (!yarg_nil(iarg)) {
                    ctx->print = ygets_l(iarg);
                }
            } else if (index == upper_index) {
                if (!yarg_nil(iarg)) {
                    set_bound(obj, ctx->upper, iarg);
                }
            } else {
                y_error("unsupported keyword");
            }
        }
    }
    if (drop > 0) {
        yarg_drop(drop);
    }
}

void Y_lbfgsb_stop(
    int argc)
{
    if (argc != 2) {
        y_error("usage: lbfgsb_stop(ctx, reason)");
    }
    context* obj = get_context(argc - 1);
    const char* reason = ygets_q(argc - 2);
    int valid = 0;
    if (reason != NULL) {
        switch (reason[0]) {
        case 'C':
            valid = (strncmp(reason, "CONVERGENCE:", 12) == 0);
            break;
        case 'S':
            valid = (strncmp(reason, "STOP:", 5) == 0);
            break;
        case 'W':
            valid = (strncmp(reason, "WARNING:", 8) == 0);
            break;
        case 'E':
            valid = (strncmp(reason, "ERROR:", 6) == 0);
            break;
        }
    }
    if (!valid) {
        y_error("argument `reason` must start with `\"CONVERGENCE:\", "
                "`\"STOP:\", `\"WARNING:\", or  `\"ERROR:\"");
    }
    ypush_long(lbfgsb_set_task(obj->ctx, reason));
}

void Y_lbfgsb_reset(
    int argc)
{
    static long bounds_index = -1L;
    if (bounds_index == -1L) {
        bounds_index = yget_global("bounds", 0);
    }
    int flags = 0;
    int drop = 0;
    context* obj = NULL;
    for (int iarg = argc - 1; iarg >= 0; --iarg) {
        long index = yarg_key(iarg);
        if (index < 0) {
            // Positional argument.
            if (obj == NULL) {
                obj = get_context(iarg);
                drop = iarg;
            } else {
                y_error("too many arguments");
            }
        } else {
            // Keyword argument.
            --iarg;
            if (index == bounds_index) {
                if (yarg_true(iarg)) {
                    flags |= 1;
                }
            } else {
                y_error("unsupported keyword");
            }
        }
    }
    if (obj == NULL) {
        y_error("missing L-BFGS-B context");
    }
    lbfgsb_reset(obj->ctx, flags);
    if (drop > 0) {
        yarg_drop(drop);
    }
}

void Y_lbfgsb_iterate(
    int argc)
{
    if (argc != 4) {
        y_error("usage: lbfgsb_iterate(ctx, x, f, g)");
    }

    // Get context (1st argument).
    context* obj = get_context(argc - 1);

    // Get x (2nd argument).
    int x_iarg = argc - 2;
    long x_index = yget_ref(x_iarg);
    if (x_index < 0) {
        y_error("variables `x` must not be a temporary expression");
    }
    long x_dims[Y_DIMSIZE], x_ntot;
    int x_type = Y_VOID;
    double* x = ygeta_any(x_iarg, &x_ntot, x_dims, &x_type);
    if (!same_dims(x_dims, obj->dims)) {
        y_error("variables `x` have incompatible dimensions");
    }
    if (x_type < Y_CHAR || x_type > Y_DOUBLE) {
        y_error("variables `x` have non-real type");
    }

    // Get f (3rd argument).
    int f_iarg = argc - 3;
    long f_index = yget_ref(f_iarg);
    if (f_index < 0) {
        y_error("function value `f` must not be a temporary expression");
    }
    double f;
    int f_type = yarg_typeid(f_iarg);
    if (f_type >= Y_CHAR && f_type <= Y_DOUBLE && yarg_rank(f_iarg) == 0) {
        f = ygets_d(f_iarg);
    } else {
        f = 0.0; // FIXME: could be NaN
        if (obj->ctx->task != LBFGSB_START) {
            y_error("function value is undefined");
        }
        if (f_type != Y_VOID) {
            y_error("function value `f` must be initialized or a number");
        }
    }
    if (isnan(f)) {
        if (obj->ctx->task != LBFGSB_START) {
            y_error("function value is NaN");
        }
        f = 0.0; // FIXME: could be NaN
    }

    // Get g (4th argument).
    int g_iarg = argc - 4;
    long g_index = yget_ref(g_iarg);
    if (g_index < 0) {
        y_error("gradient `g` must not be a temporary expression");
    }
    long g_dims[Y_DIMSIZE], g_ntot;
    int g_type = Y_VOID;
    double* g = ygeta_any(g_iarg, &g_ntot, g_dims, &g_type);
    if (!same_dims(g_dims, obj->dims)) {
        y_error("gradient `g` has incompatible dimensions");
    }
    if (g_type < Y_CHAR || g_type > Y_DOUBLE) {
        y_error("gradient `g` has non-real type");
    }

    // All arguments have been checked. Convert inputs if needed and redefine
    // caller's variables.
    if (x_type != Y_DOUBLE) {
        x = ygeta_coerce(x_iarg, x, x_ntot, x_dims, x_type, Y_DOUBLE);
        yput_global(x_index, x_iarg);
    }
    if (g_type != Y_DOUBLE) {
        g = ygeta_coerce(g_iarg, g, g_ntot, g_dims, g_type, Y_DOUBLE);
        yput_global(g_index, g_iarg);
    }

    // Call L-BFGS-B iterator.
    long task = lbfgsb_iterate(obj->ctx, x, &f, g);

    // Redefine output `f` as needed.
    ypush_double(f);
    yput_global(f_index, 0);

    // Push result.
    ypush_long(task);
}

static void define_long(const char* name, long value)
{
    ypush_long(value);
    yput_global(yget_global(name, 0), 0);
    yarg_drop(1);
}

static void define_double(const char* name, double value)
{
    ypush_double(value);
    yput_global(yget_global(name, 0), 0);
    yarg_drop(1);
}

void Y_lbfgsb_init(
    int argc)
{
#define DEFINE_LONG(id) define_long(#id, id)
    DEFINE_LONG(LBFGSB_START);
    DEFINE_LONG(LBFGSB_FG);
    DEFINE_LONG(LBFGSB_NEW_X);
    DEFINE_LONG(LBFGSB_CONVERGENCE);
    DEFINE_LONG(LBFGSB_STOP);
    DEFINE_LONG(LBFGSB_WARNING);
    DEFINE_LONG(LBFGSB_ERROR);
#undef DEFINE_LONG
    define_double("LBFGSB_INFINITY", (double)INFINITY);
    define_double("LBFGSB_NAN", (double)NAN);
    ypush_nil();
}
