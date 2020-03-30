/*=================================================================
 *
 * tvden.c    A Mathlab .MEX file
 *
 *      Perform compute-intensive part of tvdenoise_oct() function.
 *      This is a 3-dimension minimization. This code uses native
 *      Windows function calls, and is written to the Matlab R2018a
 *      C API.
 *
 * The calling syntax is:
 *
 *      tvden(A, iter, C, dt, ncpu);
 *
 * where:
 *      - A is output, must be a pre-allocated matrix of identical dim as
 *        C. The content will be overwritten.
 *      - iter is number of iterations, type double
 *      - C is f_lambda, a 3D matrix of type single
 *      - dt is a scalar of type single
 *      - ncpu is number of threads to use, 1..16
 *
 * Compile:
 *     mex -R2018a tvden.c
 *
 * Tested on Matlab R2018a, Windows 10, using mingw64 compiler v5.3.
 *
 * Copyright 2019 Nokia Inc.
 *
 * @author Mike Scheutzow <mjs973@gmail.com>
 * @version 1.0
 *
 *=================================================================*/

/* needed to enable COND prototypes */
#define _WIN32_WINNT 0x600
#define WINVER 0x600

#include <math.h>
#include <string.h>
#include <windows.h>
#include <process.h>
#include "matrix.h"
#include "mex.h"

typedef struct inpArg {
    /* inputs */
    size_t iter;
    size_t h, w;
    size_t nfrm;
    size_t elem2d;              /* h*w */
    size_t elem3d;              /* h*w*nfrm */
    size_t frm_per_thread;      /* nfrm/nthread */
    const mxSingle *f_lambda;   /* 3d inp */
    const mxSingle *dt;         /* scalar inp */
} InpArg;

/* current system time */
typedef struct mTime {
    int sec;
    int millisec;
} MTime;

#define MAXCPU 16

struct thrInp;
typedef void (*WORKER)(struct thrInp *);

/* parameters for a single worker thread */
typedef struct thrInp {
    size_t id;      /* 0 to (MAXCPU-1) */
    HANDLE handle;  /* thread handle, 0=closed */

    /* set fields, then send signal to cond_wake. threads can
     * also allowed to read st->globInp. */
    WORKER func;       /* function to execute */
    volatile size_t running;    /* 1 if thread is busy, 0 if idle */

    /* input parameters for worker thread */
    mxSingle *a;        /* out frm[0] */
    const mxSingle *b;  /* inp frm[0] */
    const mxSingle *c;  /* inp frm[0] */
    const mxSingle *d;  /* inp frm[0] */
    size_t first_frm;   /* first frame idx for this thread, 0-based */
    size_t last_frm;    /* frame idx to stop at (don't process it) */
    size_t dir;         /* spatdiff direction, 0 if unused */
    MTime time_finish;  /* timestamp of worker complete */
} ThrInp;

/* persistent state while loaded */
typedef struct tvState {
    /* number of usable thrInp[] threads */
    size_t nthread;         /* 0 => state is uninitialized */
    size_t shutdown_registered;
    volatile size_t quit_flag;  /* 1=> stop threads */
    size_t verbose;
    const InpArg *globInp;

    /* start signal, from main->worker */
    SRWLOCK mutex_wake;
    CONDITION_VARIABLE cond_wake;
    /* end signal, from worker->main */
    SRWLOCK mutex_done;
    CONDITION_VARIABLE cond_done;

    /* args for each worker thread */
    ThrInp thrInp[MAXCPU];
} TvState;

/* storage for persistent state */
static TvState my_state;

static TvState *st = &my_state;

static void exec_threads(void);
static void mt_shutdown(void);

#if 0
/* Flush the mexPrintf buffer to cmd win. from libmwservices.lib,
 * Windows-only. No longer compiles in R2018a. */
extern bool ioFlush(void);
#else
/**
 * Flush the mexPrintf buffer to cmd win. This function can only be called
 * from main thread */
static void ioFlush(void) {
    mexEvalString("drawnow;");
}
#endif

/** get current timestamp from OS */
static void my_gettime(MTime *t) {
    SYSTEMTIME s;

    /* GMT, not local time */
    GetSystemTime(&s);
    t->sec = s.wSecond;
    t->millisec = s.wMilliseconds;
}

/** sleep in current thread */
static void my_msleep(size_t millis) {
    SleepEx(millis, TRUE);
}

/* row diff, top-to-bottom, bottom is zero */
static void d2_spatdiff_down(
           mxSingle  *a,   /* 2d output, column order */
           const mxSingle *b,   /* 2d input, column order */
           size_t    m,     /* height */
           size_t    n      /* width */
           )
{
    int i;
    int max = (m*n) - 1;

    /* mexPrintf("tvden:spatdiff: m=%u n=%u\n", m, n); */

    /* out(row) = in(row+1) - in(row), bottom row is garbage */
    for(i = 0; i < max; i++) {
        a[i] = b[i+1] - b[i];
    }
    /* set bottom row to 0 */
    max = m*n;
    for(i = m - 1; i < max; i += m) {
        a[i] = 0;
    }
}

/* row diff, bottom-to-top, top is zero */
static void d2_spatdiff_up(
           mxSingle  *a,   /* 2d output, column order */
           const mxSingle *b,   /* 2d input, column order */
           size_t    m,     /* height */
           size_t    n      /* width */
           )
{
    int i;
    int max = m*n;

    /* out(row) = in(row) - in(r-1), top row is garbage */
    for(i = 1; i < max; i++) {
        a[i] = b[i] - b[i-1];
    }
    /* set top row to 0 */
    for(i = 0; i < max; i += m) {
        a[i] = 0;
    }
}

/* col diff, left-to-right, right col is zero */
static void d2_spatdiff_right(
           mxSingle *a,   /* 2d output, column order */
           const mxSingle *b,   /* 2d input, column order */
           size_t    m,     /* height */
           size_t    n      /* width */
           )
{
    int i;
    int max = m*(n-1);

    /* out(col) = in(col+1) - in(col), right col is garbage */
    for(i = 0; i < max; i++) {
        a[i] = b[i+m] - b[i];
    }
    /* set right col to 0 */
    max = m*n;
    for(i = m*(n-1); i < max; i++) {
        a[i] = 0;
    }
}

/* col diff, right-to-left, left col is zero */
static void d2_spatdiff_left(
           mxSingle *a,   /* 2d output, column order */
           const mxSingle *b,   /* 2d input, column order */
           size_t    m,     /* height */
           size_t    n      /* width */
           )
{
    int i;
    int max = m*n;

    /* out(col) = in(col) - in(col-1), left col is garbage */
    for(i = m; i < max; i++) {
        a[i] = b[i] - b[i-m];
    }
    /* set left col to 0 */
    for(i = 0; i < m; i++) {
        a[i] = 0;
    }
}

/* single - process a contiguous set of frames, specified by arg */
static void d3_spatdiff(
        mxSingle *a,            /* 3d out */
        const mxSingle *b,      /* 3d inp */
        size_t dir,             /* direction, 1..4 */
        size_t num_frm,         /* number of frames to process */
        size_t h,               /* frame height */
        size_t w                /* frame width */
        ) {
    size_t elem_per_frm = h * w;
    size_t i;

    /* i is frame number, 0-based */
    for(i = 0; i < num_frm; i++) {
        /* process one frame */
        switch(dir) {
            case 2:  d2_spatdiff_up   (a, b, h, w);  break;
            case 3:  d2_spatdiff_right(a, b, h, w);  break;
            case 4:  d2_spatdiff_left (a, b, h, w);  break;
            default: d2_spatdiff_down (a, b, h, w);  break;
        }
        a += elem_per_frm;
        b += elem_per_frm;
    } // for
}

/* worker */
static void d3_spatdiff_worker(ThrInp *arg) {
    const InpArg *inp = st->globInp;
    mxSingle *a = arg->a + (arg->first_frm * inp->elem2d);
    const mxSingle *b = arg->b + (arg->first_frm * inp->elem2d);
    size_t num_frm = arg->last_frm - arg->first_frm;

    d3_spatdiff(a, b, arg->dir, num_frm, inp->h, inp->w);
}

/* multi-thread dispatch */
static void d3_spatdiff_mt(
        mxSingle *a,        /* 3d out */
        const mxSingle *b,  /* 3d inp */
        size_t dir,         /* 1..4 */
        size_t nfrm,        /* ignored */
        size_t h,           /* ignored */
        size_t w            /* ignored */
        ) {
    const InpArg *inp = st->globInp;
    ThrInp *t;
    size_t frm = 0; /* current frame */
    size_t i;

    if (st->nthread < 1 || inp->nfrm < 1 || inp->frm_per_thread < 1) {
        mexPrintf("tvden:d3_spatdiff_mt: invalid parameter\n");
        return;
    }

    /* partition the work into groups of consecutive frames */
    for(i = 0; i < st->nthread; i++) {
        t = &st->thrInp[i];
        t->func = d3_spatdiff_worker;
        t->a = a;
        t->b = b;
        t->c = 0;   /* unused */
        t->d = 0;   /* unused */
        t->first_frm = frm;
        frm += inp->frm_per_thread;
        t->last_frm = frm;
        t->dir = dir;
    }
    /* assign left over frames to last thread */
    st->thrInp[st->nthread-1].last_frm = inp->nfrm;

    /* do the processing */
    exec_threads();
}

/* single; A = sum(B^2 + C^2, 3); */
static void d3_sum_squares(
        mxSingle *a,        /* 2d out */
        const mxSingle *b,  /* 3d inp */
        const mxSingle *c,  /* 3d inp */
        size_t work_cnt,    /* consecutive elem to process, <= elem_per_frm */
        size_t nfrm,
        size_t elem_per_frm
        ) {
    size_t i, j, n;

    /* do frame[0] */
    for(i = 0; i < work_cnt; i++) {
        a[i] = b[i]*b[i] + c[i]*c[i];
    }
    /* accumulate data from frame[1] to [N-1] */
    for(j = 1; j < nfrm; j++) {
        /* starting offset in multi-frame buffer */
        n = j * elem_per_frm;
        for(i = 0; i < work_cnt; i++, n++) {
            a[i] += b[n]*b[n] + c[n]*c[n];
        }
    }
}

/* single; A = 1 + D*sqrt(A) */
static void d2_sqroot_plus(
        mxSingle *a,          /* 2d in_out, column order */
        const mxSingle *d,    /* scalar single */
        size_t num_elem       /* elem to process */
        ) {
    size_t i;

    for(i = 0; i < num_elem; i++) {
        a[i] = 1 + d[0]*sqrtf(a[i]);
    }
}

/* single; A = 1 + D*sqrt(sum(B^2 + C^2, 3)) */
static void sum_of_sqr(
        mxSingle *a,         /* 2d out */
        const mxSingle *b,   /* 3d inp */
        const mxSingle *c,   /* 3d inp */
        const mxSingle *d,   /* scalar inp */
        size_t work_cnt,     /* must be <= elem_per_frm */
        size_t nfrm,
        size_t elem_per_frm
        ) {
    d3_sum_squares(a, b, c, work_cnt, nfrm, elem_per_frm);
    d2_sqroot_plus(a, d, work_cnt);
}

/* worker */
static void sum_of_sqr_worker(ThrInp *arg) {
    const InpArg *inp = st->globInp;
    /* in this function, these are col idx, not frame idx */
    size_t first_col = arg->first_frm;
    size_t last_col = arg->last_frm;
    mxSingle *a = arg->a + (first_col*inp->h);
    const mxSingle *b = arg->b + (first_col*inp->h);
    const mxSingle *c = arg->c + (first_col*inp->h);
    size_t work_cnt = (last_col - first_col) * inp->h;

    sum_of_sqr(a, b, c, arg->d, work_cnt, inp->nfrm, inp->elem2d);
}

/* multi-thread dispatch */
static void sum_of_sqr_mt(
        mxSingle *a,         /* 2d out */
        const mxSingle *b,   /* 3d inp */
        const mxSingle *c,   /* 3d inp */
        const mxSingle *d,   /* scalar inp */
        size_t work_cnt,     /* ignored */
        size_t nfrm,         /* ignored */
        size_t elem_per_frm  /* ignored */
        ) {
    const InpArg *inp = st->globInp;
    ThrInp *t;
    size_t col = 0; /* current col */
    size_t i;
    size_t col_per_thread = inp->w / st->nthread;

    if (st->nthread < 1 || inp->w < 1 || col_per_thread < 1) {
        mexPrintf("tvden:sum_of_sqr_mt: invalid parameter\n");
        return;
    }

    /* partition the work into vertical slices, down all frames */
    for(i = 0; i < st->nthread; i++) {
        t = &st->thrInp[i];
        t->func = sum_of_sqr_worker;
        t->a = a;
        t->b = b;
        t->c = c;
        t->d = d;
        t->first_frm = col;
        t->last_frm = col + col_per_thread;
        t->dir = 0;
        col += col_per_thread;
    }
    /* assign left over columns to last thread */
    st->thrInp[st->nthread-1].last_frm = inp->w;

    /* do the processing */
    exec_threads();
}
/* A = (A + (B*C)) / D */
static void d2_tvden(
           mxSingle  *a,        /* 2d output, column order */
           const mxSingle *b,   /* scalar */
           const mxSingle *c,   /* 2d input, column order */
           const mxSingle *d,   /* 2d input divisor */
           size_t num_elem      /* in frame */
           ) {
    int i;

    /* A = (A + (B*C)) / D */
    for(i = 0; i < num_elem; i++) {
        a[i] = (a[i] + (b[0] * c[i])) / d[i];
    }
}

/* single; A = (A + (B*C)) / D */
static void d3_tvden(
           mxSingle  *a,       /* 3d output, column order */
           const mxSingle *b,  /* scalar */
           const mxSingle *c,  /* 3d input, column order */
           const mxSingle *d,  /* 2d input divisor */
           size_t  num_elem,   /* in one frame */
           size_t  nfrm        /* number of frames */
        ) {
    size_t i;

    for(i = 0; i < nfrm; i++) {
        d2_tvden(a, b, c, d, num_elem);
        a += num_elem;
        c += num_elem;
    }
}

/* worker */
static void d3_tvden_worker(ThrInp *arg) {
    const InpArg *inp = st->globInp;
    mxSingle *a = arg->a + (arg->first_frm * inp->elem2d);
    const mxSingle *c = arg->c + (arg->first_frm * inp->elem2d);
    size_t num_frm = arg->last_frm - arg->first_frm;

    d3_tvden(a, arg->b, c, arg->d, inp->elem2d, num_frm);
}

/* multi-thread dispatch */
static void d3_tvden_mt(
        mxSingle *a,         /* 3d out */
        const mxSingle *b,   /* scalar inp */
        const mxSingle *c,   /* 3d inp */
        const mxSingle *d,   /* 2d inp */
        size_t elem_per_frm, /* ignored */
        size_t nfrm          /* ignored */
        ) {
    const InpArg *inp = st->globInp;
    ThrInp *t;
    size_t frm = 0; /* current frame */
    size_t i;

    if (st->nthread < 1 || inp->nfrm < 1 || inp->frm_per_thread < 1) {
        mexPrintf("tvden:d3_tvden_mt: invalid parameter\n");
        return;
    }

    /* partition the work into groups of consecutive frames */
    for(i = 0; i < st->nthread; i++) {
        t = &st->thrInp[i];
        t->func = d3_tvden_worker;
        t->a = a;
        t->b = b;
        t->c = c;
        t->d = d;
        t->first_frm = frm;
        frm += inp->frm_per_thread;
        t->last_frm = frm;
        t->dir = 0;
    }
    /* assign left over frames to last thread */
    st->thrInp[st->nthread-1].last_frm = inp->nfrm;

    /* do the processing */
    exec_threads();
}

/* single; A = B + C */
static void d3_mat_add(
        mxSingle *a,        /* 3d out */
        const mxSingle *b,  /* 3d inp */
        const mxSingle *c,  /* 3d inp */
        size_t num_elem
        ) {
    size_t i;

    for(i = 0; i < num_elem; i++) {
        a[i] = b[i] + c[i];
    }
}

/* worker */
static void d3_mat_add_worker(ThrInp *arg) {
    size_t elem2d = st->globInp->elem2d;
    mxSingle *a = arg->a + (arg->first_frm * elem2d);
    const mxSingle *b = arg->b + (arg->first_frm * elem2d);
    const mxSingle *c = arg->c + (arg->first_frm * elem2d);
    size_t num_cell = (arg->last_frm - arg->first_frm) * elem2d;

    d3_mat_add(a, b, c, num_cell);
}

/* multi-thread dispatch */
static void d3_mat_add_mt(
        mxSingle *a,        /* 3d out */
        const mxSingle *b,  /* 3d inp */
        const mxSingle *c,  /* 3d inp */
        size_t num_elem     /* ignored */
        ) {
    const InpArg *inp = st->globInp;
    ThrInp *t;
    size_t curfrm = 0; /* current frame */
    size_t i;

    /* partition the work into groups of consecutive frames */
    for(i = 0; i < st->nthread; i++) {
        t = &st->thrInp[i];
        t->func = d3_mat_add_worker;
        t->a = a;
        t->b = b;
        t->c = c;
        t->d = 0;
        t->first_frm = curfrm;
        t->last_frm = curfrm + inp->frm_per_thread;
        t->dir = 0;
        curfrm += inp->frm_per_thread;
    }
    /* assign left over frames to last thread */
    st->thrInp[st->nthread-1].last_frm = inp->nfrm;

    /* do the processing */
    exec_threads();
}

/**
 * Hand-coded 'for' loop in tv_denoise_oct(), avoiding data copies as
 * much as possible.
 *
 * @divp  output buffer
 * @inp   read-only input parameters
 */
static void process_input(mxSingle *divp, const InpArg *inp) {
    mxSingle *z;
    mxSingle *z1;
    mxSingle *z2;
    mxSingle *denom;
    mxSingle *p1;
    mxSingle *p2;
    mxSingle *p1_diff;
    mxSingle *p2_diff;
    size_t i;

    //z = zeros(size(f), 'single');
    z = mxCalloc(inp->elem3d, sizeof(mxSingle));
    //z1 = zeros(size(f), 'single');
    z1 = mxCalloc(inp->elem3d, sizeof(mxSingle));
    //z2 = zeros(size(f), 'single');
    z2 = mxCalloc(inp->elem3d, sizeof(mxSingle));
    //denom = zeros(N(1), N(2), 'single');
    denom = mxCalloc(inp->elem2d, sizeof(mxSingle));
    //p1 = zeros(size(f), 'single');
    p1 = mxCalloc(inp->elem3d, sizeof(mxSingle));
    //p2 = zeros(size(f), 'single');
    p2 = mxCalloc(inp->elem3d, sizeof(mxSingle));
    //p1_diff = zeros(size(f), 'single');
    p1_diff = mxCalloc(inp->elem3d, sizeof(mxSingle));
    //p2_diff = zeros(size(f), 'single');
    p2_diff = mxCalloc(inp->elem3d, sizeof(mxSingle));
    //divp = zeros(size(f), 'single');
    memset(divp, 0, inp->elem3d*sizeof(mxSingle));

//    for i=1:iters
    for(i = 0; i < inp->iter; i++) {
#if 1
//      %z = divp - f_lambda_s;
        d3_mat_add_mt(z, divp, inp->f_lambda, inp->elem3d);
//      %z1 = z(:,ir,:) - z;
        d3_spatdiff_mt(z1, z, 3, inp->nfrm, inp->h, inp->w);
//      %z2 = z(id,:,:) - z;
        d3_spatdiff_mt(z2, z, 1, inp->nfrm, inp->h, inp->w);
//      %denom = 1 + dt*sqrt(sum(z1.^2 + z2.^2,3));
        sum_of_sqr_mt(denom, z1, z2, inp->dt, inp->elem2d, inp->nfrm, inp->elem2d);
//      %p1 = (p1 + dt*z1) ./ denom;
        d3_tvden_mt(p1, inp->dt, z1, denom, inp->elem2d, inp->nfrm);
//      %p2 = (p2 + dt*z2) ./ denom;
        d3_tvden_mt(p2, inp->dt, z2, denom, inp->elem2d, inp->nfrm);
//      %p1_diff = p1 - p1(:,il,:);
        d3_spatdiff_mt(p1_diff, p1, 4, inp->nfrm, inp->h, inp->w);
//      %p2_diff = p2 - p2(iu,:,:);
        d3_spatdiff_mt(p2_diff, p2, 2, inp->nfrm, inp->h, inp->w);
//      %divp = p1_diff + p2_diff;
        d3_mat_add_mt(divp, p1_diff, p2_diff, inp->elem3d);
#endif
//    end
    } /* for */
    mxFree(z);
    mxFree(z1);
    mxFree(z2);
    mxFree(denom);
    mxFree(p1);
    mxFree(p2);
    mxFree(p1_diff);
    mxFree(p2_diff);
}

/** 
 * Worker thread that sleeps until it is signaled that there is work
 * to do. Main thread sends a broadcast to all worker threads.
 */
static unsigned int __stdcall worker_thread(void *args) {
    ThrInp *arg = args;
    size_t run;

    AcquireSRWLockExclusive(&st->mutex_wake);

    while(st->quit_flag == 0) {

        /* don't sleep if run flag has been set again */
        if (arg->running == 0) {
            SleepConditionVariableSRW(&st->cond_wake, &st->mutex_wake,
                INFINITE, 0 /*EXCLUSIVE*/);
        }
        /* note: mutex_wake is locked here */

        /* check for thread shutdown */
        if (st->quit_flag) {
            arg->running = 0;
            break;
        }

        /* main thread sets 'running' flag after arg is valid */
        if (arg->running) {
            /* unlock while working */
            ReleaseSRWLockExclusive(&st->mutex_wake);

            /* execute work function */
            if (arg->func) {
                (*arg->func)(arg);
            }
            //my_gettime(&arg->time_finish);

            AcquireSRWLockExclusive(&st->mutex_done);
            arg->running = 0;
            ReleaseSRWLockExclusive(&st->mutex_done);

            /* wake up main thread */
            WakeAllConditionVariable(&st->cond_done);

            /* lock mutex_wake again */
            AcquireSRWLockExclusive(&st->mutex_wake);
        }
    } /* while */

    ReleaseSRWLockExclusive(&st->mutex_wake);
    return 0;
}

/**
 * Launch a background thread.
 * @return 1 on success, 0 on error.
 */
static int start_thread(ThrInp *arg) {
    int rv = 1;
    size_t i;

    /* background thread */
    arg->handle = (HANDLE) _beginthreadex(0, 0, worker_thread,
        (void *)arg, 0, 0);
    if (arg->handle) {
        if (st->verbose) {
            mexPrintf("tvden: thread create handle=%x\n",
                    (size_t)arg->handle);
        }
    } else {
        rv = 0;
    }
    return rv;
}

/**
 * Wait for specified thread to shut down. Caller should already have
 * set quit_flag & woke up the threads.
 */
static void wait_for_thread(ThrInp *arg) {

    if (arg->handle) {
        WaitForSingleObject(arg->handle, INFINITE);
        CloseHandle(arg->handle);
        if (st->verbose) {
            mexPrintf("tvden: thread destroyed handle=%x\n", (size_t)arg->handle);
        }
    }
}

/**
 * Wake all threads. Each thread checks to see if there is work for it to
 * do. Pass 1 for normal dispatch, 0 to only wake-up threads (eg for
 * shutdown).
 */
static void wake_all_threads(size_t mode) {

    if (st->verbose) {
        mexPrintf("tvden: wake_all_threads\n");
        ioFlush();
    }

    if (st->nthread > 1) {
        if (mode) {
            size_t i;

            /* set run flag so worker thread knows ThrInp is good. run flag
             * is necessary because spurious wake-ups can happen. */
            AcquireSRWLockExclusive(&st->mutex_wake);
            for(i = 1; i < st->nthread; i++) {
                st->thrInp[i].running = 1;
            }
            ReleaseSRWLockExclusive(&st->mutex_wake);
        }
        /* send signal to wake up all threads */
        WakeAllConditionVariable(&st->cond_wake);
    }
    if (mode) {
        ThrInp *arg = &st->thrInp[0];

        /* run thread0 in main thread, run flag is not set */
        if (arg->func) {
            (*arg->func)(arg);
        }
    }
}

/**
 * Execute the worker threads, and wait for them to finish executing.
 *(to clear the run flag.) This function times out after ~60 seconds.
 */
static void exec_threads(void) {
    size_t i;
    size_t retry = 3;
    MTime t, t2;

    wake_all_threads(1);

    if (st->nthread < 2) {
        /* we're done */
        return;
    }

    if (st->verbose) {
        mexPrintf("tvden: begin exec_threads\n");
        ioFlush();
    }

    AcquireSRWLockExclusive(&st->mutex_done);

    /* we only get here if bg work was initiated */
    while(retry > 0) {
        BOOL n;

        /* check if worker(s) have already finished */
        for(i = 1; i < st->nthread; i++) {
            if (st->thrInp[i].running) {
                break;
            }
        } /* for */

        if (i == st->nthread) {
            /* all threads are finished */
            break;
        }

        /* unlock and wait for signal */
        n = SleepConditionVariableSRW(&st->cond_done, &st->mutex_done,
                20*1000 /*millisec or INFINITE */, 0 /*EXCLUSIVE*/);

        /* note: mutex_done is locked here */
        if (n == 0) {
            mexPrintf("tvden: warn: waiting for bg thread to complete\n");
            ioFlush();
            retry--;
        }
    } /* while */

    ReleaseSRWLockExclusive(&st->mutex_done);

    if (st->verbose) {
        mexPrintf("tvden: end wait_for_exec\n");
        ioFlush();
    }
}

/**
 * Create the background threads.
 *
 * @inp  read-only function inputs & parameters
 * @nthread  the total number of threads desired.
 *
 * It's possible user changes the number of threads on successive
 * calls; we optimize for case where it's the same.
 */
static void mt_setup(const InpArg *inp, size_t nthread) {

    st->globInp = inp;

    if (st->nthread > 0 && st->nthread != nthread) {
        /* user changed # of threads, so stop existing threads */
        mt_shutdown();
    }

    if (st->shutdown_registered == 0) {
        /* called if user unloads tvden func. from command window */
        mexAtExit(mt_shutdown);
        st->shutdown_registered = 1;
    }

    /* do not re-init state if > 0 */
    if (st->nthread == 0) {
        size_t i;

        InitializeConditionVariable(&st->cond_wake);
        InitializeSRWLock(&st->mutex_wake);

        InitializeConditionVariable(&st->cond_done);
        InitializeSRWLock(&st->mutex_done);

        for(i = 1; i < nthread; i++) {
            ThrInp *arg = &st->thrInp[i];
            arg->id = i;
            /* create a worker thread */
            start_thread(arg);
        }
        st->nthread = nthread;
    } else if (st->verbose) {
        mexPrintf("tvden: Reusing existing %u threads\n", st->nthread);
    }
    ioFlush();
}

/** destroy all the background threads */
static void mt_shutdown(void) {
    int i;

    if (st->nthread > 0) {
        st->quit_flag = 1;
        wake_all_threads(0);
        for(i = 1; i < st->nthread; i++) {
            /* wait for thread to exit */
            wait_for_thread(&st->thrInp[i]);
        }
        st->nthread = 0;
        st->quit_flag = 0;
    }
}

void mexFunction( int nlhs, mxArray *plhs[],
          int nrhs, const mxArray*prhs[] )

{
    const mxArray *a_out;   /* 3D output, single, MUST BE ALLOCATED */
    const mxArray *b_in;    /* iter, scalar input, single */
    const mxArray *c_in;    /* 3D input, single */
    const mxArray *d_in;    /* dt, scalar input, single */
    const mxArray *e_in;    /* ncpu, scalar, single */
    mxSingle *divp;
    const mxDouble *b;
    const mxSingle *c;
    const mxSingle *d;
    mwSize num_dim;
    const mwSize *dims;
    size_t i;
    size_t nthread;
    InpArg inp;

    if (nrhs != 5) {
        mexErrMsgIdAndTxt( "MATLAB:tvden:invalidNumInputs",
                "Exactly 5 input arguments required.");
        /* implied return */
    } else if (nlhs != 0) {
        mexErrMsgIdAndTxt( "MATLAB:tvden:maxlhs",
                "Exactly 0 output arguments required.");
        /* implied return */
    }

    a_out = (mxArray *)prhs[0];
    b_in = prhs[1];
    c_in = prhs[2];
    d_in = prhs[3];
    e_in = prhs[4];

    if( !mxIsSingle(a_out) || mxIsComplex(a_out) || mxIsSparse(a_out)) {
      mexErrMsgIdAndTxt( "MATLAB:tvden:invalidArg",
              "First argument must be a matrix of single.");
      /*implied return */
    }

    /* check to make sure the b argument is a scalar */
    if( !mxIsDouble(b_in) || mxIsComplex(b_in) || mxIsSparse(b_in)) {
      mexErrMsgIdAndTxt( "MATLAB:tvden:invalidArg",
              "Second argument must be a scalar double.");
      /*implied return */
    }

    /* check to make sure the c argument is a real matrix */
    if( !mxIsSingle(c_in) || mxIsComplex(c_in) || mxIsSparse(c_in)) {
      mexErrMsgIdAndTxt( "MATLAB:tvden:invalidArg",
              "Third argument must be a matrix of single.");
      /*implied return */
    }

    /* check to make sure the d argument is a real matrix */
    if( !mxIsSingle(d_in) || mxIsComplex(d_in) || mxIsSparse(d_in)) {
      mexErrMsgIdAndTxt( "MATLAB:tvden:invalidArg",
              "Fourth argument must be a matrix of single.");
      /*implied return */
    }

    /* check to make sure the e argument is a scalar */
    if( !mxIsDouble(e_in) || mxIsComplex(e_in) || mxIsSparse(e_in)) {
      mexErrMsgIdAndTxt( "MATLAB:tvden:invalidArg",
              "Fifth argument must be a scalar double.");
      /*implied return */
    }

    /* 2 or 3 */
    num_dim = num_dim = mxGetNumberOfDimensions(a_out);

    /* get length of each dimension */
    dims = mxGetDimensions(a_out);

    {
        /* range check c_in */
        const size_t *tmp = mxGetDimensions(c_in);
        for(i = 0; i < num_dim; i++) {
            if (tmp[i] != dims[i]) {
                mexErrMsgIdAndTxt( "MATLAB:tvden:invalidArg",
                    "Third argument dims do not match first.");
            }
            /*implied return */
        }
    }

    #if 0
    st->verbose = 1;
    #endif

    /* data dimensions - may be 0 */
    inp.h = dims[0];
    inp.w = dims[1];
    inp.nfrm = (num_dim > 2) ? dims[2] : 1 ;
    inp.elem2d = inp.h * inp.w;
    inp.elem3d = inp.elem2d * inp.nfrm;
    if (inp.elem3d < 1) {
        /* nothing to do */
        mexPrintf("tvden: input array has no data\n");
        return;
    }

    nthread = (size_t)mxGetDoubles(e_in)[0];
    /* range check */
    if (nthread < 1) nthread = 1;
    else if (nthread > MAXCPU) nthread = MAXCPU;

    if (inp.nfrm <= nthread) {
        nthread = 1;
    }
    inp.frm_per_thread = inp.nfrm / nthread;

    if (st->verbose) {
        mexPrintf("tvden: elem3d=%u nfrm=%u nthread=%u\n", inp.elem3d,
            inp.nfrm, nthread);
    }

    /* get ptr to data array */
    divp = mxGetSingles(a_out);
    inp.iter = mxGetDoubles(b_in)[0];
    inp.f_lambda = mxGetSingles(c_in);
    inp.dt = mxGetSingles(d_in);

    mt_setup(&inp, nthread);

    /* write result to divp */
    process_input(divp, &inp);
    ioFlush();
}
