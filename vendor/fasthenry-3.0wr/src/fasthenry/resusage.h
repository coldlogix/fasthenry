/* # ***** sort to /src/header
   # ***** */
/* header where rusage and time structs are defined */

/* SRW ** Revised this:
 * getrusage is preferred, as it computes cpu time rather than wall-clock
 * time.  If getrusage is not available (as in MinGW) use gettimeofday.
 */

#include <sys/time.h>

static double dtime = 0.0;
static long sectime, utime;

#ifndef NO_RUSAGE
#include <sys/resource.h>

struct rusage timestuff;

/* SRW */
#ifdef SRW0814
#define SRWSECONDS \
static double srw_seconds() \
{ \
    getrusage(RUSAGE_SELF, &timestuff); \
    return (timestuff.ru_utime.tv_sec + 1.0e-6*timestuff.ru_utime.tv_usec); \
}
#endif

#define starttimer getrusage(RUSAGE_SELF, &timestuff); \
    sectime = timestuff.ru_utime.tv_sec; \
    utime = timestuff.ru_utime.tv_usec

#define stoptimer getrusage(RUSAGE_SELF, &timestuff); \
    dtime = (double)(timestuff.ru_utime.tv_sec - sectime) \
        + 1.0e-6*(double)(timestuff.ru_utime.tv_usec - utime)

#else /* NO_RUSAGE */

static struct timeval ru_tv;
static struct timezone ru_tz;

/* SRW */
#ifdef SRW0814
#define SRWSECONDS \
static double srw_seconds() \
{ \
    gettimeofday(&ru_tv, &ru_tz); \
    return (ru_tv.tv_sec + 1.0e-6*ru_tv.tv_usec); \
}
#endif

#define starttimer gettimeofday(&ru_tv, &ru_tz); \
    sectime = ru_tv.tv_sec; \
    utime = ru_tv.tv_usec

#define stoptimer gettimeofday(&ru_tv, &ru_tz); \
    dtime = (double)(ru_tv.tv_sec - sectime) \
        + 1.0e-6*(double)(ru_tv.tv_usec - utime)

#endif /* NO_RUSAGE */

#define DUMPRSS			/*  */


#ifdef notdef
/* SRW ** Below find original code. */
#ifdef FOUR
#define NOTOTHER 1
#include <sys/time.h>
#include <sys/resource.h>
struct rusage timestuff;
#endif

#ifdef FIVE
#define NOTOTHER 1
#include <sys/types.h>
#include <sys/param.h>
#include <sys/times.h>
struct tms timestuff;
#endif

/* define macros for time and resident memory usage checks */

static double dtime = 0.0;
static long sectime, utime;

#ifdef NOTOTHER

#ifdef FOUR			/* 4.2,3BSD (tested: Sun4, IBM6000, DEC5000) */
/* SRW */
#ifdef SRW0814
#define SRWSECONDS \
static double srw_seconds() \
{ \
    getrusage(RUSAGE_SELF, &timestuff); \
    return (timestuff.ru_utime.tv_sec + 1.0e-6*timestuff.ru_utime.tv_usec); \
}
#endif

#define starttimer getrusage(RUSAGE_SELF, &timestuff); \
sectime = timestuff.ru_utime.tv_sec; \
utime = timestuff.ru_utime.tv_usec
#define stoptimer getrusage(RUSAGE_SELF, &timestuff); \
dtime = (double)(timestuff.ru_utime.tv_sec - sectime) \
        + 1.0e-6*(double)(timestuff.ru_utime.tv_usec - utime)
#define DUMPRSS			/*  */
#endif /* FOUR */

#ifdef FIVE			/* for System V (tested: HP300) */
/* SRW */
#ifdef SRW0814
#define SRWSECONDS \
static double srw_seconds() \
{ \
    times(&timestuff); \
    return (timestuff.tms_utime/(double)HZ); \
}
#endif

#define starttimer times(&timestuff); \
utime = timestuff.tms_utime
#define stoptimer times(&timestuff); \
dtime = (timestuff.tms_utime)-utime; \
dtime /= HZ
#define DUMPRSS			/*  */
#endif /* FIVE */

#else				/* default - no timers */

#define starttimer		/*  */
#define stoptimer		/*  */
#define DUMPRSS			/*  */

#endif /* NOTOTHER */
#endif /* notdef */
