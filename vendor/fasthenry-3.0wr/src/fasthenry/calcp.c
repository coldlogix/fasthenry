/* this calls the routine to calculate the filament-filament
   interaction exactly */

#include "calcp.h"
#include "mutual.h"

static int num2nd=0, num4th=0, numexact=0;
static int num2ndsav=0, num4thsav=0, numexactsav=0;

/* SRW */

double calcp(charge *pchg1, charge *pchg2, double *pfd)
/* double *pfd;   left over from fastcap */
{

  if (pfd != NULL)
    fprintf(stderr, "calcp: I don't know what to do with pfd!=NULL\n");

  if (pchg1->fil->filnumber == pchg2->fil->filnumber)
    /* self term */
#if SUPERCON == ON
    { double tmp = selfterm(pchg1->fil);
      struct Segment *seg1 = pchg1->fil->segm;
      if (seg1->lambda != 0.0)
        tmp += seg1->r2*pchg1->fil->length/pchg1->fil->area;
      return (tmp);
    }
#else
    return selfterm(pchg1->fil);
#endif
  else
    /* calculate mutual inductance of the two filaments */
    return mutual(pchg1->fil, pchg2->fil);
}

/* from the fastcap calcp */
void dumpnums(int flag, int size)
{
  double total;

  if(flag == ON) {		/* if first call */
    num2ndsav = num2nd;
    num4thsav = num4th;
    numexactsav = numexact;
  }
  else {
    total = num2ndsav + num4thsav + numexactsav;
#if MULDAT == ON
    fprintf(stdout, "Potential coefficient counts\n multipole only:\n");
    fprintf(stdout,
	    "  2nd order: %d %.3g%%; 4th: %d %.3g%%; Integral: %d %.3g%%\n",
	    num2nd, 100*(num2ndsav/total), num4th, 100*(num4thsav/total),
	    numexact, 100*(numexactsav/total));
#endif
    total = num2nd + num4th + numexact;
#if MULDAT == ON
    fprintf(stdout, " multipole plus adaptive:\n");
    fprintf(stdout,
	    "  2nd order: %d %.3g%%; 4th: %d %.3g%%; Integral: %d %.3g%%\n",
	    num2nd, 100*(num2nd/total), num4th, 100*(num4th/total),
	    numexact, 100*(numexact/total));
#endif
    fprintf(stdout, "Percentage of multiplies done by multipole: %.3g%%\n",
	    100*(size*size - total)/(size*size));
    if(size*size == total)
	fprintf(stdout, "Warning: no multipole acceleration\n");
  }
}

double tilelength(charge *nq)
{
  return nq->max_diag;
}
