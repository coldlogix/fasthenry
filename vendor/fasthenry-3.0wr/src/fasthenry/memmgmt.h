#ifndef __memmgmt_H__
#define __memmgmt_H__

#ifdef __memmgmt_C__
#define EXTERN
#else
#define EXTERN extern
#endif

#include <stdio.h>


/***********************************************************************
  macros for allocation with checks for NULL pntrs and 0 byte requests
  - also keep an allocated memory count
  - CALLOC() is used when the memory must be zeroed
    its core should be either calloc() or ualloc() (not as fast but
    more space efficient, no free list - uses sbrk() and never frees)
  - MALLOC() used when memory can be anything
    core should be malloc() or ualloc()
***********************************************************************/

/* SRW - Default now is to use calloc/malloc rather than ualloc.  The
 * sbrk function is deprecated in many operating systems. 11/16/13
 */

/* counts of memory usage by multipole matrix type */


#ifdef MATTDEBUG
extern long membins[1001];
#endif

/* types of memory usage by multipole matrix type */
#define AQ2M 0
#define AQ2L 1
#define AQ2P 2
#define AL2L 3
#define AM2M 4
#define AM2L 5
#define AM2P 6
#define AL2P 7
#define AQ2PD 8
#define AMSC 9
#define IND 10

EXTERN unsigned long memAREA[12];
#define memAREA_Elements (sizeof(memAREA)/sizeof(unsigned long))

#define memcount (memAREA[0])
#define memQ2M (memAREA[1+AQ2M])
#define memQ2L (memAREA[1+AQ2L])
#define memQ2P (memAREA[1+AQ2P])
#define memL2L (memAREA[1+AL2L])
#define memM2M (memAREA[1+AM2M])
#define memM2L (memAREA[1+AM2L])
#define memM2P (memAREA[1+AM2P])
#define memL2P (memAREA[1+AL2P])
#define memQ2PD (memAREA[1+AQ2PD])
#define memMSC (memAREA[1+AMSC])
#define memIND (memAREA[1+IND])

#define ON 1
#define OFF 0

void memmgmt_init();

void ALLOC_func(int, void**,int,int,int,int, const char*, int);
#define CALLOC(PNTR,NUM,TYP,FLAG,MTYP) ALLOC_func(0, (void**)&PNTR,NUM,sizeof(TYP),FLAG,MTYP,__FILE__,__LINE__ )
#define MALLOC(PNTR,NUM,TYP,FLAG,MTYP) ALLOC_func(1, (void**)&PNTR,NUM,sizeof(TYP),FLAG,MTYP,__FILE__,__LINE__ )

void DUMPALLOCSIZ(FILE*);

char* MattAlloc_func(int, int, const char*, int);
#define MattAlloc(number,size) MattAlloc_func(number,size,__FILE__,__LINE__)

#undef EXTERN

#endif
