#define __memmgmt_C__

#include "memmgmt.h"
#include <stdlib.h>

#include "uglieralloc.h"


#define NO_SBRK
#ifdef NO_SBRK
#define sbrk(x) 0L
#define CALCORE(NUM, TYPE) calloc(NUM,TYPE)
#define MALCORE(NUM, TYPE) malloc(NUM*TYPE)
#else
#define CALCORE(NUM, TYPE) ualloc(((NUM)*(TYPE))
#define MALCORE(NUM, TYPE) ualloc(((NUM)*(TYPE))
#endif


char* MattAlloc_func(int number, int size, const char* filename, int linenumber)
{
  char *blah;

  blah=NULL;
  ALLOC_func(0, (void**)&blah,number*size,sizeof(char),ON,IND,filename,linenumber );

  return blah;
}

void memmgmt_init()
{
    int i;
    for (i=0; i<memAREA_Elements ;i++)
    {
        memAREA[i]=0;
    }
}

void DUMPALLOCSIZ(FILE* fop)
{
  fprintf(fop,
		"Total Memory Allocated: %ld kilobytes (brk = 0x%lx)\n",
		memcount/1024, (long)sbrk(0));
  fprintf(fop, " Q2M  matrix memory allocated: %7.ld kilobytes\n",
		memQ2M/1024);
  memcount = memQ2M;
  fprintf(fop, " Q2L  matrix memory allocated: %7.ld kilobytes\n",
		memQ2L/1024);
  memcount += memQ2L;
  fprintf(fop, " Q2P  matrix memory allocated: %7.ld kilobytes\n",
		memQ2P/1024);
  memcount += memQ2P;
  fprintf(fop, " L2L  matrix memory allocated: %7.ld kilobytes\n",
		memL2L/1024);
  memcount += memL2L;
  fprintf(fop, " M2M  matrix memory allocated: %7.ld kilobytes\n",
		memM2M/1024);
  memcount += memM2M;
  fprintf(fop, " M2L  matrix memory allocated: %7.ld kilobytes\n",
		memM2L/1024);
  memcount += memM2L;
  fprintf(fop, " M2P  matrix memory allocated: %7.ld kilobytes\n",
		memM2P/1024);
  memcount += memM2P;
  fprintf(fop, " L2P  matrix memory allocated: %7.ld kilobytes\n",
		memL2P/1024);
  memcount += memL2P;
  fprintf(fop, " Q2PD matrix memory allocated: %7.ld kilobytes\n",
		memQ2PD/1024);
  memcount += memQ2PD;
  fprintf(fop, " Miscellaneous mem. allocated: %7.ld kilobytes\n",
		memMSC/1024);

  memcount += memMSC;
  fprintf(fop, " Inductance mem. allocated: %7.ld kilobytes\n",
		memIND/1024);
  memcount += memIND;
  fprintf(fop, " Total memory (check w/above): %7.ld kilobytes\n",
		memcount/1024);

}

void ALLOC_func(int alloctype, void** PNTR,int NUM,int TYPESIZE,int FLAG ,int MTYP,const char* filename, int linenumber)
{
  if((NUM)*TYPESIZE==0)
  {
    fprintf(stderr,
		     "zero element request in file `%s' at line %d\n",
		     __FILE__, __LINE__);
	return;
  }

  /* Check if there is already an allocation */
  if ((*PNTR)!=NULL)
  {
    free(*PNTR);
    (*PNTR)=NULL;
  }

  /* Allcate memory */
  if (alloctype==0)
  {
    (*PNTR)=(void*)CALCORE(NUM, TYPESIZE);
  }
  else
  {
    (*PNTR)=(void*)MALCORE(NUM, TYPESIZE);
  }

  /* Error handling in case memory was not granted */
  if((*PNTR)==0)
  {
    fprintf(stderr,
            "\nOut of memory in file `%s' at line %d\n",
                    filename, linenumber);
    fprintf(stderr, " (NULL pointer on %ld byte request)\n",
                 (NUM)*TYPESIZE);
                DUMPALLOCSIZ(stderr);

    fflush(stderr);
    fflush(stdout);
    if(FLAG == ON)
    {
      exit(1);
    }
    return;
  }

  /* Update statistics */
  memcount += ((NUM)*TYPESIZE);
  if ( (MTYP+1)< memAREA_Elements  )
  {
    memAREA[MTYP+1]+= ((NUM)*TYPESIZE);
  }
  else
  {
    fprintf(stderr, "ALLOC: unknown memory type %d\n", MTYP);
    exit(1);
  }
}

