#ifndef __Precond_H__
#define __Precond_H__

#include "induct.h"


unsigned long PreCondSize(PRE_ELEMENT** Precond, int size);
void indPrecond(ssystem*, SYS*, double);
void multPrecond(PRE_ELEMENT**, CX*, CX*, int);
// MELEMENT *getnext(MELEMENT*, int*);
// void cx_invert_dup(CX**, int, DUPS*);
// void mark_dup_mesh(MELEMENT**, int*, int, DUPS*, int*);
// void dumpPrecond(PRE_ELEMENT**, int, char*);
void indPrecond_direct(ssystem*, SYS*, double);

#endif
