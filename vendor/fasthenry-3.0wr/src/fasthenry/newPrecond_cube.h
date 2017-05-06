#ifndef __newPrecond_cube_H__
#define __newPrecond_cube_H__

#include "induct.h"

/* newPrecond_cube.c */
PRE_ELEMENT** initPrecond(int);
void indPrecond(ssystem*, SYS*, double);
void multPrecond(PRE_ELEMENT**, CX*, CX*, int);
// MELEMENT *getnext(MELEMENT*, int*);
// void cx_invert_dup(CX**, int, DUPS*);
// void mark_dup_mesh(MELEMENT**, int*, int, DUPS*, int*);
// void dumpPrecond(PRE_ELEMENT**, int);
void indPrecond_direct(ssystem*, SYS*, double);

#endif
