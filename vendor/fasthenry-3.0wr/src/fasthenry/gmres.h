#ifndef __gmres_H__
#define __gmres_H__

#include "induct.h"
#include "sparse/spMatrix.h"

typedef CX (*gmres_cb1)(CX*, CX*, int);
typedef void (*gmres_cb2)(CX*, ssystem*, CX*, int, charge*, double, double*,
    SYS*);
void cx_invert(CX**, int);
int gmres(CX**, CX*, CX*, gmres_cb1, gmres_cb2, int, int, double, ssystem*,
     charge*, double, double*, SYS*, int);
CX inner(CX*, CX*, int);
void directmatvec(CX*, ssystem*, CX*, int, charge*, double, double*, SYS*);

#endif
