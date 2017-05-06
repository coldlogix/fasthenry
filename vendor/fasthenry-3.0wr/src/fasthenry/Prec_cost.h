#ifndef __Prec_cost_H__
#define __Prec_cost_H__

#include "induct.h"

/* Prec_cost.c */
double OneCubeCost(cube*****, int, int, int, int, int, double*);
double ratio_of_divided_segs(double, charge*, SYS*);
int is_gp_charge(charge*);
// void add_to_counts(cube*, int, int*****, int*****);

#endif
