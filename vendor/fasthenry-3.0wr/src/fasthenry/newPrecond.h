#ifndef __newPrecond_H__
#define __newPrecond_H__

#include "induct.h"

void choose_and_setup_precond(SYS*);
// void get_selfs(SYS*);
void fill_spPre(ssystem*, SYS*, double);
void create_sparMatrix(SYS*);
// void fill_bySegment(ssystem*, SYS*, double);
// void fill_diagL(ssystem*, SYS*, double);
void fill_diagR(SYS*);
double shift_mutual(FILAMENT*, FILAMENT*, ssystem*);

#endif
