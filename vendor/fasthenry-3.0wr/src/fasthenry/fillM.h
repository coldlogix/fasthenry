#ifndef __fillM_H__
#define __fillM_H__

#include "induct.h"

void fillM(SYS*);
MELEMENT *insert_in_list(SYS*, MELEMENT*, MELEMENT*);
MELEMENT *make_melement(SYS*, int, FILAMENT*, int);
void fill_b(EXTERNAL*, CX*);
void extractYcol(CX**, CX*, EXTERNAL*, EXTERNAL*);
char *get_a_name(PSEUDO_SEG*);
void makegrids(SYS*, CX*, int, int);

#endif
