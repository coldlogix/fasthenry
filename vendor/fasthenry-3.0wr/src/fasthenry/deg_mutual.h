#ifndef __deg_mutual_H__
#define __deg_mutual_H__

#include "induct.h"

/* deg_mutual.c */
enum degen_type find_deg_dims(FILAMENT*);
double compute_for_degenerate(FILAMENT*, FILAMENT*, int, double*, double*,
    enum degen_type, enum degen_type, double);
// void setup_tape_to_tape(FILAMENT*, FILAMENT*, int, double*, double*,
//     enum degen_type, enum degen_type, FILAMENT*, FILAMENT*, double**, double**);
// double do_tape_to_brick(FILAMENT*, FILAMENT*, int, double*, double*,
//     enum degen_type, enum degen_type);

#endif
