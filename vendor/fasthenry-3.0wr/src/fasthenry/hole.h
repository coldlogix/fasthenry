#ifndef __hole_H__
#define __hole_H__

#include "induct.h"

/* hole.c */
int is_next_word(char*, char*);
int is_hole(NODES*);
HoleList *make_holelist(SYS*, HoleList*, char*, double, double, double, double, int*);

int skipspace(char*);
int eos(char);
// void hole_error(char*, char*, HoleList*);
int is_one_of(char, char*);
// void delete_node(NODES*);
void make_holes(HoleList*, GROUNDPLANE*);
// void hole_point(HoleList*, GROUNDPLANE*, double, double, double, double);
// void hole_rect(HoleList*, GROUNDPLANE*, double, double, double, double);
// void hole_circle(HoleList*, GROUNDPLANE*, double, double, double, double);
// void hole_user1(HoleList*, GROUNDPLANE*, double, double, double, double);
// void hole_user2(HoleList*, GROUNDPLANE*, double, double, double, double);
// void hole_user3(HoleList*, GROUNDPLANE*, double, double, double, double);
// void hole_user4(HoleList*, GROUNDPLANE*, double, double, double, double);
// void hole_user5(HoleList*, GROUNDPLANE*, double, double, double, double);
// void hole_user6(HoleList*, GROUNDPLANE*, double, double, double, double);
// void hole_user7(HoleList*, GROUNDPLANE*, double, double, double, double);

#endif
