#ifndef __findpaths_H__
#define __findpaths_H__

#include "induct.h"

/* findpaths.c */
NODES *getrealnode(NODES*);
NODES *getothernode(NODES*, seg_ptr);
int is_normal(NODES*);
int is_gp_node(NODES*);
int is_node_in_list(NODES*, NPATH*);
int is_orignode_in_list(NODES*, NPATH*);
NPATH *add_node_to_list(SYS*, NODES*, NPATH*);
NODES *get_node_from_name(char*, SYS*);
int equivnodes(char*, SYS*);
PSEUDO_NODE *create_pn(SYS* sys, char*, NODES*);
void make_equiv(NODES*, NODES*);
void append_pnlist(PSEUDO_NODE*, SYS*);
void add_to_connected_segs(SYS*, NODES*, SEGMENT*, PSEUDO_SEG*);
void remove_from_connected_segs(SYS*, NODES*, SEGMENT*, PSEUDO_SEG*);
double mag(double, double, double);
double dotp(double, double, double, double, double, double);
NODES *find_nearest_gpnode(double, double, double, GROUNDPLANE*, int*, int*);
SPATH *add_seg_to_list(SYS*, seg_ptr, SPATH*);
PSEUDO_SEG *make_pseudo_seg(SYS*,NODES*, NODES*, char);
SPATH *make_new_fake_segs(SYS*, NODES*, NPATH*, SPATH*);
EXTERNAL *add_to_external_list(EXTERNAL*, EXTERNAL*);
EXTERNAL *make_external(SYS*, PSEUDO_SEG*, int, char*, char*, char*);
EXTERNAL *get_external_from_portname(char*, SYS*);
EXTERNAL *get_next_ext(EXTERNAL*);
void make_trees(SYS*);
int count_tree_meshes(TREE*);
int count_externals(EXTERNAL*);
void find_hole_meshes(SYS*);

#endif
