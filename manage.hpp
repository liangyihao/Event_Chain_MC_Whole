/*
Used to manage systematic variables
*/

#ifndef MANAGE
#define MANAGE
#include "basic.hpp"
extern double Lx,Ly,Lz;

int Create_Type();//Create Type, return it's id
int Create_Bead(int type_id,double4 x);//Create Bead, return its secondary id. If failed, return -1.

#endif
