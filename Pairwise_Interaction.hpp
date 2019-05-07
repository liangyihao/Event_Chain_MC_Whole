#ifndef PAIRWISE_INTERACTION_ECMC
#define PAIRWISE_INTERACTION_ECMC
#include"basic.hpp"
double Event_Time_Hard_Sphere(double4 X1,double4 X2,int axis_index,double*Params,double Max_Event_T);
double Event_Time_Spring(double4 X1,double4 X2,int axis_index,double*Params,double Max_Event_T);

void Create_Hard_Sphere_Interaction_Between_Types(int Type_id1,int Type_id2,double d);
void Create_Spring_Interaction_Between_Beads(int2 Bead_id1,int2 Bead_id2,double k,double l0);
#endif
