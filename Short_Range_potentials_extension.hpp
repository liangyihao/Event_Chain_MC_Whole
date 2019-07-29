#ifndef LJ_POTENTIAL_ECMC
#define LJ_POTENTIAL_ECMC

#include "basic.hpp"

double Event_Time_LJ_Potential(double4 X1, double4 X2, int axis_index, double*Params, double Max_Event_T);

void Create_LJ_Interaction_Between_Types(int Type_id1, int Type_id2, bool Using_CellList1, bool Using_CellList2, double sigma, double epsilon, double rcut);

double Event_Time_Gauss_Potential(double4 X1, double4 X2, int axis_index, double*Params, double Max_Event_T);

void Create_Gauss_Interaction_Between_Types(int Type_id1, int Type_id2, bool Using_CellList1, bool Using_CellList2, double sigma, double epsilon, double rcut);
#endif