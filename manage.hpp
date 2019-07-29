/*
Author: Yihao Liang
liangyihaosjtu@gmail.com
This code is for Event Chain Monte Carlo for pairwise interacting many body system
*/
/*
Used to manage system variables
*/

#ifndef MANAGE
#define MANAGE
#include "basic.hpp"
#include "CellVetoList.hpp"
extern double Lx,Ly,Lz;
extern vector<double>valence_of_type;
extern vector<double>valence_list;//record all the non-zero valence of beads
extern vector<CellVetoList*> Cell_Veto_Lists;//for all valences, you need to create a Cell-Veto List

int Create_Type();//Create Type, return it's id
int Create_Bead(int type_id,double4 x);//Create Bead, return its secondary id. If failed, return -1.
int Create_Charged_Type(double valence);//Create Charged Type, return it's id
int Create_Charged_Bead(int type_id,double4 x);//Create Charged Bead, return its secondary id. If failed, return -1.
void initialize_sys_charge();
void Run(char*InputFileName);
#endif
