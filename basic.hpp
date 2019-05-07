/*
Author: Yihao Liang
liangyihaosjtu@gmail.com
This code is for Event Chain Monte Carlo for pairwise interacting many body system
*/

#ifndef BASIC_ECMC
#define BASIC_ECMC

#include <vector>
#include "public.hpp"
#include "CellList.hpp"

using namespace std;

class Parameter_List{
    public:
    double*data;
    //data[0]: Lx
    //data[1]: Ly
    //data[2]: Lz
    ///...
    //When calling function for Event Time, just pass the address of data[0]
    Parameter_List(double Lx,double Ly,double Lz) {
        data=new double [8];
        data[0]=Lx;
        data[1]=Ly;
        data[2]=Lz;
    }
};

extern vector<Bead_Type> Types;
extern vector<double(*)(double4,double4,int,double*,double)> Event_Time_Generator_List;
extern vector<Parameter_List> Param_Lists;
extern int2 Active_Bead;
extern double MAX_SHORT_INTERACTION_RANGE;
extern CellList* Global_Cell_List_Pointer;
void Monte_Carlo(double Stop_Clock,int axis);
#endif
