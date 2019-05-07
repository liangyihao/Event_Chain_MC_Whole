#ifndef EWALD_SUM_FACTOR_HPP
#define EWALD_SUM_FACTOR_HPP
#include "public.hpp"

class Ewald_Sum{
    private:
        double alpha;
        double Lx,Ly,Lz;
        double B[100][100][100];
        int cut_off_n,cut_off_m;
        double accuracy;
        //For computation of event rate, D_Potential, See Paper: M. F. Faulkner, L. Qin, A. C. Maggs, W. Krauth All-atom computations with irreversible Markov chains Journal of Chemical Physics 149, 064113 (2018)
    public:
        double Potential(double4 X,double4 Y);
        double Self_E(double4 X);
        double D_Potential(double4 X,double4 Y,int axis);
        void Reset(double Lx,double Ly,double Lz,double accuracy);
        Ewald_Sum(double Lx,double Ly,double Lz,double accuracy);
        //double D_Potential_Check(double4 X,double4 Y,int axis);
        double Event_Time_Colomb(double4 X1,double4 X2,int axis_index,double Bjerrum_Length, double Max_Event_T);
};
#endif
