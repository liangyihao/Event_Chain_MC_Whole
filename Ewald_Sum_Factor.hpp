/*
Author: Yihao Liang
liangyihaosjtu@gmail.com
This code is for Event Chain Monte Carlo for pairwise interacting many body system
*/
#ifndef EWALD_SUM_FACTOR_HPP
#define EWALD_SUM_FACTOR_HPP
#include "public.hpp"

class Ewald_Sum{
    private:
        double alpha;
        double Lx,Ly,Lz;
        double B[30][30][300];
        
        double C[30][30][30];//Accelarated coefficients

        int cut_off_n,cut_off_m;
        double accuracy;
        //For computation of event rate, D_Potential, See Paper: M. F. Faulkner, L. Qin, A. C. Maggs, W. Krauth All-atom computations with irreversible Markov chains Journal of Chemical Physics 149, 064113 (2018)
    public:

        double Potential(double4 X,double4 Y);
        double Self_E(double4 X);
        double D_Potential(double4 X,double4 Y,int axis);
        //double4 Gradient_D_Potential(double4 X,double4 Y,int axis);
        void Reset(double Lx,double Ly,double Lz,double accuracy);
        Ewald_Sum(double Lx,double Ly,double Lz,double accuracy);
        
        double D_Potential_Check(double4 X,double4 Y,int axis);

        double Event_Time_Colomb(double4 X1,double4 X2,int axis_index,double Bjerrum_Length, double Max_Event_T);


        double get_accuracy(){return accuracy;}
        double get_Lx(){return Lx;}
        double get_Ly(){return Ly;}
        double get_Lz(){return Lz;}
};
#endif
