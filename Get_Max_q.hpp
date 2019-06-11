/*
Author: Yihao Liang
liangyihaosjtu@gmail.com
This code is for Event Chain Monte Carlo for pairwise interacting many body system
*/
#ifndef GET_MAX_Q
#define GET_MAX_Q
#include "Ewald_Sum_Factor.hpp"
double get_max_q(double xb1, double xb2, double yb1, double yb2, double zb1, double zb2, Ewald_Sum&ES,int axis);

double get_min_q(double xb1, double xb2, double yb1, double yb2, double zb1, double zb2, Ewald_Sum&ES,int axis);

#endif
