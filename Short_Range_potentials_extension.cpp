// Author: Boran Ma
// Email: boran.ma@northwestern.edu

// pairwise interaction described by (shifted) L-J potential
// pairwise interaction described by Gaussian potential

// Modified by Yihao Liang

#include <vector>
#include <cmath>
#include <iostream>
#include "Random_Number.hpp"
#include "basic.hpp"
#include "Short_Range_potentials_extension.hpp"
#include "manage.hpp"
using namespace std;

inline double distance_sqr_Calc(double s, double drr2, double dxx, double Lx)
{
	double distOri2 = (dxx-s)*(dxx-s)+drr2;
	double distIma2 = (Lx+dxx-s)*(Lx+dxx-s)+drr2;
	return min(distOri2, distIma2);
}

inline double VLJCalc(double epsilon, double sigma, double r_sqr)
{
	double temp;
	temp=pow(sigma*sigma/r_sqr,3.0);
	return 4.0*epsilon*(temp*temp-temp);
}


double Event_Time_LJ_Potential(double4 X1, double4 X2, int axis_index, double *Params, double Max_Event_T)
{
    //Move X1 along axis_index direction, get reject event time of spring interaction with X2, with parameters Params.
    //If the event time is less than Max_Event_T and less than period in the moving direction, return real event time.
    //If the event time is greater than Max_Event_T, return 2*(period of moving direction), to reduce the time cost. However, in this function, we just return real value.
    //If the event time is greater than one period of moving direction, return 2*(period of moving direction)

	double Lx, Ly, Lz, sigma, epsilon, rcut;

	Lx = Params[0];
	Ly = Params[1];
	Lz = Params[2];
	sigma = Params[3];
	epsilon = Params[4];
	rcut = Params[5];

	double rmin = pow(2.0, 1.0/6.0)*sigma;
	double shift = 4.0*epsilon*(pow(sigma/rcut, 12.0)-pow(sigma/rcut, 6.0));

	double dxx, dyy, dzz;
	double Lxx, Lyy, Lzz;

	int sign = axis_index/abs(axis_index);

	if(abs(axis_index)==1)
	{
		Lxx = Lx;
		Lyy = Ly;
		Lzz = Lz;
		dxx = X2.x - X1.x;
		dyy = X2.y - X1.y;
		dzz = X2.z - X1.z;
	}

	if(abs(axis_index)==2)
	{
		Lxx = Ly;
		Lyy = Lz;
		Lzz = Lx;
		dxx = X2.y - X1.y;
		dyy = X2.z - X1.z;
		dzz = X2.x - X1.x;
	}

	if(abs(axis_index)==3)
	{
		Lxx = Lz;
		Lyy = Lx;
		Lzz = Ly;
		dxx = X2.z - X1.z;
		dyy = X2.x - X1.x;
		dzz = X2.y - X1.y;
	}

	if(dxx < -Lxx/2){dxx += Lxx;}
	if(dxx > +Lxx/2){dxx -= Lxx;}
	if(dyy < -Lyy/2){dyy += Lyy;}
	if(dyy > +Lyy/2){dyy -= Lyy;}
	if(dzz < -Lzz/2){dzz += Lzz;}
	if(dzz > +Lzz/2){dzz -= Lzz;}

	if(dyy < 1E-17 && dzz < 1E-17)
	{
		dyy += 1E-17;
	}

	double drr2 = dyy*dyy + dzz*dzz;
	double drr = sqrt(dxx*dxx + drr2);
	dxx *= sign;
	double xfrmin = sqrt(rmin*rmin - drr2);
	double xfrcut = sqrt(rcut*rcut - drr2);

	double q = -log(1 - Uniform_Random());
	double s1, s2;
	double dr_1_sqr, dr_2_sqr;
	double q1, q2;

	if(drr2 < rmin*rmin)
	{
		if(dxx <= 0)
		{
			if(drr < rmin)
			{
				s1 = dxx + xfrmin;
				s2 = dxx + xfrcut;
				dr_1_sqr = distance_sqr_Calc(s1, drr2, dxx, Lxx);
				dr_2_sqr = distance_sqr_Calc(s2, drr2, dxx, Lxx);
				q1 = 0;
				q2 = VLJCalc(epsilon, sigma, dr_2_sqr) - VLJCalc(epsilon, sigma, dr_1_sqr) + q1;

				if(q2 < q)
				{
					s1 = Lxx + dxx - xfrmin;
					s2 = Lxx + dxx;
					dr_1_sqr = distance_sqr_Calc(s1, drr2, dxx, Lxx);
					dr_2_sqr = distance_sqr_Calc(s2, drr2, dxx, Lxx);
					q1 = q2;
					q2 = VLJCalc(epsilon, sigma, dr_2_sqr) - VLJCalc(epsilon, sigma, dr_1_sqr) + q1;
				}
			}

			else if(drr >= rmin && drr < rcut)
			{
				s1 = 0;
				s2 = dxx + xfrcut;
				dr_1_sqr = distance_sqr_Calc(s1, drr2, dxx, Lxx);
				dr_2_sqr = distance_sqr_Calc(s2, drr2, dxx, Lxx);
				q1 = 0;
				q2 = VLJCalc(epsilon, sigma, dr_2_sqr) - VLJCalc(epsilon, sigma, dr_1_sqr) + q1;

				if(q2 < q)
				{
					s1 = Lxx + dxx - xfrmin;
					s2 = Lxx + dxx;
					dr_1_sqr = distance_sqr_Calc(s1, drr2, dxx, Lxx);
					dr_2_sqr = distance_sqr_Calc(s2, drr2, dxx, Lxx);
					q1 = q2;
					q2 = VLJCalc(epsilon, sigma, dr_2_sqr) - VLJCalc(epsilon, sigma, dr_1_sqr) + q1;

					if(q2 < q)
					{
						s1 = Lxx + dxx + xfrmin;
						s2 = Lxx;
						dr_1_sqr = distance_sqr_Calc(s1, drr2, dxx, Lxx);
						dr_2_sqr = distance_sqr_Calc(s2, drr2, dxx, Lxx);
						q1 = q2;
						q2 = VLJCalc(epsilon, sigma, dr_2_sqr) - VLJCalc(epsilon, sigma, dr_1_sqr) + q1;
					}
				}
			}

			else
			{
				s1 = Lxx + dxx - xfrmin;
				s2 = Lxx + dxx;
				dr_1_sqr = distance_sqr_Calc(s1, drr2, dxx, Lxx);
				dr_2_sqr = distance_sqr_Calc(s2, drr2, dxx, Lxx);
				q1 = 0;
				q2 = VLJCalc(epsilon, sigma, dr_2_sqr) - VLJCalc(epsilon, sigma, dr_1_sqr) + q1;

				if(q2 < q)
				{
					s1 = Lxx + dxx + xfrmin;
					s2 = Lxx + dxx + xfrcut;
					dr_1_sqr = distance_sqr_Calc(s1, drr2, dxx, Lxx);
					dr_2_sqr = distance_sqr_Calc(s2, drr2, dxx, Lxx);
					q1 = q2;
					q2 = VLJCalc(epsilon, sigma, dr_2_sqr) - VLJCalc(epsilon, sigma, dr_1_sqr) + q1;
				}
			}
		}

		else
		{
			if(drr < rmin)
			{
				s1 = 0;
				s2 = dxx;
				dr_1_sqr = distance_sqr_Calc(s1, drr2, dxx, Lxx);
				dr_2_sqr = distance_sqr_Calc(s2, drr2, dxx, Lxx);
				q1 = 0;
				q2 = VLJCalc(epsilon, sigma, dr_2_sqr) - VLJCalc(epsilon, sigma, dr_1_sqr) + q1;

				if(q2 < q)
				{
					s1 = dxx + xfrmin;
					s2 = dxx + xfrcut;
					dr_1_sqr = distance_sqr_Calc(s1, drr2, dxx, Lxx);
					dr_2_sqr = distance_sqr_Calc(s2, drr2, dxx, Lxx);
					q1 = q2;
					q2 = VLJCalc(epsilon, sigma, dr_2_sqr) - VLJCalc(epsilon, sigma, dr_1_sqr) + q1;

					if(q2 < q)
					{
						s1 = Lxx + dxx -xfrmin;
						s2 = Lxx;
						dr_1_sqr = distance_sqr_Calc(s1, drr2, dxx, Lxx);
						dr_2_sqr = distance_sqr_Calc(s2, drr2, dxx, Lxx);
						q1 = q2;
						q2 = VLJCalc(epsilon, sigma, dr_2_sqr) - VLJCalc(epsilon, sigma, dr_1_sqr) + q1;
					}
				}
			}

			else
			{
				s1 = dxx - xfrmin;
				s2 = dxx;
				dr_1_sqr = distance_sqr_Calc(s1, drr2, dxx, Lxx);
				dr_2_sqr = distance_sqr_Calc(s2, drr2, dxx, Lxx);
				q1 = 0;
				q2 = VLJCalc(epsilon, sigma, dr_2_sqr) - VLJCalc(epsilon, sigma, dr_1_sqr) + q1;

				if(q2 < q)
				{
					s1 = dxx + xfrmin;
					s2 = dxx + xfrcut;
					dr_1_sqr = distance_sqr_Calc(s1, drr2, dxx, Lxx);
					dr_2_sqr = distance_sqr_Calc(s2, drr2, dxx, Lxx);
					q1 = q2;
					q2 = VLJCalc(epsilon, sigma, dr_2_sqr) - VLJCalc(epsilon, sigma, dr_1_sqr) + q1;
				}
			}
		}
	}

	else if((drr2 >= rmin*rmin) && (drr2 < rcut*rcut))
	{
		if(dxx <= 0)
		{
			if(drr < rcut)
			{
				s1 = 0;
				s2 = dxx + xfrcut;
				dr_1_sqr = distance_sqr_Calc(s1, drr2, dxx, Lxx);
				dr_2_sqr = distance_sqr_Calc(s2, drr2, dxx, Lxx);
				q1 = 0;
				q2 = VLJCalc(epsilon, sigma, dr_2_sqr) - VLJCalc(epsilon, sigma, dr_1_sqr) + q1;

				if(q2 < q)
				{
					s1 = Lxx + dxx;
					s2 = Lxx;
					dr_1_sqr = distance_sqr_Calc(s1, drr2, dxx, Lxx);
					dr_2_sqr = distance_sqr_Calc(s2, drr2, dxx, Lxx);
					q1 = q2;
					q2 = VLJCalc(epsilon, sigma, dr_2_sqr) - VLJCalc(epsilon, sigma, dr_1_sqr) + q1;
				}
			}

			else
			{
				s1 = Lxx + dxx;
				s2 = Lxx + dxx + xfrcut;
				dr_1_sqr = distance_sqr_Calc(s1, drr2, dxx, Lxx);
				dr_2_sqr = distance_sqr_Calc(s2, drr2, dxx, Lxx);
				q1 = 0;
				q2 = VLJCalc(epsilon, sigma, dr_2_sqr) - VLJCalc(epsilon, sigma, dr_1_sqr) + q1;
			}
		}

		else
		{
			s1 = dxx;
			s2 = dxx + xfrcut;
			dr_1_sqr = distance_sqr_Calc(s1, drr2, dxx, Lxx);
			dr_2_sqr = distance_sqr_Calc(s2, drr2, dxx, Lxx);
			q1 = 0;
			q2 = VLJCalc(epsilon, sigma, dr_2_sqr) - VLJCalc(epsilon, sigma, dr_1_sqr) + q1;
		}
	}

	else
	{
		return 2*Lxx;
	}

	if(q2<q){return 2*Lxx;}


	//Find solution directly
	double deltaV = q - q1 +VLJCalc(epsilon, sigma, dr_1_sqr);
	double R0_sqr,R1_sqr,R_sqr;
	double temp1=sqrt(1+deltaV/epsilon);
	R0_sqr = pow(2.0/(1+temp1),1.0/3.0)*sigma*sigma;
	R1_sqr = pow(2.0/(1-temp1),1.0/3.0)*sigma*sigma;

	if((R0_sqr-dr_1_sqr)*(R0_sqr-dr_2_sqr)<=0)R_sqr=R0_sqr;
	if((R1_sqr-dr_1_sqr)*(R1_sqr-dr_2_sqr)<=0)R_sqr=R1_sqr;

	double S;
	double temp2;
	temp2=sqrt(R_sqr-drr2);
	S = dxx - temp2;
	if((S>=s1) && (S<=s2))return S;

	S = dxx + temp2;
	if((S>=s1) && (S<=s2))return S;

	S = Lx + dxx - temp2;
	if((S>=s1) && (S<=s2))return S;

	S = Lx + dxx + temp2;
	if((S>=s1) && (S<=s2))return S;

	return 2*Lxx;
}

void Create_LJ_Interaction_Between_Types(int Type_id1, int Type_id2, bool Using_CellList1, bool Using_CellList2, double sigma, double epsilon, double rcut)
{
	if((sigma > Lx/3) || (sigma > Ly/3) || (sigma > Lz/3))
	{
		cout << "Lj Potential Error: interaction distance should be smaller than 1/3 of system size." << endl;
		exit(1);
	}

	if(max(Type_id1, Type_id2) >= Types.size())
	{
		cout << "LJ Potential Error: at least one of the types does not exist." << endl;
		exit(1);
	}

	Parameter_List Params(Lx, Ly, Lz);
	Params.data[3] = sigma;
	Params.data[4] = epsilon;
	Params.data[5] = rcut;
	Parameter_List_For_Short_Range_Interaction.push_back(Params);
	double*data;
	data=Parameter_List_For_Short_Range_Interaction[Parameter_List_For_Short_Range_Interaction.size()-1].data;
	//register interaction
	Short_Range_Interaction_Between_Types*SR;
	SR=new Short_Range_Interaction_Between_Types(Type_id1,Type_id2,&(Types[Type_id1].X),&(Types[Type_id2].X),Event_Time_LJ_Potential,data,rcut,Using_CellList1,Using_CellList2);
	Short_Range_Interaction_Between_Types_List.push_back(SR);
	int Interaction_Global_ID = Short_Range_Interaction_Between_Types_List.size()-1;

    //connect two types
    Types[Type_id1].Interactions_with_Types.push_back(Interaction_Global_ID);
    if(Type_id1==Type_id2)return;
    Types[Type_id2].Interactions_with_Types.push_back(Interaction_Global_ID);
    return;
}



inline double VGaussCalc(double epsilon, double sigma, double r_sqr)
{
	return epsilon*exp(-0.5*r_sqr/(sigma*sigma));
}

double Event_Time_Gauss_Potential(double4 X1, double4 X2, int axis_index, double *Params, double Max_Event_T)
{
    //Move X1 along axis_index direction, get reject event time of spring interaction with X2, with parameters Params.
    //If the event time is less than Max_Event_T and less than period in the moving direction, return real event time.
    //If the event time is greater than Max_Event_T, return 2*(period of moving direction), to reduce the time cost. However, in this function, we just return real value.
    //If the event time is greater than one period of moving direction, return 2*(period of moving direction)

	double Lx, Ly, Lz, sigma, epsilon, rcut;

	Lx = Params[0];
	Ly = Params[1];
	Lz = Params[2];
	sigma = Params[3];
	epsilon = Params[4];
	rcut = Params[5];

	//double shift = epsilon*exp(-0.5*rcut*rcut/sigma/sigma);

	double dxx, dyy, dzz;
	double Lxx, Lyy, Lzz;

	int sign = axis_index/abs(axis_index);

	if(abs(axis_index)==1)
	{
		Lxx = Lx;
		Lyy = Ly;
		Lzz = Lz;
		dxx = X2.x - X1.x;
		dyy = X2.y - X1.y;
		dzz = X2.z - X1.z;
	}

	if(abs(axis_index)==2)
	{
		Lxx = Ly;
		Lyy = Lz;
		Lzz = Lx;
		dxx = X2.y - X1.y;
		dyy = X2.z - X1.z;
		dzz = X2.x - X1.x;
	}

	if(abs(axis_index)==3)
	{
		Lxx = Lz;
		Lyy = Lx;
		Lzz = Ly;
		dxx = X2.z - X1.z;
		dyy = X2.x - X1.x;
		dzz = X2.y - X1.y;
	}

	if(dxx < -Lxx/2){dxx += Lxx;}
	if(dxx > +Lxx/2){dxx -= Lxx;}
	if(dyy < -Lyy/2){dyy += Lyy;}
	if(dyy > +Lyy/2){dyy -= Lyy;}
	if(dzz < -Lzz/2){dzz += Lzz;}
	if(dzz > +Lzz/2){dzz -= Lzz;}

	if(dyy < 1E-17 && dzz < 1E-17)
	{
		dyy += 1E-17;
	}

	double drr2 = dyy*dyy + dzz*dzz;
	double drr = dxx*dxx + drr2;
	dxx *= sign;
	double xfrcut = sqrt(rcut*rcut - drr2);

	double q = -log(1 - Uniform_Random());
	double s1, s2;
	double dr_1_sqr, dr_2_sqr;
	double q1, q2;


	if(epsilon > 0)
	{
		if(drr < rcut*rcut)
		{
			if(dxx > 0)
			{
				s1 = 0;
				s2 = dxx;
				dr_1_sqr = distance_sqr_Calc(s1, drr2, dxx, Lxx);
				dr_2_sqr = distance_sqr_Calc(s2, drr2, dxx, Lxx);
				q1 = 0;
				q2 = VGaussCalc(epsilon, sigma, dr_2_sqr) - VGaussCalc(epsilon, sigma, dr_1_sqr) + q1;

				if(q2 < q)
				{
					s1 = Lxx + dxx - xfrcut;
					s2 = Lxx;
					dr_1_sqr = distance_sqr_Calc(s1, drr2, dxx, Lxx);
					dr_2_sqr = distance_sqr_Calc(s2, drr2, dxx, Lxx);
					q1 = q2;
					q2 = VGaussCalc(epsilon, sigma, dr_2_sqr) - VGaussCalc(epsilon, sigma, dr_1_sqr) + q1;
				}
			}

			else
			{
				s1 = Lxx + dxx - xfrcut;
				s2 = Lxx + dxx;
				dr_1_sqr = distance_sqr_Calc(s1, drr2, dxx, Lxx);
				dr_2_sqr = distance_sqr_Calc(s2, drr2, dxx, Lxx);
				q1 = 0;
				q2 = VGaussCalc(epsilon, sigma, dr_2_sqr) - VGaussCalc(epsilon, sigma, dr_1_sqr) + q1;
			}	
		}

		else if(drr2 < rcut*rcut)
		{
			if(dxx > 0)
			{
				s1 = dxx - xfrcut;
				s2 = dxx;
				dr_1_sqr = distance_sqr_Calc(s1, drr2, dxx, Lxx);
				dr_2_sqr = distance_sqr_Calc(s2, drr2, dxx, Lxx);
				q1 = 0;
				q2 = VGaussCalc(epsilon, sigma, dr_2_sqr) - VGaussCalc(epsilon, sigma, dr_1_sqr) + q1;
			}

			else
			{
				s1 = Lxx + dxx - xfrcut;
				s2 = Lxx + dxx;
				dr_1_sqr = distance_sqr_Calc(s1, drr2, dxx, Lxx);
				dr_2_sqr = distance_sqr_Calc(s2, drr2, dxx, Lxx);
				q1 = 0;
				q2 = VGaussCalc(epsilon, sigma, dr_2_sqr) - VGaussCalc(epsilon, sigma, dr_1_sqr) + q1;
			}
		}

		else
		{
			return 2*Lxx;
		}
	}

	if(epsilon < 0)
	{
		if(drr < rcut*rcut)
		{
			if(dxx > 0)
			{
				s1 = dxx;
				s2 = dxx + xfrcut;
				dr_1_sqr = distance_sqr_Calc(s1, drr2, dxx, Lxx);
				dr_2_sqr = distance_sqr_Calc(s2, drr2, dxx, Lxx);
				q1 = 0;
				q2 = VGaussCalc(epsilon, sigma, dr_2_sqr) - VGaussCalc(epsilon, sigma, dr_1_sqr) + q1;
			}

			else
			{
				s1 = 0;
				s2 = dxx + xfrcut;
				dr_1_sqr = distance_sqr_Calc(s1, drr2, dxx, Lxx);
				dr_2_sqr = distance_sqr_Calc(s2, drr2, dxx, Lxx);
				q1 = 0;
				q2 = VGaussCalc(epsilon, sigma, dr_2_sqr) - VGaussCalc(epsilon, sigma, dr_1_sqr) + q1;

				if(q2 < q)
				{
					s1 = Lxx + dxx;
					s2 = Lxx;
					dr_1_sqr = distance_sqr_Calc(s1, drr2, dxx, Lxx);
					dr_2_sqr = distance_sqr_Calc(s2, drr2, dxx, Lxx);
					q1 = q2;
					q2 = VGaussCalc(epsilon, sigma, dr_2_sqr) - VGaussCalc(epsilon, sigma, dr_1_sqr) + q1;
				}
			}
		}

		else if(drr2 < rcut*rcut)
		{
			if(dxx > 0)
			{
				s1 = dxx;
				s2 = dxx + xfrcut;
				dr_1_sqr = distance_sqr_Calc(s1, drr2, dxx, Lxx);
				dr_2_sqr = distance_sqr_Calc(s2, drr2, dxx, Lxx);
				q1 = 0;
				q2 = VGaussCalc(epsilon, sigma, dr_2_sqr) - VGaussCalc(epsilon, sigma, dr_1_sqr) + q1;
			}

			else
			{
				s1 = Lxx + dxx;
				s2 = Lxx + dxx + xfrcut;
				dr_1_sqr = distance_sqr_Calc(s1, drr2, dxx, Lxx);
				dr_2_sqr = distance_sqr_Calc(s2, drr2, dxx, Lxx);
				q1 = 0;
				q2 = VGaussCalc(epsilon, sigma, dr_2_sqr) - VGaussCalc(epsilon, sigma, dr_1_sqr) + q1;
			}
		}

		else
		{
			return 2*Lxx;
		}
	}

	if(q2<q){return 2*Lxx;}

	//Find solution directly
	double deltaV = q - q1 +VGaussCalc(epsilon, sigma, dr_1_sqr);

	double R_sqr;
	R_sqr = -2.0*log(deltaV/epsilon)*sigma*sigma;
	
	double S;
	double temp;
	temp=sqrt(R_sqr-drr2);
	S = dxx - temp;
	if(S>=s1 && S<=s2)return S;

	S = dxx + temp;
	if(S>=s1 && S<=s2)return S;

	S = Lx + dxx - temp;
	if(S>=s1 && S<=s2)return S;

	S = Lx + dxx + temp;
	if(S>=s1 && S<=s2)return S;

	return 2*Lxx;

}

void Create_Gauss_Interaction_Between_Types(int Type_id1, int Type_id2, bool Using_CellList1, bool Using_CellList2, double sigma, double epsilon, double rcut)
{
	if((sigma > Lx/3) || (sigma > Ly/3) || (sigma > Lz/3))
	{
		cout << "Gauss Potential Error: interaction distance should be smaller than 1/3 of system size." << endl;
		exit(1);
	}

	if(max(Type_id1, Type_id2) >= Types.size())
	{
		cout << "Gauss Potential Error: at least one of the types does not exist." << endl;
		exit(1);
	}

	Parameter_List Params(Lx, Ly, Lz);
	Params.data[3] = sigma;
	Params.data[4] = epsilon;
	Params.data[5] = rcut;
	Parameter_List_For_Short_Range_Interaction.push_back(Params);
	double*data;
	data=Parameter_List_For_Short_Range_Interaction[Parameter_List_For_Short_Range_Interaction.size()-1].data;
	//register interaction
	Short_Range_Interaction_Between_Types*SR;
	SR=new Short_Range_Interaction_Between_Types(Type_id1,Type_id2,&(Types[Type_id1].X),&(Types[Type_id2].X),Event_Time_Gauss_Potential,data,rcut,Using_CellList1,Using_CellList2);
	Short_Range_Interaction_Between_Types_List.push_back(SR);
	int Interaction_Global_ID = Short_Range_Interaction_Between_Types_List.size()-1;

	//connect two types
	Types[Type_id1].Interactions_with_Types.push_back(Interaction_Global_ID);
	if(Type_id1 == Type_id2)return;
	Types[Type_id2].Interactions_with_Types.push_back(Interaction_Global_ID);
	return;
}
