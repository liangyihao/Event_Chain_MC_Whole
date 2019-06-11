/*
Author: Yihao Liang
liangyihaosjtu@gmail.com
This code is for checking
*/
#include <cmath>
#include "CellVetoList.hpp"
#include "Get_Max_q.hpp"
#include "public.hpp"
using namespace std;
void CellVetoList::check_print() {
	//Step1: check cell and iwc
	cout<<"Cell Veto List checking"<<endl;
	int2 ids;
	int3 IWC_mem,IWC_dyn;
	double4 X;
	for(int i=0;i<NC_x;i++)
		for(int j=0;j<NC_y;j++)
			for(int k=0;k<NC_z;k++)
				for(int l=0;l<Veto_Cells[i][j][k].particle_list.size();l++){
					ids=Veto_Cells[i][j][k].particle_list[l];
					X=(*Types_pointer)[ids.x].X[ids.y];
					IWC_mem=InWhichVetoCell[ids.x][ids.y];
					IWC_dyn=In_Which_Veto_Cell(X);
					if((IWC_mem.x!=i)||(IWC_mem.y!=j)||(IWC_mem.z!=k)){
						cout<<"Error Cell information inconsistent"<<endl;
						exit(0);
					}
					if((IWC_dyn.x!=i)||(IWC_dyn.y!=j)||(IWC_dyn.z!=k)){
						cout<<"Warning Cell information inconsistent"<<endl;
						cout<<"X=("<<X.x<<","<<X.y<<","<<X.z<<","<<X.w<<")"<<endl;
						cout<<"correct cell position:"<<IWC_dyn.x<<' '<<IWC_dyn.y<<' '<<IWC_dyn.z<<')'<<endl;
						cout<<"current cell position:"<<i<<' '<<j<<' '<<k<<endl<<endl;
					}
					if(l>=Num_Particle_Per_Cell){//check if it is in the exception particle list
						int m;
						for(m=0;m<Exception_Particle_List.size();m++)if((Exception_Particle_List[m].x==ids.x)&&(Exception_Particle_List[m].y==ids.y))break;
						if(m==Exception_Particle_List.size()){
							cout<<"Error: Exception Particle not in the list"<<endl;
						}
					}
				}
	for(int k=0;k<Exception_Particle_List.size();k++){
		ids=Exception_Particle_List[k];
		IWC_mem=InWhichVetoCell[ids.x][ids.y];
		int l;
		for(l=0;Veto_Cells[IWC_mem.x][IWC_mem.y][IWC_mem.z].particle_list.size();l++){
			if((Veto_Cells[IWC_mem.x][IWC_mem.y][IWC_mem.z].particle_list[l].x==ids.x)&&(Veto_Cells[IWC_mem.x][IWC_mem.y][IWC_mem.z].particle_list[l].y==ids.y))break;
		}
		if(l<Num_Particle_Per_Cell){
			cout<<"Error: Exception Particle incorrect"<<endl;
		}
	}
	cout<<"# of Exception Particles: "<<Exception_Particle_List.size()<<endl;
	cout<<"Q:"<<Qx_tot<<' '<<Qy_tot<<' '<<Qz_tot<<endl;
	cout<<"Cell Veto List checking finished"<<endl;
}

void CellVetoList::check_rate() {
	cout<<"checking rate"<<endl;
	/*
	vector<vector<vector<double> > >qx_max_stat;
	vector<vector<vector<double> > >qy_max_stat;
	vector<vector<vector<double> > >qz_max_stat;
	
	qx_max_stat.resize(NC_x);
	qy_max_stat.resize(NC_x);
	qz_max_stat.resize(NC_x);
	for(int Ix=0;Ix<NC_x;Ix++){
		qx_max_stat[Ix].resize(NC_y);
		qy_max_stat[Ix].resize(NC_y);
		qz_max_stat[Ix].resize(NC_y);
		for(int Iy=0;Iy<NC_y;Iy++){
			qx_max_stat[Ix][Iy].resize(NC_z);
			qy_max_stat[Ix][Iy].resize(NC_z);
			qz_max_stat[Ix][Iy].resize(NC_z);
			for(int Iz=0;Iz<NC_z;Iz++){
				qx_max_stat[Ix][Iy][Iz]=0;
				qy_max_stat[Ix][Iy][Iz]=0;
				qz_max_stat[Ix][Iy][Iz]=0;
			}
		}
	}
\
	double4 X,Y;
	int3 IWC_X,IWC_Y,IWC_R;
	int axis;
	double q,q_bound;
	X.w=1;Y.w=1;
	for(int k=0;k<10000000;k++){
		if(k%100000==0)cout<<k/100000<<"%"<<endl;
		X.x=Uniform_Random()*Lx;
		X.y=Uniform_Random()*Ly;
		X.z=Uniform_Random()*Lz;
		
		Y.x=Uniform_Random()*Lx;
		Y.y=Uniform_Random()*Ly;
		Y.z=Uniform_Random()*Lz;

		IWC_X=In_Which_Veto_Cell(X);
		IWC_Y=In_Which_Veto_Cell(Y);

		IWC_R.x=IWC_Y.x-IWC_X.x;
		IWC_R.y=IWC_Y.y-IWC_X.y;
		IWC_R.z=IWC_Y.z-IWC_X.z;

		if(IWC_R.x<-NC_x/2)IWC_R.x+=NC_x;
		if(IWC_R.x>+NC_x/2)IWC_R.x-=NC_x;

		if(IWC_R.y<-NC_y/2)IWC_R.y+=NC_y;
		if(IWC_R.y>+NC_y/2)IWC_R.y-=NC_y;

		if(IWC_R.z<-NC_z/2)IWC_R.z+=NC_z;
		if(IWC_R.z>+NC_z/2)IWC_R.z-=NC_z;

		if((abs(IWC_R.x)<=1)&&(abs(IWC_R.y)<=1)&&(abs(IWC_R.z)<=1))continue;
		if(IWC_R.x<0)IWC_R.x+=NC_x;
		if(IWC_R.x>=NC_x)IWC_R.x-=NC_x;

		if(IWC_R.y<0)IWC_R.y+=NC_y;
		if(IWC_R.y>=NC_y)IWC_R.y-=NC_y;

		if(IWC_R.z<0)IWC_R.z+=NC_z;
		if(IWC_R.z>=NC_z)IWC_R.z-=NC_z;

		axis=1;
		q=ES->D_Potential(X,Y,axis);
		if(q<0)q=0;
		if(q>qx_max_stat[IWC_R.x][IWC_R.y][IWC_R.z])qx_max_stat[IWC_R.x][IWC_R.y][IWC_R.z]=q;
		q_bound=qx_max[IWC_R.x][IWC_R.y][IWC_R.z];
		if(q>q_bound){
			cout<<q<<' '<<q_bound<<endl;
			exit(0);
		}

		axis=2;
		q=ES->D_Potential(X,Y,axis);
		if(q<0)q=0;
		if(q>qy_max_stat[IWC_R.x][IWC_R.y][IWC_R.z])qy_max_stat[IWC_R.x][IWC_R.y][IWC_R.z]=q;
		q_bound=qy_max[IWC_R.x][IWC_R.y][IWC_R.z];
		if(q>q_bound){
			cout<<q<<' '<<q_bound<<endl;
			exit(0);
		}

		axis=3;
		q=ES->D_Potential(X,Y,axis);
		if(q<0)q=0;
		if(q>qz_max_stat[IWC_R.x][IWC_R.y][IWC_R.z])qz_max_stat[IWC_R.x][IWC_R.y][IWC_R.z]=q;
		q_bound=qz_max[IWC_R.x][IWC_R.y][IWC_R.z];
		if(q>q_bound){
			cout<<q<<' '<<q_bound<<endl;
			exit(0);
		}
	}
	for(int i=0;i<NC_x;i++)
		for(int j=0;j<NC_y;j++)
			for(int k=0;k<NC_z;k++){
				cout<<i<<','<<j<<','<<k<<endl;
				cout<<"qx"<<qx_max_stat[i][j][k]<<' '<<qx_max[i][j][k]<<endl;
				cout<<"qy"<<qy_max_stat[i][j][k]<<' '<<qy_max[i][j][k]<<endl;
				cout<<"qz"<<qz_max_stat[i][j][k]<<' '<<qz_max[i][j][k]<<endl;
				cout<<endl;
			}
	*/
	for(int i=0;i<NC_x;i++)
		for(int j=0;j<NC_y;j++)
			for(int k=0;k<NC_z;k++){
				int I,J,K;
				I=-i;
				J=-j;
				K=-k;
				while(I<0)I+=NC_x;
				while(I>NC_x)I-=NC_x;
				
				while(J<0)J+=NC_y;
				while(J>NC_y)J-=NC_y;
				
				while(K<0)K+=NC_z;
				while(K>NC_z)K-=NC_z;

				//qx_max
				cout<<"At "<<i<<' '<<j<<' '<<k<<endl;
				if(abs(qx_max[i][j][k]-qx_max[i][j][K])>EPSILON)cout<<"qx_max Not symmetric, value:"<<qx_max[i][j][k]<<' '<<qx_max[i][j][K]<<endl;
				if(abs(qx_max[i][j][k]-qx_max[i][J][k])>EPSILON)cout<<"qx_max Not symmetric, value:"<<qx_max[i][j][k]<<' '<<qx_max[i][J][k]<<endl;
				if(abs(qx_max[i][j][k]-qx_max[i][J][K])>EPSILON)cout<<"qx_max Not symmetric, value:"<<qx_max[i][j][k]<<' '<<qx_max[i][J][K]<<endl;

				//qy_max
				if(abs(qy_max[i][j][k]-qy_max[i][j][K])>EPSILON)cout<<"qy_max Not symmetric, value:"<<qy_max[i][j][k]<<' '<<qy_max[i][j][K]<<endl;
				if(abs(qy_max[i][j][k]-qy_max[I][j][k])>EPSILON)cout<<"qy_max Not symmetric, value:"<<qy_max[i][j][k]<<' '<<qy_max[I][j][k]<<endl;
				if(abs(qy_max[i][j][k]-qy_max[I][j][K])>EPSILON)cout<<"qy_max Not symmetric, value:"<<qy_max[i][j][k]<<' '<<qy_max[I][j][K]<<endl;

				//qz_max
				if(abs(qz_max[i][j][k]-qz_max[I][j][k])>EPSILON)cout<<"qz_max Not symmetric, value:"<<qz_max[i][j][k]<<' '<<qz_max[I][j][k]<<endl;
				if(abs(qz_max[i][j][k]-qz_max[i][J][k])>EPSILON)cout<<"qz_max Not symmetric, value:"<<qz_max[i][j][k]<<' '<<qz_max[i][J][k]<<endl;
				if(abs(qz_max[i][j][k]-qz_max[I][J][k])>EPSILON)cout<<"qz_max Not symmetric, value:"<<qz_max[i][j][k]<<' '<<qz_max[I][J][k]<<endl;

			}
	cout<<"checking rate finished"<<endl;
	exit(0);
}

void CellVetoList::calcu_min_q() {

	//Now construct min event rate
	qx_min.clear();
	qy_min.clear();
	qz_min.clear();

	//For x-axis
	qx_min.resize(NC_x);
	qy_min.resize(NC_x);
	qz_min.resize(NC_x);
	for(int Ix=0;Ix<NC_x;Ix++){
		qx_min[Ix].resize(NC_y);
		qy_min[Ix].resize(NC_y);
		qz_min[Ix].resize(NC_y);
		for(int Iy=0;Iy<NC_y;Iy++){
			qx_min[Ix][Iy].resize(NC_z);
			qy_min[Ix][Iy].resize(NC_z);
			qz_min[Ix][Iy].resize(NC_z);
		}
	}
	int Ix;
	#pragma omp parallel for
	for(Ix=0;Ix<NC_x;Ix++){
		for(int Iy=0;Iy<NC_y;Iy++)
			for(int Iz=0;Iz<NC_z;Iz++){

				int Ixr=Ix,Iyr=Iy,Izr=Iz;
				if(Ixr>NC_x/2)Ixr-=NC_x;
				if(Iyr>NC_y/2)Iyr-=NC_y;
				if(Izr>NC_z/2)Izr-=NC_z;
				if((abs(Ixr)<=1)&&(abs(Iyr)<=1)&&(abs(Izr)<=1)){
					qx_min[Ix][Iy][Iz]=0;
					qy_min[Ix][Iy][Iz]=0;
					qz_min[Ix][Iy][Iz]=0;
					continue;
				}

				double xb1,xb2,yb1,yb2,zb1,zb2;
				xb1=(Ix-1)*dLx;
				xb2=(Ix+1)*dLx;
				yb1=(Iy-1)*dLy;
				yb2=(Iy+1)*dLy;
				zb1=(Iz-1)*dLz;
				zb2=(Iz+1)*dLz;
				qx_min[Ix][Iy][Iz]=get_min_q(xb1, xb2, yb1, yb2, zb1, zb2, *ES, 1);
				qy_min[Ix][Iy][Iz]=get_min_q(xb1, xb2, yb1, yb2, zb1, zb2, *ES, 2);
				qz_min[Ix][Iy][Iz]=get_min_q(xb1, xb2, yb1, yb2, zb1, zb2, *ES, 3);
			}
	}
}
