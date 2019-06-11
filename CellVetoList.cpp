/*
Author: Yihao Liang
liangyihaosjtu@gmail.com
This code is for Event Chain Monte Carlo for pairwise interacting many body system
*/
//Cell-Veto List
//deal with long-range interactions
//range:[0,Lx],[0,Ly],[0,Lz]
#include <cmath>
#include "CellVetoList.hpp"
#include "Get_Max_q.hpp"
#include "public.hpp"
#include "omp.h"

void CellVetoList::Insert(int2 ids, int3 CellID) {
	Veto_Cells[CellID.x][CellID.y][CellID.z].particle_list.push_back(ids);
	if(Veto_Cells[CellID.x][CellID.y][CellID.z].particle_list.size()>Num_Particle_Per_Cell){
		Exception_Particle_List.push_back(ids);
	}
}
void CellVetoList::Delete(int2 ids, int3 CellID) {
		Veto_Cell*VC;
		VC=&Veto_Cells[CellID.x][CellID.y][CellID.z];
		int Num_Particle_In_This_Cell;
		int Num_Particle_Exception;
		Num_Particle_In_This_Cell=VC->particle_list.size();
		Num_Particle_Exception=Exception_Particle_List.size();

		for(int k=0;k<Num_Particle_In_This_Cell;k++)
			if((VC->particle_list[k].x==ids.x)&&(VC->particle_list[k].y==ids.y)) {
				
				if(k>=Num_Particle_Per_Cell){//remove ids from exception particle list
					for(int l=0;l<Num_Particle_Exception;l++){
						if((Exception_Particle_List[l].x==ids.x)&&(Exception_Particle_List[l].y==ids.y)) {
							Exception_Particle_List[l]=Exception_Particle_List[Num_Particle_Exception-1];
							Exception_Particle_List.pop_back();
							break;
						}
					}
				}


				if((k<Num_Particle_Per_Cell)&&(Num_Particle_In_This_Cell>Num_Particle_Per_Cell)){//remove last bead from exception particle list
					for(int l=0;l<Num_Particle_Exception;l++){
						if((Exception_Particle_List[l].x==VC->particle_list[Num_Particle_In_This_Cell-1].x)&&(Exception_Particle_List[l].y==VC->particle_list[Num_Particle_In_This_Cell-1].y)) {
							Exception_Particle_List[l]=Exception_Particle_List[Num_Particle_Exception-1];
							Exception_Particle_List.pop_back();
							break;
						}
					}
				}

				VC->particle_list[k]=VC->particle_list[Num_Particle_In_This_Cell-1];
				VC->particle_list.pop_back();
				return;
			}
}


CellVetoList::CellVetoList(double Lx, double Ly, double Lz, double valence, vector<Bead_Type>*Types_pointer) {
	this->Types_pointer=Types_pointer;
	NParticles=0;
	vector<int>type_id_list;
	type_id_list.clear();
	ES=new Ewald_Sum(Lx,Ly,Lz,1E-10);
	//All the charged particles with assigned valence should be listed in thise Cell-Veto-List, now check
	InWhichVetoCell.resize(Types_pointer->size());
    for(int type_id=0;type_id<Types_pointer->size();type_id++){
		InWhichVetoCell[type_id].clear();
        if((*Types_pointer)[type_id].X.size()==0)continue;
		if(abs((*Types_pointer)[type_id].X[0].w-valence)>EPSILON)continue;
		type_id_list.push_back(type_id);
		for(int k=0;k<(*Types_pointer)[type_id].X.size();k++){
            
			if(abs((*Types_pointer)[type_id].X[k].w-valence)>EPSILON){
                cout<<"Error, valence doesn't match "<<endl;
				exit(-1);
            }
			NParticles++;
        }
    }
	this->valence=valence;
	//Now generate grid
	double dL;
	dL=pow(Lx*Ly*Lz/NParticles,1.0/3.0);
	NC_x=(int)(Lx/dL);//for x>NC_x*dL, it is absorbed into last cell
	NC_y=(int)(Ly/dL);//for y>NC_y*dL, it is absorbed into last cell
	NC_z=(int)(Lz/dL);//for z>NC_z*dL, it is absorbed into last cell

	if(NC_x<5)NC_x=5;
	if(NC_y<5)NC_y=5;
	if(NC_z<5)NC_z=5;

	this->Lx=Lx;
	this->Ly=Ly;
	this->Lz=Lz;

	this->dLx=Lx/NC_x;
    this->dLy=Ly/NC_y;
    this->dLz=Lz/NC_z;

	this->Num_Particle_Per_Cell=1;

	//Now construct Cell List
	//vector<vector<int3> >InWhichCell;
	//InWhichCell[type_id][bead_id] stores the in which cell each bead is
	//vector<vector<vector<Cell> > >Cells;
	//If there are more than Num_Particle_Per_Cell particles in a cell, the excessive particle should be put into Exception_Particle_List
	//Veto_Cell Empty_Cell;
	Veto_Cells.resize(NC_x);
	for(int k=0;k<NC_x;k++){
		Veto_Cells[k].resize(NC_y);
		for(int l=0;l<NC_y;l++) {
			Veto_Cells[k][l].resize(NC_z);
		}
	}
	Exception_Particle_List.clear();

	int2 ids;
	int3 IWC;
	double4 X;
	for(int k=0;k<type_id_list.size();k++){
		int type_id=type_id_list[k];
		for(int bead_id=0;bead_id<(*Types_pointer)[type_id].X.size();bead_id++) {
			ids.x=type_id;
			ids.y=bead_id;
			X=(*Types_pointer)[type_id].X[bead_id];
			IWC=In_Which_Veto_Cell(X);
			InWhichVetoCell[type_id].push_back(IWC); 
			Insert(ids, IWC);
		}
	}

	//Now construct max event rate
	qx_max.clear();
	qy_max.clear();
	qz_max.clear();
	Qx_tot=0;
	Qy_tot=0;
	Qz_tot=0;
	vector<double>qx_max_hist(NC_x*NC_y*NC_z);
	vector<double>qy_max_hist(NC_x*NC_y*NC_z);
	vector<double>qz_max_hist(NC_x*NC_y*NC_z);
	//For x-axis
	qx_max.resize(NC_x);
	for(int Ix=0;Ix<NC_x;Ix++){
		qx_max[Ix].resize(NC_y);
		for(int Iy=0;Iy<NC_y;Iy++){
			qx_max[Ix][Iy].resize(NC_z);
		}
	}
	int Ix;
	double Qx_tot_temp=0;
	#pragma omp parallel for reduction (+:Qx_tot_temp)
	for(Ix=0;Ix<NC_x;Ix++){
		for(int Iy=0;Iy<NC_y;Iy++)
			for(int Iz=0;Iz<NC_z;Iz++){

				int Ixr=Ix,Iyr=Iy,Izr=Iz;
				if(Ixr>NC_x/2)Ixr-=NC_x;
				if(Iyr>NC_y/2)Iyr-=NC_y;
				if(Izr>NC_z/2)Izr-=NC_z;
				if((abs(Ixr)<=1)&&(abs(Iyr)<=1)&&(abs(Izr)<=1)){
					qx_max[Ix][Iy][Iz]=0;
					continue;
				}

				double xb1,xb2,yb1,yb2,zb1,zb2;
				xb1=(Ix-1)*dLx;
				xb2=(Ix+1)*dLx;
				yb1=(Iy-1)*dLy;
				yb2=(Iy+1)*dLy;
				zb1=(Iz-1)*dLz;
				zb2=(Iz+1)*dLz;
				qx_max[Ix][Iy][Iz]=get_max_q(xb1, xb2, yb1, yb2, zb1, zb2, *ES, 1);
				Qx_tot_temp+=qx_max[Ix][Iy][Iz];
				int temp;
				int3 I;
				I.x=Ix;I.y=Iy;I.z=Iz;
				temp=Int3_To_Int(I);
				qx_max_hist[temp]=qx_max[Ix][Iy][Iz];
			}
	}
	Qx_tot=Qx_tot_temp;
	cout<<"Qx_tot="<<Qx_tot<<endl;
	FG_x=new Frequency_Generator(qx_max_hist);
	//For y-axis
	qy_max.resize(NC_x);
	for(int Ix=0;Ix<NC_x;Ix++){
		qy_max[Ix].resize(NC_y);
		for(int Iy=0;Iy<NC_y;Iy++){
			qy_max[Ix][Iy].resize(NC_z);
		}
	}

	double Qy_tot_temp=0;
	#pragma omp parallel for reduction(+:Qy_tot_temp)
	for(Ix=0;Ix<NC_x;Ix++){
		for(int Iy=0;Iy<NC_y;Iy++)
			for(int Iz=0;Iz<NC_z;Iz++){


				int Ixr=Ix,Iyr=Iy,Izr=Iz;
				if(Ixr>NC_x/2)Ixr-=NC_x;
				if(Iyr>NC_y/2)Iyr-=NC_y;
				if(Izr>NC_z/2)Izr-=NC_z;
				if((abs(Ixr)<=1)&&(abs(Iyr)<=1)&&(abs(Izr)<=1)){
					qy_max[Ix][Iy][Iz]=0;
					continue;
				}

				double xb1,xb2,yb1,yb2,zb1,zb2;
				xb1=(Ix-1)*dLx;
				xb2=(Ix+1)*dLx;
				yb1=(Iy-1)*dLy;
				yb2=(Iy+1)*dLy;
				zb1=(Iz-1)*dLz;
				zb2=(Iz+1)*dLz;
				qy_max[Ix][Iy][Iz]=get_max_q(xb1, xb2, yb1, yb2, zb1, zb2, *ES, 2);
				Qy_tot_temp+=qy_max[Ix][Iy][Iz];
				int temp;
				int3 I;
				I.x=Ix;I.y=Iy;I.z=Iz;
				temp=Int3_To_Int(I);
				qy_max_hist[temp]=qy_max[Ix][Iy][Iz];
			}
	}
	Qy_tot=Qy_tot_temp;
	cout<<"Qy_tot="<<Qy_tot<<endl;

	FG_y=new Frequency_Generator(qy_max_hist);
	//For z-axis
	qz_max.resize(NC_x);
	for(int Ix=0;Ix<NC_x;Ix++){
		qz_max[Ix].resize(NC_y);
		for(int Iy=0;Iy<NC_y;Iy++){
			qz_max[Ix][Iy].resize(NC_z);
		}
	}
	double Qz_tot_temp=0;
	#pragma omp parallel for reduction(+:Qz_tot_temp)
	for(Ix=0;Ix<NC_x;Ix++){
		for(int Iy=0;Iy<NC_y;Iy++)
			for(int Iz=0;Iz<NC_z;Iz++){

				int Ixr=Ix,Iyr=Iy,Izr=Iz;
				if(Ixr>NC_x/2)Ixr-=NC_x;
				if(Iyr>NC_y/2)Iyr-=NC_y;
				if(Izr>NC_z/2)Izr-=NC_z;
				if((abs(Ixr)<=1)&&(abs(Iyr)<=1)&&(abs(Izr)<=1)){
					qz_max[Ix][Iy][Iz]=0;
					continue;
				}

				double xb1,xb2,yb1,yb2,zb1,zb2;
				xb1=(Ix-1)*dLx;
				xb2=(Ix+1)*dLx;
				yb1=(Iy-1)*dLy;
				yb2=(Iy+1)*dLy;
				zb1=(Iz-1)*dLz;
				zb2=(Iz+1)*dLz;
				qz_max[Ix][Iy][Iz]=get_max_q(xb1, xb2, yb1, yb2, zb1, zb2, *ES, 3);
				Qz_tot_temp+=qz_max[Ix][Iy][Iz];
				int temp;
				int3 I;
				I.x=Ix;I.y=Iy;I.z=Iz;
				temp=Int3_To_Int(I);
				qz_max_hist[temp]=qz_max[Ix][Iy][Iz];
			}
	}
	Qz_tot=Qz_tot_temp;
	cout<<"Qz_tot="<<Qz_tot<<endl;
	FG_z=new Frequency_Generator(qz_max_hist);
	#ifdef DEBUG
		calcu_min_q();
	#endif
}


void CellVetoList::Update(int2 const ids, double4 const NX) {
	int type_id,bead_id;
	type_id=ids.x;
	bead_id=ids.y;

	//check if it is in this cell-veto list
	if(InWhichVetoCell[ids.x].size()==0)return;/*{
		cout<<"Error, bead "<<type_id<<","<<bead_id<<" is not in this cell-veto list"<<endl;
		exit(0);
	}*/
	//move it
	int3 IWC1,IWC2;

	IWC1=InWhichVetoCell[ids.x][ids.y];
	IWC2=In_Which_Veto_Cell(NX);

	if((IWC1.x==IWC2.x)&&(IWC1.y==IWC2.y)&&(IWC1.z==IWC2.z))return;

	Delete(ids, IWC1);
	Insert(ids, IWC2);
	InWhichVetoCell[ids.x][ids.y]=IWC2;
}


void CellVetoList::update_max_num_per_cell(){
	vector<int>Num_Exception;
	int num;
	int max_num=-1;
	double S=0;
	int tot_num=0;
	for(int i=0;i<NC_x;i++)
		for(int j=0;j<NC_y;j++)
			for(int k=0;k<NC_z;k++){
				num=Veto_Cells[i][j][k].particle_list.size();
				if(num>max_num)max_num=num;
				while(num>=Num_Exception.size())Num_Exception.push_back(0);
				for(int l=0;l<=num;l++)Num_Exception[l]+=(num-l);
				
				if(num!=0)S+=-num*log(num);
				tot_num+=num;
			}
	S=S/tot_num;//+log(tot_num);
	cout<<"Homogeneity:"<<S<<endl;
	int Old_Num_Particle_Per_Cell;
	Old_Num_Particle_Per_Cell=Num_Particle_Per_Cell;
	if(Num_Exception[max_num-1]<2)
		Num_Particle_Per_Cell=max_num-1;
	else
		Num_Particle_Per_Cell=max_num;

	//if(Num_Particle_Per_Cell==0)Num_Particle_Per_Cell=1;//For debug
	cout<<"Change max number of particle per cell to "<<Num_Particle_Per_Cell<<endl;

	if(Old_Num_Particle_Per_Cell<Num_Particle_Per_Cell){
		int2 ids;
		int3 IWC;
		int m;
		vector<int2>*particle_list_p;
		for(int l=0;l<Exception_Particle_List.size();l++){
			ids=Exception_Particle_List[l];
			IWC=InWhichVetoCell[ids.x][ids.y];
			particle_list_p=&Veto_Cells[IWC.x][IWC.y][IWC.z].particle_list;
			for(m=0;m<particle_list_p->size();m++)if(((*particle_list_p)[m].x==ids.x)&&((*particle_list_p)[m].y==ids.y))break;
			if(m<Num_Particle_Per_Cell){
				//remove from Exception_Particle_List
				Exception_Particle_List[l]=Exception_Particle_List[Exception_Particle_List.size()-1];
				Exception_Particle_List.pop_back();
				l--;
			}
		}
	}else if(Old_Num_Particle_Per_Cell>Num_Particle_Per_Cell){
		int2 ids;
		vector<int2>*particle_list_p;

		Exception_Particle_List.clear();
		for(int i=0;i<NC_x;i++)
			for(int j=0;j<NC_y;j++)
				for(int k=0;k<NC_z;k++){
					particle_list_p=&Veto_Cells[i][j][k].particle_list;
					if(particle_list_p->size()<=Num_Particle_Per_Cell)continue;
					for(int l=Num_Particle_Per_Cell;l<particle_list_p->size();l++){
						ids=(*particle_list_p)[l];
						Exception_Particle_List.push_back(ids);
					}
				}
	}

}
