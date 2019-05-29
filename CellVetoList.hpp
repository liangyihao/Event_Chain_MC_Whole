/*
Author: Yihao Liang
liangyihaosjtu@gmail.com
This code is for Event Chain Monte Carlo for pairwise interacting many body system
*/
//Cell-Veto List
//deal with long-range interactions
//range:[0,Lx],[0,Ly],[0,Lz]
//For Each Valance, we need to build a Cell-Veto List seperately
#ifndef CELL_VETO_LIST
#define CELL_VETO_LIST

#include <vector>
#include <iostream>
#include "public.hpp"
#include "Random_Number.hpp"
#include "Ewald_Sum_Factor.hpp"
using namespace std;

typedef struct Veto_Cell_Struct{
	vector<int2> particle_list;//.x: type-id                .y: bead-id
	Veto_Cell_Struct(){
		particle_list.clear();
	}
}Veto_Cell;

class CellVetoList{
private:
	vector<Bead_Type>* Types_pointer;
	double Lx,Ly,Lz;
	double dLx,dLy,dLz;//size of each cell
	int NC_x,NC_y,NC_z;//Number of cells per dimension
	int NParticles;//Number of particles this cell-veto list deals with
	double valence;
	vector<vector<vector<Veto_Cell> > >Veto_Cells;
	/*
	Veto_Cells[i][j][k] access the cell in 
	i*dLx<= x < (i+1)*dLx
	j*dLy<= y < (j+1)*dLy
	k*dLz<= z < (k+1)*dLz
	*/
	vector<vector<int3> >InWhichVetoCell;
	/*
	InWhichCell[type_id][bead_id] stores the in which cell each bead is. 
	If InWhichCell[type_id].size()==0, that means this cell doesn't record this type.
	*/
    int Num_Particle_Per_Cell;
    vector<int2>Exception_Particle_List;

    vector<vector<vector<double> > >qx_max;
	vector<vector<vector<double> > >qy_max;
	vector<vector<vector<double> > >qz_max;
    double Qx_tot,Qy_tot,Qz_tot;
    Frequency_Generator*FG_x;//To assign pre-event to cell
	Frequency_Generator*FG_y;//To assign pre-event to cell
	Frequency_Generator*FG_z;//To assign pre-event to cell
	Ewald_Sum*ES;//Ewald sum calculator
	
	int3 In_Which_Veto_Cell(double4 X){//Compute the cell ids of input coordinate X
		int3 Id;
		Id.x=X.x/dLx;
		Id.y=X.y/dLy;
		Id.z=X.z/dLz;
		if(Id.x==NC_x)Id.x=0;
		if(Id.y==NC_y)Id.y=0;
		if(Id.z==NC_z)Id.z=0;
		return Id;
	}
	int3 Int_To_Int3(int I){
		int3 temp;
		temp.x=I%NC_x;
		I=I/NC_x;
		temp.y=I%NC_y;
		temp.z=I/NC_y;
		return temp;
	}

	int Int3_To_Int(int3 I){
		int temp;
		temp=I.x+I.y*NC_x+I.z*NC_x*NC_y;
		return temp;
	}

	void Insert(int2 ids, int3 CellID);
	void Delete(int2 ids, int3 CellID);

	void Get_Colomb_Event_Exceptional_Particles(int2 id_active_particle, int axis, double Bjerrum_Length, double&time, int2&id_next_active_bead);
	void Get_Colomb_Event_Cell_Veto(int2 id_active_particle, int axis, double Bjerrum_Length, double&time, int2&id_next_active_bead);

public:
    CellVetoList(double Lx, double Ly, double Lz, double valence, vector<Bead_Type>*Types_pointer);
	void Update(int2 const ids, double4 const NX);
	void Get_Colomb_Event(int2 id_active_particle, int axis, double Bjerrum_Length, double&time, int2&id_next_active_bead);
	//void print();
	//void Adjust();//If the number of exception particles is big enough, try to increase the maximum number of particles per cell.
};


#endif
