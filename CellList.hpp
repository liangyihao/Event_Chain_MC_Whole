//deal with short-range interactions
//range:[0,Lx],[0,Ly],[0,Lz]
#ifndef CELL_LIST
#define CELL_LIST

#include <vector>
#include <iostream>
#include "public.hpp"
using namespace std;

typedef struct CellStruct{
	int type_id;//type id of particles in this cell
	vector<int>particle_list;//secondary id of particles in this cell
	CellStruct(int type_id){
		this->type_id=type_id;
		particle_list.clear();
	}

	void Insert(int bead_id) {
		particle_list.push_back(bead_id);
		return;
}
	void Delete(int bead_id) {
		for(int k=0;k<particle_list.size();k++)
			if(particle_list[k]==bead_id) {
				particle_list[k]=particle_list[particle_list.size()-1];
				particle_list.pop_back();
				return;
			}
		cout<<"Deletion failed. The particle ("<<type_id<<","<<bead_id<<") not found."<<endl;
		return;
}

}Cell;

class CellList{
private:
//One cell list only store 1 type of bead and serve only one short range interaction instance
	int type_id;
	double (*Gen)(double4,double4,int,double*,double);
	double *Params;

	vector<double4>*X_pointer;
	double Lx,Ly,Lz;
	double dL;//size of each cell
	int NC_x,NC_y,NC_z;//Number of cells per dimension
	vector<int3>InWhichCell;//InWhichCell[bead_id] stores the in which cell the bead (type_id,bead_is) is
	vector<vector<vector<Cell> > >Cells;
	/*
	CellList[i][j][k] access the cell in 
	i*Lx<= x < (i+1)*Lx
	j*Ly<= y < (j+1)*Ly
	k*Lz<= z < (k+1)*Lz
	*/
	void Event_with_Cell(TwoBody_Event&Event, int3 Cell_IJK, int2 Active_Bead, double4 X_Active_Bead, int axis);
	int3 In_Which_Cell(int2 ids,double4 X_ids) {
		if(ids.x==type_id){
			return InWhichCell[ids.y];
		}else{
			int3 temp;
			temp.x=X_ids.x/dL;
			temp.y=X_ids.y/dL;
			temp.z=X_ids.z/dL;
			
			if(temp.x==NC_x)temp.x=NC_x-1;
			if(temp.y==NC_y)temp.y=NC_y-1;
			if(temp.z==NC_z)temp.z=NC_z-1;

			return temp;
		}
	}

public:
	CellList(double Lx, double Ly, double Lz, double dL, int type_id, vector<double4>*X_pointer, double (*Gen)(double4,double4,int,double*,double),double *Params);
	void Get_Event(TwoBody_Event&Event, int2 Active_Bead, double4 X_Active_Bead, int axis);
	void Update(int2 ids, double4 New_X);
	void print();
};
#endif
