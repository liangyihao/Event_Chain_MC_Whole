//deal with short-range interactions
//range:[0,Lx],[0,Ly],[0,Lz]
#ifndef CELL_LIST
#define CELL_LIST

#include <vector>
#include <iostream>
#include "public.hpp"
using namespace std;

typedef struct CellStruct{
	vector<vector<int> >particle_list;
	/*
		secondary id of particles in this cell
		particle_list[type_id] stores all the particles of type type_id in current cell
	*/
	CellStruct(int type_num){
		particle_list.resize(type_num);
		for(int k=0;k<type_num;k++)particle_list[k].clear();
	}

	void Insert(int2 ids);
	void Delete(int2 ids);
}Cell;

class CellList{
private:
	vector<Bead_Type>* Types_pointer;
	double Lx,Ly,Lz;
	double dL;//size of each cell
	int NC_x,NC_y,NC_z;//Number of cells per dimension
	vector<vector<int3> >InWhichCell;//InWhichCell[type_id][bead_id] stores the in which cell each bead is
	vector<vector<vector<Cell> > >Cells;
	/*
	CellList[i][j][k] access the cell in 
	i*Lx<= x < (i+1)*Lx
	j*Ly<= y < (j+1)*Ly
	k*Lz<= z < (k+1)*Lz
	*/
public:
	CellList(double Lx, double Ly, double Lz, double dL, vector<Bead_Type>*Types_pointer);
	void Move(int2 ids, double s, int axis);
	int3 In_Which_Cell(int2 ids);
	vector<int>* Particle_List_In_Cell(int3 Cell_IJK, int type);
	int NCell_x(){return NC_x;}
	int NCell_y(){return NC_y;}
	int NCell_z(){return NC_z;}
	double Cell_size(){return dL;}
	void print();
};
#endif
