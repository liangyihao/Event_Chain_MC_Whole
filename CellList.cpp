/*
Author: Yihao Liang
liangyihaosjtu@gmail.com
This code is for Event Chain Monte Carlo for pairwise interacting many body system
*/
#include "CellList.hpp"
#include "public.hpp"

void CellStruct::Insert(int2 ids) {
		particle_list[ids.x].push_back(ids.y);
		return;
}
void CellStruct::Delete(int2 ids) {
		for(int k=0;k<particle_list[ids.x].size();k++)
			if(particle_list[ids.x][k]==ids.y) {
				particle_list[ids.x][k]=particle_list[ids.x][particle_list[ids.x].size()-1];
				particle_list[ids.x].pop_back();
				return;
			}
		cout<<"Deletion failed. The particle ("<<ids.x<<","<<ids.y<<") not found."<<endl;
		return;
}


CellList::CellList(double Lx, double Ly, double Lz, double dL, vector<Bead_Type>*Types_pointer) {
	this->Lx=Lx;
	this->Ly=Ly;
	this->Lz=Lz;
	this->dL=dL;
	this->Types_pointer=Types_pointer;

	NC_x=(int)(Lx/dL);//for x>NC_x*dL, it is absorbed into last cell
	NC_y=(int)(Ly/dL);//for y>NC_y*dL, it is absorbed into last cell
	NC_z=(int)(Lz/dL);//for z>NC_z*dL, it is absorbed into last cell

	//Now construct Cell List
	//vector<vector<int3> >InWhichCell;
	//InWhichCell[type_id][bead_id] stores the in which cell each bead is
	//vector<vector<vector<Cell> > >Cells;
	int type_num=Types_pointer->size();
	
	Cell Empty_Cell(type_num);
	Cells.resize(NC_x);
	for(int k=0;k<NC_x;k++){
		Cells[k].resize(NC_y);
		for(int l=0;l<NC_y;l++) {
			Cells[k][l].clear();
			for(int m=0;m<NC_z;m++) {
				Cells[k][l].push_back(Empty_Cell);
			}
		}
	}

	int2 ids;
	int Ix,Iy,Iz;
	double4 X;
	InWhichCell.resize(type_num);
	for(int type_id=0;type_id<type_num;type_id++){
		InWhichCell[type_id].resize( (*Types_pointer)[type_id].X.size() );
		for(int bead_id=0;bead_id<InWhichCell[type_id].size();bead_id++) {
			X=(*Types_pointer)[type_id].X[bead_id];
			Ix=X.x/dL;
			Iy=X.y/dL;
			Iz=X.z/dL;
			
			if(Ix==NC_x)Ix=NC_x-1;
			if(Iy==NC_y)Iy=NC_y-1;
			if(Iz==NC_z)Iz=NC_z-1;

			InWhichCell[type_id][bead_id].x=Ix;
			InWhichCell[type_id][bead_id].y=Iy;
			InWhichCell[type_id][bead_id].z=Iz;
			ids.x=type_id;
			ids.y=bead_id;
			Cells[Ix][Iy][Iz].Insert(ids);
		}
	}

}


void CellList::Move(int2 ids, double s, int axis) {
	double4 X;
	int3 IWC1,IWC2;

	bool update_tag=false;

	X=(*Types_pointer)[ids.x].X[ids.y];
	IWC1=InWhichCell[ids.x][ids.y];
	IWC2=IWC1;

	if(abs(axis)==1) {
		X.x+=(axis>0?1:-1)*s;
		while(X.x>Lx)X.x-=Lx;
		while(X.x<0 )X.x+=Lx;
		IWC2.x=X.x/dL;
		if(IWC2.x==NC_x)IWC2.x=NC_x-1;
		if(IWC1.x!=IWC2.x)update_tag=true;
	}
	if(abs(axis)==2) {
		X.y+=(axis>0?1:-1)*s;
		while(X.y>Ly)X.y-=Ly;
		while(X.y<0 )X.y+=Ly;
		IWC2.y=X.y/dL;
		if(IWC2.y==NC_y)IWC2.y=NC_y-1;
		if(IWC1.y!=IWC2.y)update_tag=true;
	}
	if(abs(axis)==3) {
		X.z+=(axis>0?1:-1)*s;
		while(X.z>Lz)X.z-=Lz;
		while(X.z<0 )X.z+=Lz;
		IWC2.z=X.z/dL;
		if(IWC2.z==NC_z)IWC2.z=NC_z-1;
		if(IWC1.z!=IWC2.z)update_tag=true;
	}

	(*Types_pointer)[ids.x].X[ids.y]=X;

	if(update_tag){
		Cells[IWC1.x][IWC1.y][IWC1.z].Delete(ids);
		Cells[IWC2.x][IWC2.y][IWC2.z].Insert(ids);
		InWhichCell[ids.x][ids.y]=IWC2;
	}

}

int3 CellList::In_Which_Cell(int2 ids) {
	return InWhichCell[ids.x][ids.y];
}

vector<int>* CellList::Particle_List_In_Cell(int3 Cell_IJK, int type) {
	while(Cell_IJK.x<0)Cell_IJK.x+=NC_x;
	while(Cell_IJK.x>=NC_x)Cell_IJK.x-=NC_x;

	while(Cell_IJK.y<0)Cell_IJK.y+=NC_y;
	while(Cell_IJK.y>=NC_y)Cell_IJK.y-=NC_y;
	
	while(Cell_IJK.z<0)Cell_IJK.z+=NC_z;
	while(Cell_IJK.z>=NC_z)Cell_IJK.z-=NC_z;

	return &(Cells[Cell_IJK.x][Cell_IJK.y][Cell_IJK.z].particle_list[type]);
}

void CellList::print() {
	//Step1: Print Cell List
	
	bool mark;
	for(int i=0;i<NC_x;i++)
		for(int j=0;j<NC_y;j++)
			for(int k=0;k<NC_z;k++) {
				mark=false;
				for(int type_id=0;type_id<Cells[i][j][k].particle_list.size();type_id++)
					for(int uu=0;uu<Cells[i][j][k].particle_list[type_id].size();uu++)
						{
							if(!mark){
								cout<<"In cell ("<<i<<','<<j<<','<<k<<")"<<endl;
								mark=true;
							}
							cout<<"type "<<type_id<<" id: "<<Cells[i][j][k].particle_list[type_id][uu]<<endl;
						}
			}		
	cout<<endl<<endl;

	//Step2: Print "In which cell"
	int2 ids;
	int3 i3;
	for(int type_id=0;type_id<Types_pointer->size();type_id++)
		for(int bead_id=0;bead_id<(*Types_pointer)[type_id].X.size();bead_id++){
			cout<<"Type "<<type_id<<" id: "<<bead_id<<endl;
			cout<<"Position "<<(*Types_pointer)[type_id].X[bead_id].x<<' '<<(*Types_pointer)[type_id].X[bead_id].y<<' '<<(*Types_pointer)[type_id].X[bead_id].z<<endl;
			ids.x=type_id;ids.y=bead_id;
			i3=In_Which_Cell(ids);
			cout<<"In cell "<<i3.x<<' '<<i3.y<<' '<<i3.z<<endl<<endl;
		}
	cout<<endl;
}