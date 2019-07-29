/*
Author: Yihao Liang
liangyihaosjtu@gmail.com
This code is for Event Chain Monte Carlo for pairwise interacting many body system
*/
#include "CellList.hpp"
#include "public.hpp"
#include <cmath>

void CellList::Event_with_Cell(TwoBody_Event&Event, int3 Cell_IJK, int2 Active_Bead, double4 X_Active_Bead, int axis){
    while(Cell_IJK.x<0)Cell_IJK.x+=NC_x;
	while(Cell_IJK.x>=NC_x)Cell_IJK.x-=NC_x;

	while(Cell_IJK.y<0)Cell_IJK.y+=NC_y;
	while(Cell_IJK.y>=NC_y)Cell_IJK.y-=NC_y;
	
	while(Cell_IJK.z<0)Cell_IJK.z+=NC_z;
	while(Cell_IJK.z>=NC_z)Cell_IJK.z-=NC_z;

    double4 Y;
    vector<int>*Particle_List_Pointer;
	Particle_List_Pointer=&(Cells[Cell_IJK.x][Cell_IJK.y][Cell_IJK.z].particle_list);
    int bead_id;
    double t;

    for(int k=0;k<(*Particle_List_Pointer).size();k++) {
        bead_id=(*Particle_List_Pointer)[k];
        if((Active_Bead.x==type_id)&&(Active_Bead.y==bead_id)) continue;
        Y=(*X_pointer)[bead_id];
        t=Gen(X_Active_Bead,Y,axis,Params,Event.Event_Time);
        if(t<Event.Event_Time){
            Event.Event_Time=t;
            Event.Target_Bead.x=type_id;
            Event.Target_Bead.y=bead_id;
        }
    }
}

void CellList::Get_Event(TwoBody_Event&Event, int2 Active_Bead, double4 X_Active_Bead, int axis) {
    //Get Event within time Lx if axis==1(or Ly if axis==2, or Lz if axis==3)
    //If event time out of bound L(=Lx or Ly or Lz),let time=2*L
    //Normally, time is the event time and id_next_active_bead is the id of next active bead

    double4 Y;    
    double (*Gen)(double4,double4,int,double*,double);
    double *Params;

    //Type-Type interactions
    int3 IWC;
    IWC=In_Which_Cell(Active_Bead,X_Active_Bead);
    
    //Run first layer
    int3 Cell_ID;
    for(Cell_ID.x=IWC.x-1;Cell_ID.x<=IWC.x+1;Cell_ID.x++)
        for(Cell_ID.y=IWC.y-1;Cell_ID.y<=IWC.y+1;Cell_ID.y++)
            for(Cell_ID.z=IWC.z-1;Cell_ID.z<=IWC.z+1;Cell_ID.z++){
                Event_with_Cell(Event, Cell_ID, Active_Bead, X_Active_Bead, axis);
            }

    //Run additional layers, if needed
    double Position_1D,Size_1D;
    int N_Layers_1D;//Number of layers on selected dimension

    int3 Cell_In_One_Layer[3][3];
    int3 dC;
    int sign;
    sign=(axis>0)?1:-1;
    if(abs(axis)==1) {
        N_Layers_1D=NC_x;
        Position_1D=X_Active_Bead.x;
        Size_1D=Lx;
        dC.x=sign;
        dC.y=0;
        dC.z=0;
        for(int i=-1;i<=+1;i++)
            for(int j=-1;j<=+1;j++) {
                Cell_In_One_Layer[i+1][j+1].x=IWC.x+sign;
                Cell_In_One_Layer[i+1][j+1].y=IWC.y+i;
                Cell_In_One_Layer[i+1][j+1].z=IWC.z+j;
            }
    }
    if(abs(axis)==2) {
        N_Layers_1D=NC_y;
        Position_1D=X_Active_Bead.y;
        Size_1D=Ly;
        dC.x=0;
        dC.y=sign;
        dC.z=0;
        for(int i=-1;i<=+1;i++)
            for(int j=-1;j<=+1;j++) {
                Cell_In_One_Layer[i+1][j+1].x=IWC.x+i;
                Cell_In_One_Layer[i+1][j+1].y=IWC.y+sign;
                Cell_In_One_Layer[i+1][j+1].z=IWC.z+j;
            }
    }
    if(abs(axis)==3) {
        N_Layers_1D=NC_z;
        Position_1D=X_Active_Bead.z;
        Size_1D=Lz;
        dC.x=0;
        dC.y=0;
        dC.z=sign;
        for(int i=-1;i<=+1;i++)
            for(int j=-1;j<=+1;j++) {
                Cell_In_One_Layer[i+1][j+1].x=IWC.x+i;
                Cell_In_One_Layer[i+1][j+1].y=IWC.y+j;
                Cell_In_One_Layer[i+1][j+1].z=IWC.z+sign;
            }
    }

    double target_position_1D;
    int Ie,Is,dI;
    int N_layers_Visited=3;
    int N_layers_needed;
    
    target_position_1D=Position_1D+Event.Event_Time*sign;
    if(target_position_1D>Size_1D)target_position_1D-=Size_1D;
    if(target_position_1D<0      )target_position_1D+=Size_1D;
    Is=Position_1D/dL;
    Ie=target_position_1D/dL;
    if(Is==N_Layers_1D)Is--;
    if(Ie==N_Layers_1D)Ie--;
    dI=(Ie-Is)*sign;
    if(dI>0){
        N_layers_needed=dI+3;
    }else if(dI<0){
        N_layers_needed=min(N_Layers_1D,N_Layers_1D+3+dI);
    }else {//case: Ie==Is
        N_layers_needed=((Event.Event_Time>Size_1D/2)?N_Layers_1D:3);
    }
    if(Event.Event_Time>Size_1D)N_layers_needed=N_Layers_1D;
    
    while(N_layers_needed>N_layers_Visited) {
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++) {
                Cell_In_One_Layer[i][j].x+=dC.x;
                Cell_In_One_Layer[i][j].y+=dC.y;
                Cell_In_One_Layer[i][j].z+=dC.z;
                //Event_with_Cell(axis, T, NA, Cell_In_One_Layer[i][j]);
                Event_with_Cell(Event, Cell_In_One_Layer[i][j], Active_Bead, X_Active_Bead, axis);
            }
        N_layers_Visited++;

        //target_position_1D=Position_1D+T*sign;
        target_position_1D=Position_1D+Event.Event_Time*sign;
        if(target_position_1D>Size_1D)target_position_1D-=Size_1D;
        if(target_position_1D<0      )target_position_1D+=Size_1D;
        Is=Position_1D/dL;
        Ie=target_position_1D/dL;
        if(Is==N_Layers_1D)Is--;
        if(Ie==N_Layers_1D)Ie--;
        dI=(Ie-Is)*sign;
        if(dI>0){
            N_layers_needed=dI+3;
        }else if(dI<0){
            N_layers_needed=min(N_Layers_1D,N_Layers_1D+3+dI);
        }else {//case: Ie==Is
            //N_layers_needed=((T>Size_1D/2)?N_Layers_1D:3);
            N_layers_needed=((Event.Event_Time>Size_1D/2)?N_Layers_1D:3);
        }
        //if(T>Size_1D)N_layers_needed=N_Layers_1D;
        if(Event.Event_Time>Size_1D)N_layers_needed=N_Layers_1D;
    }
}



CellList::CellList(double Lx, double Ly, double Lz, double dL, int type_id, vector<double4>*X_pointer, double (*Gen)(double4,double4,int,double*,double),double *Params) {
	this->Lx=Lx;
	this->Ly=Ly;
	this->Lz=Lz;
	this->dL=dL;
    this->type_id=type_id;
    this->Gen=Gen;
    this->Params=Params;
    this->X_pointer=X_pointer;
	NC_x=(int)(Lx/dL);//for x>NC_x*dL, it is absorbed into last cell
	NC_y=(int)(Ly/dL);//for y>NC_y*dL, it is absorbed into last cell
	NC_z=(int)(Lz/dL);//for z>NC_z*dL, it is absorbed into last cell

	//Now construct Cell List
	//vector<vector<int3> >InWhichCell;
	//InWhichCell[type_id][bead_id] stores the in which cell each bead is
	//vector<vector<vector<Cell> > >Cells;	
	Cell Empty_Cell(type_id);
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
	InWhichCell.resize(X_pointer->size());
	for(int bead_id=0;bead_id<InWhichCell.size();bead_id++) {
			X=(*X_pointer)[bead_id];
			Ix=X.x/dL;
			Iy=X.y/dL;
			Iz=X.z/dL;
			
			if(Ix==NC_x)Ix=NC_x-1;
			if(Iy==NC_y)Iy=NC_y-1;
			if(Iz==NC_z)Iz=NC_z-1;

			InWhichCell[bead_id].x=Ix;
			InWhichCell[bead_id].y=Iy;
			InWhichCell[bead_id].z=Iz;

			Cells[Ix][Iy][Iz].Insert(bead_id);
	}

}


void CellList::Update(int2 ids, double4 New_X) {
    if(ids.x!=type_id)return;
	int3 IWC1,IWC2;

	bool update_tag=false;

	IWC1=InWhichCell[ids.y];
	IWC2=IWC1;


	IWC2.x=New_X.x/dL;
	if(IWC2.x==NC_x)IWC2.x=NC_x-1;
	if(IWC1.x!=IWC2.x)update_tag=true;

	IWC2.y=New_X.y/dL;
	if(IWC2.y==NC_y)IWC2.y=NC_y-1;
	if(IWC1.y!=IWC2.y)update_tag=true;

	IWC2.z=New_X.z/dL;
	if(IWC2.z==NC_z)IWC2.z=NC_z-1;
	if(IWC1.z!=IWC2.z)update_tag=true;

	if(update_tag){
		Cells[IWC1.x][IWC1.y][IWC1.z].Delete(ids.y);
		Cells[IWC2.x][IWC2.y][IWC2.z].Insert(ids.y);
		InWhichCell[ids.y]=IWC2;
	}

}

void CellList::print() {
	//Step1: Print Cell List
    cout<<"Type id:"<<type_id<<endl;
	bool mark;
	for(int i=0;i<NC_x;i++)
		for(int j=0;j<NC_y;j++)
			for(int k=0;k<NC_z;k++) {
				mark=false;
					for(int uu=0;uu<Cells[i][j][k].particle_list.size();uu++)
						{
							if(!mark){
								cout<<"In cell ("<<i<<','<<j<<','<<k<<")"<<endl;
								mark=true;
							}
							cout<<" id: "<<Cells[i][j][k].particle_list[uu]<<endl;
						}
			}		
	cout<<endl<<endl;

	//Step2: Print "In which cell"
	int3 i3;
    int2 ids;
	for(int bead_id=0;bead_id<X_pointer->size();bead_id++){
			cout<<"Position "<<(*X_pointer)[bead_id].x<<' '<<(*X_pointer)[bead_id].y<<' '<<(*X_pointer)[bead_id].z<<endl;
			ids.x=type_id;
            ids.y=bead_id;
            i3=In_Which_Cell(ids,(*X_pointer)[bead_id]);
			cout<<"In cell "<<i3.x<<' '<<i3.y<<' '<<i3.z<<endl<<endl;
	}
	cout<<endl;
}
