/*
Author: Yihao Liang
liangyihaosjtu@gmail.com
This code is for Event Chain Monte Carlo for pairwise interacting many body system
*/
#include<vector>
#include<iostream>
//#include"omp.h"
#include"basic.hpp"
#include"manage.hpp"
#include"Random_Number.hpp"
#include<cmath>
using namespace std;
/*
//DEBUG///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Get_Event_BF(double&time, int2&id_next_active_bead, int axis);//FOR DEBUG
//DEBUG///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
*/
vector<Bead_Type> Types;
vector<double(*)(double4,double4,int,double*,double)> Event_Time_Generator_List;
vector<Parameter_List> Param_Lists;
int2 Active_Bead;
double MAX_SHORT_INTERACTION_RANGE=0;
CellList* Global_Cell_List_Pointer;

void Event_with_Cell(int axis, double&time, int2&id_veto_bead, int3 Cell_ID){
    Bead_Type*Type1;
    Bead_Type*Type2;
    double4 X,Y;

    double (*Gen)(double4,double4,int,double*,double);
    double *Params;
    Type1=&Types[Active_Bead.x];
    X=Type1->X[Active_Bead.y];
    vector<int>*Particle_List_Pointer;
    int bead_id;
    int target_type_id,interaction_id;
    double t;
    for(int z=0;z<Type1->Interactions_with_Types.size();z++){
        target_type_id=Type1->Interactions_with_Types[z].x;
        interaction_id=Type1->Interactions_with_Types[z].y;
        Type2=&Types[target_type_id];
        Gen=Event_Time_Generator_List[interaction_id];
        Params=Param_Lists[interaction_id].data;
        Particle_List_Pointer=Global_Cell_List_Pointer->Particle_List_In_Cell(Cell_ID,target_type_id);
        for(int k=0;k<(*Particle_List_Pointer).size();k++) {
            //cout<<"Tri"<<endl;
            bead_id=(*Particle_List_Pointer)[k];
            if(((Type1->index)==(Type2->index))&&(Active_Bead.y==bead_id)) continue;
            Y=Type2->X[bead_id];
            t=Gen(X,Y,axis,Params,time);
            if(t<time){
                time=t;
                id_veto_bead.x=Type2->index;
                id_veto_bead.y=bead_id;
            }
        }
    }
}

void Get_Event(double&time, int2&id_next_active_bead, int axis) {
    //Get Event within time Lx if axis==1(or Ly if axis==2, or Lz if axis==3)
    //If event time out of bound L(=Lx or Ly or Lz),let time=2*L
    //Normally, time is the event time and id_next_active_bead is the id of next active bead

    Bead_Type*Type1;
    double4 X,Y;
    int2 NA;
    NA=Active_Bead;
    double T;
    if(abs(axis)==1)T=2*Lx;
    if(abs(axis)==2)T=2*Ly;
    if(abs(axis)==3)T=2*Lz;

    double (*Gen)(double4,double4,int,double*,double);
    double *Params;
    Type1=&Types[Active_Bead.x];
    X=Type1->X[Active_Bead.y];


    //Special bead-bead interactions
    int3 ids3;
    double t;
    for(int k=0;k<Type1->Interactions_with_Beads[Active_Bead.y].size();k++){
        ids3=Type1->Interactions_with_Beads[Active_Bead.y][k];
        Y=Types[ids3.x].X[ids3.y];
        Gen=Event_Time_Generator_List[ids3.z];
        Params=Param_Lists[ids3.z].data;
        t=Gen(X,Y,axis,Params,T);
        if(t<T){
                T=t;
                NA.x=ids3.x;
                NA.y=ids3.y;
        }
    }

    //Type-Type interactions
    int3 IWC;

    IWC=Global_Cell_List_Pointer->In_Which_Cell(Active_Bead);
    
    //Run first layer
    int3 Cell_ID;
    for(Cell_ID.x=IWC.x-1;Cell_ID.x<=IWC.x+1;Cell_ID.x++)
        for(Cell_ID.y=IWC.y-1;Cell_ID.y<=IWC.y+1;Cell_ID.y++)
            for(Cell_ID.z=IWC.z-1;Cell_ID.z<=IWC.z+1;Cell_ID.z++){
                Event_with_Cell(axis, T, NA, Cell_ID);
            }

    //Run additional layers, if needed
    double Position_1D,Size_1D;
    int N_Layers_1D;//Number of layers on selected dimension

    int3 Cell_In_One_Layer[3][3];
    int3 dC;
    int sign;
    sign=(axis>0)?1:-1;
    if(abs(axis)==1) {
        N_Layers_1D=Global_Cell_List_Pointer->NCell_x();
        Position_1D=X.x;
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
        N_Layers_1D=Global_Cell_List_Pointer->NCell_y();
        Position_1D=X.y;
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
        N_Layers_1D=Global_Cell_List_Pointer->NCell_z();
        Position_1D=X.z;
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
    
    target_position_1D=Position_1D+T*sign;
    if(target_position_1D>Size_1D)target_position_1D-=Size_1D;
    if(target_position_1D<0      )target_position_1D+=Size_1D;
    Is=Position_1D/(Global_Cell_List_Pointer->Cell_size());
    Ie=target_position_1D/(Global_Cell_List_Pointer->Cell_size());
    if(Is==N_Layers_1D)Is--;
    if(Ie==N_Layers_1D)Ie--;
    dI=(Ie-Is)*sign;
    if(dI>0){
        N_layers_needed=dI+3;
    }else if(dI<0){
        N_layers_needed=min(N_Layers_1D,N_Layers_1D+3+dI);
    }else {//case: Ie==Is
        N_layers_needed=((T>Size_1D/2)?N_Layers_1D:3);
    }
    if(T>Size_1D)N_layers_needed=N_Layers_1D;
    
    
    while(N_layers_needed>N_layers_Visited) {
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++) {
                Cell_In_One_Layer[i][j].x+=dC.x;
                Cell_In_One_Layer[i][j].y+=dC.y;
                Cell_In_One_Layer[i][j].z+=dC.z;
                Event_with_Cell(axis, T, NA, Cell_In_One_Layer[i][j]);
            }
        N_layers_Visited++;

        target_position_1D=Position_1D+T*sign;
        if(target_position_1D>Size_1D)target_position_1D-=Size_1D;
        if(target_position_1D<0      )target_position_1D+=Size_1D;
        Is=Position_1D/(Global_Cell_List_Pointer->Cell_size());
        Ie=target_position_1D/(Global_Cell_List_Pointer->Cell_size());
        if(Is==N_Layers_1D)Is--;
        if(Ie==N_Layers_1D)Ie--;
        dI=(Ie-Is)*sign;
        if(dI>0){
            N_layers_needed=dI+3;
        }else if(dI<0){
            N_layers_needed=min(N_Layers_1D,N_Layers_1D+3+dI);
        }else {//case: Ie==Is
            N_layers_needed=((T>Size_1D/2)?N_Layers_1D:3);
        }
        if(T>Size_1D)N_layers_needed=N_Layers_1D;
    }


    time=T;
    id_next_active_bead=NA;
}



void Monte_Carlo(double Stop_Clock,int axis) {
    double Clock=0;
    double time,exe_time;
    double temp;
    //Need to randomly choose an active particle, implement it later
        Active_Bead.x=(int)(Uniform_Random()*Types.size());
        Active_Bead.y=(int)(Uniform_Random()*Types[Active_Bead.x].X.size());
        //cout<<"First Active bead: "<<Active_Bead.x<<" "<<Active_Bead.y<<endl;

    if((abs(axis)>3)||(axis==0)){cout<<"Error, axis id should be -3,-2,-1,1,2,3"<<endl;return;}
    int2 id_next_active_bead;
    bool go_ahead=true;
    while(go_ahead) {
        //cout<<"Active: "<<Active_Bead.x<<" "<<Active_Bead.y<<endl;
        //cout<<Types[0].X[0].x<<' '<<Types[0].X[0].y<<' '<<Types[0].X[0].z<<endl;
        //cout<<Types[1].X[0].x<<' '<<Types[1].X[0].y<<' '<<Types[1].X[0].z<<endl<<endl;
        Get_Event(time,id_next_active_bead,axis);
/*
//DEBUG//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        double time2;int2 id_next_active_bead2;//DEBUG
        Get_Event_BF(time2,id_next_active_bead2,axis);//DEBUG
        if(abs(time-time2)>0.000001){
            cout<<"Wrong in event time : Cell List::"<<time<<" Real::"<<time2<<endl;
            cout<<"Current Active::"<<Active_Bead.x<<' '<<Active_Bead.y<<endl;
            cout<<"Next Active(cell list)::"<<id_next_active_bead.x<<' '<<id_next_active_bead.y<<endl;
            cout<<"Next Active(real     )::"<<id_next_active_bead2.x<<' '<<id_next_active_bead2.y<<endl;
            Global_Cell_List_Pointer->print();
            cout<<"Direction: "<<axis<<endl;
            exit(0);
        }//DEBUG
        if(id_next_active_bead.x!=id_next_active_bead2.x){cout<<"Wrong in next bead"<<endl;exit(0);}//DEBUG
        if(id_next_active_bead.y!=id_next_active_bead2.y){cout<<"Wrong in next bead"<<endl;exit(0);}//DEBUG
//DEBUG//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
*/
        //cout<<"Direction:"<<axis<<" Clock:"<<Clock<<"Current active bead: "<<Active_Bead.x<<','<<Active_Bead.y<<" time:"<<time<<endl;
        
        if(abs(axis)==1)exe_time=min(Lx,time);//axis==1: on +x direction, axis==-1: on -x direction
        if(abs(axis)==2)exe_time=min(Ly,time);//axis==2: on +y direction, axis==-2: on -y direction
        if(abs(axis)==3)exe_time=min(Lz,time);//axis==3: on +z direction, axis==-3: on -z direction

        if(Clock+exe_time>=Stop_Clock) {
                go_ahead=false;
                exe_time=Stop_Clock-Clock;
        }
        
        Global_Cell_List_Pointer->Move(Active_Bead, exe_time, axis);
        Clock+=exe_time;
        Active_Bead=id_next_active_bead;
    }
}
