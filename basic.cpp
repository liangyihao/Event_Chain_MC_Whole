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
vector<double(*)(double4,double4,int,double*,double)> Event_Time_Generator_List_For_Bonds;
vector<Parameter_List> Param_Lists_For_Bonds;
int2 Active_Bead;
CellList* Global_Cell_List_Pointer;
vector<Short_Range_Interaction_Between_Types*>Short_Range_Interaction_Between_Types_List;


void Get_Event(double&time, int2&id_next_active_bead, int axis) {
    //Get Event within time Lx if axis==1(or Ly if axis==2, or Lz if axis==3)
    //If event time out of bound L(=Lx or Ly or Lz),let time=2*L
    //Normally, time is the event time and id_next_active_bead is the id of next active bead

    Bead_Type*Type1;
    double4 X_Active_Bead,Y;

    TwoBody_Event Event;
    Event.Target_Bead=Active_Bead;
    if(abs(axis)==1)Event.Event_Time=2*Lx;
    if(abs(axis)==2)Event.Event_Time=2*Ly;
    if(abs(axis)==3)Event.Event_Time=2*Lz;


    double (*Gen)(double4,double4,int,double*,double);
    double *Params;
    Type1=&Types[Active_Bead.x];
    X_Active_Bead=Type1->X[Active_Bead.y];

    //Special bead-bead interactions
    int3 ids3;
    double t;
    for(int k=0;k<Type1->Interactions_with_Beads[Active_Bead.y].size();k++){
        ids3=Type1->Interactions_with_Beads[Active_Bead.y][k];
        Y=Types[ids3.x].X[ids3.y];
        Gen=Event_Time_Generator_List_For_Bonds[ids3.z];
        Params=Param_Lists_For_Bonds[ids3.z].data;
        t=Gen(X_Active_Bead,Y,axis,Params,Event.Event_Time);
        if(t<Event.Event_Time){
                Event.Event_Time=t;
                Event.Target_Bead.x=ids3.x;
                Event.Target_Bead.y=ids3.y;
        }
    }

    //Type-Type interactions
    for(int k=0;k<Type1->Interactions_with_Types.size();k++){
        Short_Range_Interaction_Between_Types_List[Type1->Interactions_with_Types[k]]->Get_Event(Event,Active_Bead,X_Active_Bead,axis);
    }

    //Coloumb interaction
    for(int k=0;k<Cell_Veto_Lists.size();k++)Cell_Veto_Lists[k]->Get_Colomb_Event(Active_Bead,axis,Bjerrum_Length,Event.Event_Time,Event.Target_Bead);
    time=Event.Event_Time;
    id_next_active_bead=Event.Target_Bead;
}



void Monte_Carlo(double Stop_Clock,int axis) {
    double Clock=0;
    double time,exe_time;
    double temp;
    double4 NX;
    //Need to randomly choose an active particle, implement it later
        Active_Bead.x=(int)(Uniform_Random()*Types.size());
        Active_Bead.y=(int)(Uniform_Random()*Types[Active_Bead.x].X.size());
        //cout<<"First Active bead: "<<Active_Bead.x<<" "<<Active_Bead.y<<endl;

    if((abs(axis)>3)||(axis==0)){cout<<"Error, axis id should be -3,-2,-1,1,2,3"<<endl;return;}
    int2 id_next_active_bead;
    bool go_ahead=true;
    while(go_ahead) {
        //cout<<endl<<"Active: "<<Active_Bead.x<<" "<<Active_Bead.y<<endl;
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
        
        //Global_Cell_List_Pointer->Move(Active_Bead, exe_time, axis);
        double4 New_X;
        New_X=Types[Active_Bead.x].X[Active_Bead.y];
        if(abs(axis)==1){
            New_X.x+=exe_time*(axis/abs(axis));
            while(New_X.x>Lx)New_X.x-=Lx;
            while(New_X.x<0 )New_X.x+=Lx;
        }
        if(abs(axis)==2){
            New_X.y+=exe_time*(axis/abs(axis));
            while(New_X.y>Ly)New_X.y-=Ly;
            while(New_X.y<0 )New_X.y+=Ly;
        }
        if(abs(axis)==3){
            New_X.z+=exe_time*(axis/abs(axis));
            while(New_X.z>Lz)New_X.z-=Lz;
            while(New_X.z<0 )New_X.z+=Lz;
        }
        for(int k=0;k<Types[Active_Bead.x].Interactions_with_Types.size();k++){
            Short_Range_Interaction_Between_Types_List[Types[Active_Bead.x].Interactions_with_Types[k]]->Update(Active_Bead,New_X);
        }
        Types[Active_Bead.x].X[Active_Bead.y]=New_X;

        //!!!For all cell-veto list, update
        //NX=Types[Active_Bead.x].X[Active_Bead.y];
        for(int k=0;k<Cell_Veto_Lists.size();k++)Cell_Veto_Lists[k]->Update(Active_Bead,New_X);

        Clock+=exe_time;
        Active_Bead=id_next_active_bead;
    }
}
