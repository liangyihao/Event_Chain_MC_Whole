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

vector<Bead_Type> Types;
vector<double(*)(double4,double4,int,double*,double)> Event_Time_Generator_List_For_Bonds;
vector<Parameter_List> Param_Lists_For_Bonds;
int2 Active_Bead;
CellList* Global_Cell_List_Pointer;
vector<Short_Range_Interaction_Between_Types*>Short_Range_Interaction_Between_Types_List;

//////////////////////For Pressure Computation////////////////////////////////////////////////////////////////////////////////////////
double Pressure=0;//excessive pressure
int Pressure_Count=0;
double Event_Count_Per_Particle=0;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Get_Event(double&time, int2&id_next_active_bead, int axis){
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
    double Pressure_Chain=0;
    int Event_Count=0;
    double temp;
    double4 NX;
    //Need to randomly choose an active particle, implement it later
    int total_particle_number=0;
    for(int i=0;i<Types.size();i++)total_particle_number+=Types[i].X.size();
    int index=(int)(Uniform_Random()*total_particle_number);
    for(int j=0;j<Types.size();j++){
        if (index<Types[j].X.size()){
            Active_Bead.x=j;
            Active_Bead.y=index;
        }
        else{
            index-=Types[j].X.size();
        }
    }
    //cout<<"First Active bead: "<<Active_Bead.x<<" "<<Active_Bead.y<<endl;

    if((abs(axis)>3)||(axis==0)){cout<<"Error, axis id should be -3,-2,-1,1,2,3"<<endl;return;}
    int2 id_next_active_bead;
    bool go_ahead=true;
    int sign=axis/abs(axis);
    while(go_ahead) {
        Get_Event(time,id_next_active_bead,axis);
        
        if(abs(axis)==1)exe_time=min(Lx,time);//axis==1: on +x direction, axis==-1: on -x direction
        if(abs(axis)==2)exe_time=min(Ly,time);//axis==2: on +y direction, axis==-2: on -y direction
        if(abs(axis)==3)exe_time=min(Lz,time);//axis==3: on +z direction, axis==-3: on -z direction

        if(Clock+exe_time>=Stop_Clock) {
                go_ahead=false;
                exe_time=Stop_Clock-Clock;
        }else if((Active_Bead.x!=id_next_active_bead.x)||(Active_Bead.y!=id_next_active_bead.y)){
            //compute pressure contribution
            double dx_p;
            if(abs(axis)==1){
                dx_p=Types[id_next_active_bead.x].X[id_next_active_bead.y].x - Types[Active_Bead.x].X[Active_Bead.y].x-exe_time*sign;
                while(dx_p>+Lx/2)dx_p-=Lx;
                while(dx_p<-Lx/2)dx_p+=Lx;
                dx_p*=sign;
            }
            if(abs(axis)==2){
                dx_p=Types[id_next_active_bead.x].X[id_next_active_bead.y].y - Types[Active_Bead.x].X[Active_Bead.y].y-exe_time*sign;
                while(dx_p>+Ly/2)dx_p-=Ly;
                while(dx_p<-Ly/2)dx_p+=Ly;
                dx_p*=sign;
            }
            if(abs(axis)==3){
                dx_p=Types[id_next_active_bead.x].X[id_next_active_bead.y].z - Types[Active_Bead.x].X[Active_Bead.y].z-exe_time*sign;
                while(dx_p>+Lz/2)dx_p-=Lz;
                while(dx_p<-Lz/2)dx_p+=Lz;
                dx_p*=sign;
            }
            Pressure_Chain+=dx_p;
            Event_Count++;
        }
        
        //Global_Cell_List_Pointer->Move(Active_Bead, exe_time, axis);
        double4 New_X;
        New_X=Types[Active_Bead.x].X[Active_Bead.y];
        if(abs(axis)==1){
            New_X.x+=exe_time*sign;
            while(New_X.x>=Lx)New_X.x-=Lx;
            while(New_X.x<0 )New_X.x+=Lx;
        }
        if(abs(axis)==2){
            New_X.y+=exe_time*sign;
            while(New_X.y>=Ly)New_X.y-=Ly;
            while(New_X.y<0 )New_X.y+=Ly;
        }
        if(abs(axis)==3){
            New_X.z+=exe_time*sign;
            while(New_X.z>=Lz)New_X.z-=Lz;
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
    Pressure_Chain/=Stop_Clock;
    Pressure+=Pressure_Chain;
    Pressure_Count++;
    int N_tot=0;
    for(int k=0;k<Types.size();k++)N_tot+=Types[k].X.size();
    Event_Count_Per_Particle+=((double)Event_Count)/N_tot;
}
