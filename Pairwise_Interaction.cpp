/*
In this file, we define the event time routines for basic pairwise interactions
*/
#include <vector>
#include <cmath>
#include <iostream>
#include "Random_Number.hpp"
#include "basic.hpp"
#include "Pairwise_Interaction.hpp"
#include "manage.hpp"
using namespace std;


double Event_Time_Hard_Sphere(double4 X1,double4 X2,int axis_index,double*Params, double Max_Event_T){
    //Move X1 along axis_index direction, get reject event time of hard sphere interaction with X2, with parameters Params.
    //If the event time is less than Max_Event_T and less than period in the moving direction, return real event time.
    //If the event time is greater than Max_Event_T, return 2*(period of moving direction)
    //If the event time is greater than one period of moving direction, return 2*(period of moving direction)

    double Lx,Ly,Lz,d;

    Lx=Params[0];
    Ly=Params[1];
    Lz=Params[2];
    d=Params[3];

    double L,Lyy,Lzz,min_d2;
    double dxx,dyy,dzz;
    int sign=axis_index/abs(axis_index);
    if(abs(axis_index)==1){
        L=Lx;
        dyy=X2.y-X1.y;Lyy=Ly;
        dzz=X2.z-X1.z;Lzz=Lz;
        dxx=(X2.x-X1.x)*sign;
    }
    if(abs(axis_index)==2){
        L=Ly;
        dyy=X2.z-X1.z;Lyy=Lz;
        dzz=X2.x-X1.x;Lzz=Lx;
        dxx=(X2.y-X1.y)*sign;
    }
    if(abs(axis_index)==3){
        L=Lz;
        dyy=X2.x-X1.x;Lyy=Lx;
        dzz=X2.y-X1.y;Lzz=Ly;
        dxx=(X2.z-X1.z)*sign;
    }

    if(dyy>+Lyy/2)dyy-=Lyy;
    if(dyy<-Lyy/2)dyy+=Lyy;

    if(dzz>+Lzz/2)dzz-=Lzz;
    if(dzz<-Lzz/2)dzz+=Lzz;

    if(dxx<0)dxx+=L;
    min_d2=dyy*dyy+dzz*dzz;
    if(min_d2-d*d>-1E-20)
    //if(min_d2>d*d)
        return 2*L;
    else
        return dxx-sqrt(d*d-min_d2);

    return 2*L;
}

double Event_Time_Spring(double4 X1,double4 X2,int axis_index,double*Params, double Max_Event_T) {
    //Move X1 along axis_index direction, get reject event time of spring interaction with X2, with parameters Params.
    //If the event time is less than Max_Event_T and less than period in the moving direction, return real event time.
    //If the event time is greater than Max_Event_T, return 2*(period of moving direction)
    //If the event time is greater than one period of moving direction, return 2*(period of moving direction)

    double Lx,Ly,Lz;
    double K,l0;

    Lx=Params[0];
    Ly=Params[1];
    Lz=Params[2];
    K =Params[3];
    l0=Params[4];

    double min_dr2;
    double dxx,dyy,dzz;
    double Lxx,Lyy,Lzz;

    int sign=axis_index/abs(axis_index);

    if(abs(axis_index)==1){
        Lxx=Lx;
        Lyy=Ly;
        Lzz=Lz;
        dyy=X2.y-X1.y;
        dzz=X2.z-X1.z;
        dxx=X2.x-X1.x;
    }
    if(abs(axis_index)==2){
        Lxx=Ly;
        Lyy=Lz;
        Lzz=Lx;
        dyy=X2.z-X1.z;
        dzz=X2.x-X1.x;
        dxx=X2.y-X1.y;
    }
    if(abs(axis_index)==3){
        Lxx=Lz;
        Lyy=Lx;
        Lzz=Ly;
        dyy=X2.x-X1.x;
        dzz=X2.y-X1.y;
        dxx=X2.z-X1.z;
    }
    if(dxx<-Lxx/2) dxx+=Lxx;
    if(dxx>+Lxx/2) dxx-=Lxx;
    if(dyy<-Lyy/2) dyy+=Lyy;
    if(dyy>+Lyy/2) dyy-=Lyy;
    if(dzz<-Lzz/2) dzz+=Lzz;
    if(dzz>+Lzz/2) dzz-=Lzz;
    min_dr2=dyy*dyy+dzz*dzz;

    dxx*=sign;
    //Do transformation, make the moving along +x axis

    double q;
    q=-log(1-Uniform_Random());
    /*
    q=Uniform_Random();
    if(q==0)
        return 2*Lxx;
    else
        q=-log(q);
    */
    double s1,s2;
    double dr_1,dr_2;
    double q1,q2;

    //Set s1,s2,q1,q2, or just return 2*Lxx
    if(min_dr2>=l0*l0) {
        s2=Lxx;
        s1=(dxx>0)?dxx:0;
        dr_1=sqrt((dxx-s1)*(dxx-s1)+min_dr2);
        dr_2=sqrt((dxx-s2)*(dxx-s2)+min_dr2);
        q1=0;
        q2=0.5*K*((dr_2-l0)*(dr_2-l0)-(dr_1-l0)*(dr_1-l0))+q1;
    }else {
        double x0;
        x0=sqrt(l0*l0-min_dr2);
        if(dxx>0) {
            //1st: try 0~dxx
                s2=dxx;
                s1=(dxx-x0>0)?dxx-x0:0;
                dr_1=sqrt((dxx-s1)*(dxx-s1)+min_dr2);
                dr_2=sqrt((dxx-s2)*(dxx-s2)+min_dr2);
                q1=0;
                q2=0.5*K*((dr_2-l0)*(dr_2-l0)-(dr_1-l0)*(dr_1-l0))+q1;
            //2nd: try dxx->Lxx
                if(q2<q) {
                    s2=Lxx;
                    s1=x0+dxx;
                    dr_1=sqrt((dxx-s1)*(dxx-s1)+min_dr2);
                    dr_2=sqrt((dxx-s2)*(dxx-s2)+min_dr2);
                    q1=q2;
                    q2=0.5*K*((dr_2-l0)*(dr_2-l0)-(dr_1-l0)*(dr_1-l0))+q1;
                }
        }else {
            s2=Lxx;
            s1=(x0+dxx>0)?x0+dxx:0;
            dr_1=sqrt((dxx-s1)*(dxx-s1)+min_dr2);
            dr_2=sqrt((dxx-s2)*(dxx-s2)+min_dr2);
            q1=0;
            q2=0.5*K*((dr_2-l0)*(dr_2-l0)-(dr_1-l0)*(dr_1-l0))+q1;
        }
    }

    if(q2<q)return 2*Lxx;

    //Find solution
    double sm,qm,dr_m;
    double s10,q10,dr_10;
    s10=s1;q10=q1;dr_10=dr_1;
    double bas;
    bas=q10-0.5*K*(dr_10-l0)*(dr_10-l0);
    while((abs(q2-q1)>1E-10)||(abs(s2-s1)>1E-10)) {
        if(s1>Max_Event_T)return 2*Lxx;//This trick is for saving time
        sm=(s1+s2)/2;
        dr_m=sqrt((dxx-sm)*(dxx-sm)+min_dr2);
        qm=0.5*K*(dr_m-l0)*(dr_m-l0)+bas;
        if(qm>q){
            q2=qm;
            s2=sm;
            dr_2=dr_m;
        }else {
            q1=qm;
            s1=sm;
            dr_1=dr_m;
        }
    }

    return sm;
}


void Create_Hard_Sphere_Interaction_Between_Types(int Type_id1,int Type_id2,double d){
    if((d>Lx/10)||(d>Ly/10)||(d>Lz/10)){
        cout<<"Hard sphere Error: interaction distance should be smaller than 1/10 of system size"<<endl;
        exit(1);
    }
    //Check if Type_id1 and Type_id2 are available
    if(max(Type_id1,Type_id2)>=Types.size()){
        cout<<"Hard sphere Error: type doesn't exist"<<endl;
        exit(1);
    }

    //Register interaction range
    if(MAX_SHORT_INTERACTION_RANGE<d)MAX_SHORT_INTERACTION_RANGE=d;

    Parameter_List Param(Lx,Ly,Lz);
    Param.data[3]=d;
    //register interaction
    Event_Time_Generator_List.push_back(Event_Time_Hard_Sphere);
    Param_Lists.push_back(Param);
    int Interaction_Global_ID;
    Interaction_Global_ID=Event_Time_Generator_List.size()-1;

    //connect two types
    int2 temp;
    temp.x=Type_id2;temp.y=Interaction_Global_ID;
    Types[Type_id1].Interactions_with_Types.push_back(temp);

    if(Type_id1==Type_id2)return;
    temp.x=Type_id1;
    Types[Type_id2].Interactions_with_Types.push_back(temp);
    return;
}

void Create_Spring_Interaction_Between_Beads(int2 Bead_id1,int2 Bead_id2,double k,double l0){
    if(0.5*k*l0*l0<10) {
        cout<<"Spring Error: k too weak"<<endl;
    }
    if((l0>Lx/10)||(l0>Ly/10)||(l0>Lz/10)){
        cout<<"Spring Error: spring original length should be smaller than 1/10 of system size"<<endl;
        exit(1);
    }
    //Check if Bead_id1 and Bead_id2 are available
    if(max(Bead_id1.x,Bead_id2.x)>=Types.size()){
        cout<<"Spring Error: type doesn't exist"<<endl;
        exit(1);
    }
    if(Types[Bead_id1.x].X.size()<=Bead_id1.y){
        cout<<"Spring Error: bead 1 doesn't exist"<<endl;
        exit(1);
    }
    if(Types[Bead_id2.x].X.size()<=Bead_id2.y){
        cout<<"Spring Error: bead 2 doesn't exist"<<endl;
        exit(1);
    }
    if((Bead_id1.x==Bead_id2.x)&&(Bead_id1.y==Bead_id2.y)){
        cout<<"Spring Error: bead 1 and bead 2 should be different"<<endl;
        exit(1);
    }

    Parameter_List Param(Lx,Ly,Lz);
    Param.data[3]=k;
    Param.data[4]=l0;
    //register interaction
    Event_Time_Generator_List.push_back(Event_Time_Spring);
    Param_Lists.push_back(Param);
    int Interaction_Global_ID;
    Interaction_Global_ID=Event_Time_Generator_List.size()-1;
    
    //connect two beads
    int3 temp;
    temp.z=Interaction_Global_ID;
    int Gid1=Bead_id1.x;
    int Bid1=Bead_id1.y;
    int Gid2=Bead_id2.x;
    int Bid2=Bead_id2.y;

    temp.x=Gid2;temp.y=Bid2;
    Types[Gid1].Interactions_with_Beads[Bid1].push_back(temp);
    temp.x=Gid1;temp.y=Bid1;
    Types[Gid2].Interactions_with_Beads[Bid2].push_back(temp);

}
