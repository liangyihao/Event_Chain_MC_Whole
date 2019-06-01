#include "basic.hpp"
#include "manage.hpp"
#include "Pairwise_Interaction.hpp"
#include "Random_Number.hpp"
#include "CellList.hpp"
#include <iostream>
#include <cmath>
#define BINNUM 100
using namespace std;
const double d0=7.5;
void Check(){
    cout<<"Checking"<<endl;
    double dx,dy,dz;
    for(int k0=0;k0<Types.size();k0++)
        for(int k1=0;k1<Types.size();k1++)
            for(int l0=0;l0<Types[k0].X.size();l0++)
                for(int l1=0;l1<Types[k1].X.size();l1++){
                    if((k0==k1)&&(l0==l1))continue;
                    dx=Types[k0].X[l0].x-Types[k1].X[l1].x;
                    dy=Types[k0].X[l0].y-Types[k1].X[l1].y;
                    dz=Types[k0].X[l0].z-Types[k1].X[l1].z;
                    if(dx>Lx/2)dx-=Lx;
                    if(dx<-Lx/2)dx+=Lx;

                    if(dy>Ly/2)dx-=Ly;
                    if(dy<-Ly/2)dy+=Ly;

                    if(dz>Lz/2)dz-=Lz;
                    if(dz<-Lz/2)dz+=Lz;

                    //cout<<"checking"<<k<<' '<<l<<endl;
                    if(dx*dx+dy*dy+dz*dz<d0*d0/4.0){
                        cout<<"error"<<endl;
                        cout<<k0<<' '<<l0<<' '<<k1<<' '<<l1<<endl;
                        cout<<Types[k0].X[l0].x<<','<<Types[k0].X[l0].y<<','<<Types[k0].X[l0].z<<endl;
                        cout<<Types[k1].X[l1].x<<','<<Types[k1].X[l1].y<<','<<Types[k1].X[l1].z<<endl;
                        cout<<sqrt(dx*dx+dy*dy+dz*dz)<<endl;
                        exit(0);
                    }
                
                }

    cout<<"Checking finished"<<endl;
}

void Init(){
    Lx=100;Ly=100;Lz=100;
    const int N=1000;
    rand_init(1);
    int tid1,tid2;
    tid1=Create_Charged_Type(1);//type A, hard sphere diameter 10
    tid2=Create_Charged_Type(-1);//type A, hard sphere diameter 10

    Create_Hard_Sphere_Interaction_Between_Types(tid1,tid1,d0);
    Create_Hard_Sphere_Interaction_Between_Types(tid1,tid2,d0);
    Create_Hard_Sphere_Interaction_Between_Types(tid2,tid2,d0);

    double4 X;
    double xc,yc,zc;

    int type=0;
    for(xc=0;xc<Lx-d0;xc+=1.01*d0)
        for(yc=0;yc<Ly-d0;yc+=1.01*d0)
            for(zc=0;zc<Lz-d0;zc+=1.01*d0) {
                if((Types[0].X.size()==N)&&(type==0))type=1;
                if((Types[1].X.size()==N)&&(type==1))break;
                X.x=xc;X.y=yc;X.z=zc;
                Create_Charged_Bead(type,X);
            }
    cout<<Types[0].X.size()<<' '<<Types[1].X.size()<<endl;
    Check();//For debug
    Global_Cell_List_Pointer=new CellList(Lx, Ly, Lz, MAX_SHORT_INTERACTION_RANGE, (&Types));
    initialize_sys_charge();
}
/*
double g11[BINNUM],g12[BINNUM],g22[BINNUM];
int sample_T=0;
void sample(bool ifstat,bool ifprint) {
	if((!ifstat)&&(!ifprint))return;
    double dr=Ly/BINNUM;;
    sample_T++;
    double dx,dy,dz,r;
    double Dx,Dy,Dz;
    double4 X,Y;
    double eedist=0;
    for(int l=0;l<end_to_end_ids.size();l++){
        Dx=0;Dy=0;Dz=0;
        for(int z=end_to_end_ids[l].x;z<end_to_end_ids[l].y;z++){
            X=Types[0].X[z];
            Y=Types[0].X[z+1];
            dx=X.x-Y.x;
            dy=X.y-Y.y;
            dz=X.z-Y.z;
            if(dx>Lx/2)dx-=Lx;
            if(dx<-Lx/2)dx+=Lx;
            if(dy>Ly/2)dy-=Ly;
            if(dy<-Ly/2)dy+=Ly;
            if(dz>Lz/2)dz-=Lz;
            if(dz<-Lz/2)dz+=Lz;
            Dx+=dx;
            Dy+=dy;
            Dz+=dz;
        }
        r=sqrt(Dx*Dx+Dy*Dy+Dz*Dz);
        if(ifstat)if(r/dr<BINNUM)g[(int)(r/dr)]++;
        eedist+=r;
    }
    eedist/=end_to_end_ids.size();
    if(ifprint)cout<<"Average e-e dist over chains: "<<eedist<<endl;
}*/

int main(){
    Init();
    Global_Cell_List_Pointer->print();
    cout<<endl<<endl<<endl;

    double ifstat=false,ifprint=false;
//    for(int l=0;l<Cell_Veto_Lists.size();l++)Cell_Veto_Lists[l]->check_print();//For debug
//    for(int l=0;l<Cell_Veto_Lists.size();l++)Cell_Veto_Lists[l]->check_rate();//For debug
    
    for(int k=0;k<2000000;k++) {
        //if(k%100000==0)
        //Check();//For debug
//        cout<<"Moving along +x"<<endl;
        Monte_Carlo(800,1);
//        for(int l=0;l<Cell_Veto_Lists.size();l++)Cell_Veto_Lists[l]->check_print();//For debug

//        Global_Cell_List_Pointer->print();
        //Check();//For debug
//        cout<<"Moving along +y"<<endl;
        Monte_Carlo(800,2);
//        for(int l=0;l<Cell_Veto_Lists.size();l++)Cell_Veto_Lists[l]->check_print();//For debug

//        Global_Cell_List_Pointer->print();
//        Check();//For debug        
//        cout<<"Moving along +z"<<endl;
        Monte_Carlo(800,3);
//        for(int l=0;l<Cell_Veto_Lists.size();l++)Cell_Veto_Lists[l]->check_print();//For debug

//        Global_Cell_List_Pointer->print();
        
//        cout<<"Moving along -x"<<endl;
        Monte_Carlo(800,-1);
//        for(int l=0;l<Cell_Veto_Lists.size();l++)Cell_Veto_Lists[l]->check_print();//For debug

//        Global_Cell_List_Pointer->print();

//        cout<<"Moving along -y"<<endl;        
        Monte_Carlo(800,-2);
//        for(int l=0;l<Cell_Veto_Lists.size();l++)Cell_Veto_Lists[l]->check_print();//For debug

//        Global_Cell_List_Pointer->print();

//        cout<<"Moving along -z"<<endl;        
        Monte_Carlo(800,-3);
//        for(int l=0;l<Cell_Veto_Lists.size();l++)Cell_Veto_Lists[l]->check_print();//For debug

//        Global_Cell_List_Pointer->print();
        ifprint=false;
        
//        cout<<"k="<<k<<" finished"<<endl;
        if(k%10==0){
            cout<<"k="<<k<<" finished"<<endl;
            //Global_Cell_List_Pointer->print();
            for(int l=0;l<Cell_Veto_Lists.size();l++)Cell_Veto_Lists[l]->check_print();//For debug
            ifprint=true;
        }
        if(k>300000)ifstat=true;
        //if(k%100==0)sample(ifstat,ifprint);
        
    }
/*
    double dr=Ly/BINNUM;;
    for(int k=0;k<BINNUM;k++) {
        double r1,r2;
        r1=dr*k;
        r2=dr*(k+1);
        cout<<dr*k<<" "<<g[k]/(4*M_PI/3*(r2*r2*r2-r1*r1*r1)*sample_T)<<endl;
        //cout<<dr*(k-BINNUM/2)<<" "<<g[k]<<endl;
    }
*/
}
