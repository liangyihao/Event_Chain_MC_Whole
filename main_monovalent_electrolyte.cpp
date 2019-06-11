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
    const int N=1;
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

double g00[BINNUM],g01[BINNUM],g11[BINNUM];
double c0_x[BINNUM],c0_y[BINNUM],c0_z[BINNUM];
double c1_x[BINNUM],c1_y[BINNUM],c1_z[BINNUM];
int sample_T;
double dhx,dhy,dhz,dr;
void sample_init(){
    sample_T=0;
    dhx=Lx/BINNUM;
    dhy=Ly/BINNUM;
    dhz=Lz/BINNUM;
    dr=min(min(Lx,Ly),Lz)/BINNUM;
    for(int k=0;k<BINNUM;k++){
        g00[k]=0;
        g01[k]=0;
        g11[k]=0;

        c0_x[k]=0;
        c1_x[k]=0;
        c0_y[k]=0;
        c1_y[k]=0;
        c0_z[k]=0;
        c1_z[k]=0;
    }
}
void sample() {
    sample_T++;
    double dx,dy,dz,r;
    double4 X,Y;
    for(int i=0;i<Types[0].X.size();i++){
        X=Types[0].X[i];
        c0_x[(int)(X.x/dhx)]++;
        c0_y[(int)(X.y/dhy)]++;
        c0_z[(int)(X.z/dhz)]++;
    }

    for(int i=0;i<Types[1].X.size();i++){
        X=Types[1].X[i];
        c1_x[(int)(X.x/dhx)]++;
        c1_y[(int)(X.y/dhy)]++;
        c1_z[(int)(X.z/dhz)]++;
    }

    for(int i=0;i<Types[0].X.size();i++){
        for(int j=0;j<Types[0].X.size();j++){
            if(i==j)continue;
            X=Types[0].X[i];
            Y=Types[0].X[j];
            dx=X.x-Y.x;
            dy=X.y-Y.y;
            dz=X.z-Y.z;
            if(dx>Lx/2)dx-=Lx;
            if(dx<-Lx/2)dx+=Lx;
            if(dy>Ly/2)dy-=Ly;
            if(dy<-Ly/2)dy+=Ly;
            if(dz>Lz/2)dz-=Lz;
            if(dz<-Lz/2)dz+=Lz;
            r=sqrt(dx*dx+dy*dy+dz*dz);
            if(r/dr<BINNUM)g00[(int)(r/dr)]++;
        }
    }

    for(int i=0;i<Types[0].X.size();i++){
        for(int j=0;j<Types[1].X.size();j++){
            //if(i==j)continue;
            X=Types[0].X[i];
            Y=Types[1].X[j];
            dx=X.x-Y.x;
            dy=X.y-Y.y;
            dz=X.z-Y.z;
            if(dx>Lx/2)dx-=Lx;
            if(dx<-Lx/2)dx+=Lx;
            if(dy>Ly/2)dy-=Ly;
            if(dy<-Ly/2)dy+=Ly;
            if(dz>Lz/2)dz-=Lz;
            if(dz<-Lz/2)dz+=Lz;
            r=sqrt(dx*dx+dy*dy+dz*dz);
            if(r/dr<BINNUM)g01[(int)(r/dr)]++;
        }
    }

    for(int i=0;i<Types[1].X.size();i++){
        for(int j=0;j<Types[1].X.size();j++){
            if(i==j)continue;
            X=Types[1].X[i];
            Y=Types[1].X[j];
            dx=X.x-Y.x;
            dy=X.y-Y.y;
            dz=X.z-Y.z;
            if(dx>Lx/2)dx-=Lx;
            if(dx<-Lx/2)dx+=Lx;
            if(dy>Ly/2)dy-=Ly;
            if(dy<-Ly/2)dy+=Ly;
            if(dz>Lz/2)dz-=Lz;
            if(dz<-Lz/2)dz+=Lz;
            r=sqrt(dx*dx+dy*dy+dz*dz);
            if(r/dr<BINNUM)g11[(int)(r/dr)]++;
        }
    }
}
void Sample_Out(){
    cout<<"r g00 g01 g11"<<endl;
    for(int k=0;k<BINNUM;k++){
        double V=4*M_PI/3.0*(3*k*k+3*k+1);
        cout<<dr*(k+0.5)<<' '<<g00[k]/(V*sample_T)<<' '<<g01[k]/(V*sample_T)<<' '<<g11[k]/(V*sample_T)<<endl;
    }
    cout<<"c0_x c0_y c0_z"<<endl;
    for(int k=0;k<BINNUM;k++){
        cout<<c0_x[k]/sample_T<<' '<<c0_y[k]/sample_T<<' '<<c0_z[k]/sample_T<<endl;
    }

    cout<<"c1_x c1_y c1_z"<<endl;
    for(int k=0;k<BINNUM;k++){
        cout<<c1_x[k]/sample_T<<' '<<c1_y[k]/sample_T<<' '<<c1_z[k]/sample_T<<endl;
    }
}

int main(){
    Init();
    Global_Cell_List_Pointer->print();
    cout<<endl<<endl<<endl;
    sample_init();

//    double ifstat=false,ifprint=false;
//    for(int l=0;l<Cell_Veto_Lists.size();l++)Cell_Veto_Lists[l]->check_print();//For debug
//    for(int l=0;l<Cell_Veto_Lists.size();l++)Cell_Veto_Lists[l]->check_rate();//For debug

    for(int l=0;l<Cell_Veto_Lists.size();l++){
        Cell_Veto_Lists[l]->check_print();//For debug
        Cell_Veto_Lists[l]->update_max_num_per_cell();
        Cell_Veto_Lists[l]->check_print();//For debug
    }

    for(int k=0;k<20000*10;k++) {
        //if(k%100000==0)
        //Check();//For debug
//        cout<<"Moving along +x"<<endl;
        Monte_Carlo(33.33,1);
//        for(int l=0;l<Cell_Veto_Lists.size();l++)Cell_Veto_Lists[l]->check_print();//For debug

//        Global_Cell_List_Pointer->print();
        //Check();//For debug
//        cout<<"Moving along +y"<<endl;
        Monte_Carlo(33.33,2);
//        for(int l=0;l<Cell_Veto_Lists.size();l++)Cell_Veto_Lists[l]->check_print();//For debug

//        Global_Cell_List_Pointer->print();
//        Check();//For debug        
//        cout<<"Moving along +z"<<endl;
        Monte_Carlo(33.33,3);
//        for(int l=0;l<Cell_Veto_Lists.size();l++)Cell_Veto_Lists[l]->check_print();//For debug

//        Global_Cell_List_Pointer->print();
        
//        cout<<"Moving along -x"<<endl;
        Monte_Carlo(33.33,-1);
//        for(int l=0;l<Cell_Veto_Lists.size();l++)Cell_Veto_Lists[l]->check_print();//For debug

//        Global_Cell_List_Pointer->print();

//        cout<<"Moving along -y"<<endl;        
        Monte_Carlo(33.33,-2);
//        for(int l=0;l<Cell_Veto_Lists.size();l++)Cell_Veto_Lists[l]->check_print();//For debug

//        Global_Cell_List_Pointer->print();

//        cout<<"Moving along -z"<<endl;        
        Monte_Carlo(33.33,-3);
//        for(int l=0;l<Cell_Veto_Lists.size();l++)Cell_Veto_Lists[l]->check_print();//For debug

//        Global_Cell_List_Pointer->print();
        
//        cout<<"k="<<k<<" finished"<<endl;
        if((k%10000==0)||(k<10)){
            cout<<"k="<<k<<" finished"<<endl;
            //Global_Cell_List_Pointer->print();
            for(int l=0;l<Cell_Veto_Lists.size();l++){
                Cell_Veto_Lists[l]->update_max_num_per_cell();
                Cell_Veto_Lists[l]->check_print();//For debug
            }
        }
        if(k>100)sample();
    }
    Sample_Out();
}
