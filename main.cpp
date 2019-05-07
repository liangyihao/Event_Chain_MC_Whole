#include "basic.hpp"
#include "manage.hpp"
#include "Pairwise_Interaction.hpp"
#include "Random_Number.hpp"
#include "CellList.hpp"
#include <iostream>
#include <cmath>
#define BINNUM 1000
using namespace std;
vector<int2>end_to_end_ids;
const int Chain_Length=40;
const int Chain_Num=90;
const double d0=6.2;
const double K0=10;
const double l0=6.4;
void Check(){
    cout<<"checking"<<endl;
    double dx,dy,dz;
    for(int k=0;k<Types[0].X.size();k++)
        for(int l=0;l<Types[0].X.size();l++)
            if(k!=l){
                dx=Types[0].X[k].x-Types[0].X[l].x;
                dy=Types[0].X[k].y-Types[0].X[l].y;
                dz=Types[0].X[k].z-Types[0].X[l].z;
                if(dx>Lx/2)dx-=Lx;
                if(dx<-Lx/2)dx+=Lx;

                if(dy>Ly/2)dx-=Ly;
                if(dy<-Ly/2)dy+=Ly;

                if(dz>Lz/2)dz-=Lz;
                if(dz<-Lz/2)dz+=Lz;

                //cout<<"checking"<<k<<' '<<l<<endl;
                if(dx*dx+dy*dy+dz*dz<d0*d0/4.0){
                    cout<<"error"<<endl;
                    cout<<k<<' '<<l<<endl;
                    cout<<Types[0].X[k].x<<','<<Types[0].X[k].y<<','<<Types[0].X[k].z<<endl;
                    cout<<Types[0].X[l].x<<','<<Types[0].X[l].y<<','<<Types[0].X[l].z<<endl;
                    cout<<sqrt(dx*dx+dy*dy+dz*dz)<<endl;
                    exit(0);
                }
            }
    cout<<"Checking finished"<<endl;
}

void Init(){
    end_to_end_ids.clear();
    Lx=100;Ly=100;Lz=100;
    rand_init(1);
    int tid1;
    tid1=Create_Type();//type A, hard sphere diameter 10
    Create_Hard_Sphere_Interaction_Between_Types(tid1,tid1,d0);

    double4 X;
    double xc,yc,zc;

    int2 id,id_o;
    id.x=tid1;

    int2 end_end_Pair;

    int chain_length=0,chain_num=0;
    for(xc=0;xc<Lx-l0;xc+=l0)
        for(yc=0;yc<Ly-l0;yc+=l0)
            for(zc=0;zc<Lz-l0;zc+=l0) {
                if(chain_num==Chain_Num)break;

                id_o=id;
                X.x=xc;X.y=yc;X.z=zc;
                id.y=Create_Bead(tid1,X);
                chain_length++;
                if(chain_length==1)end_end_Pair.x=id.y;

                if(chain_length>=2)Create_Spring_Interaction_Between_Beads(id_o,id,K0,l0);
                
                if(chain_length==Chain_Length){
                    end_end_Pair.y=id.y;
                    end_to_end_ids.push_back(end_end_Pair);
                    chain_length=0;
                    chain_num++;
                }
            }
    
    for(int k=0;k<end_to_end_ids.size();k++)cout<<end_to_end_ids[k].x<<' '<<end_to_end_ids[k].y<<endl;
    
    Global_Cell_List_Pointer=new CellList(Lx, Ly, Lz, MAX_SHORT_INTERACTION_RANGE, (&Types));
}

double g[BINNUM];
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
}

int main(){
    Init();
    Global_Cell_List_Pointer->print();
    cout<<endl<<endl<<endl;
    double ifstat=false,ifprint=false; 
    for(int k=0;k<2000000;k++) {
        //if(k%100000==0)
        //Check();//For debug
//        cout<<"Moving along +x"<<endl;
        Monte_Carlo(800,1);
//        Global_Cell_List_Pointer->print();
        //Check();//For debug
//        cout<<"Moving along +y"<<endl;
        Monte_Carlo(800,2);
//        Global_Cell_List_Pointer->print();
//        Check();//For debug        
//        cout<<"Moving along +z"<<endl;
        Monte_Carlo(800,3);
//        Global_Cell_List_Pointer->print();
        
//        cout<<"Moving along -x"<<endl;
        Monte_Carlo(800,-1);
//        Global_Cell_List_Pointer->print();

//        cout<<"Moving along -y"<<endl;        
        Monte_Carlo(800,-2);
//        Global_Cell_List_Pointer->print();

//        cout<<"Moving along -z"<<endl;        
        Monte_Carlo(800,-3);
//        Global_Cell_List_Pointer->print();
        ifprint=false;
        
//        cout<<"k="<<k<<" finished"<<endl;
        if(k%10000==0){
            cout<<"k="<<k<<" finished"<<endl;
            //Global_Cell_List_Pointer->print();
            ifprint=true;
        }
        if(k>300000)ifstat=true;
        if(k%100==0)sample(ifstat,ifprint);
    }

    double dr=Ly/BINNUM;;
    for(int k=0;k<BINNUM;k++) {
        double r1,r2;
        r1=dr*k;
        r2=dr*(k+1);
        cout<<dr*k<<" "<<g[k]/(4*M_PI/3*(r2*r2*r2-r1*r1*r1)*sample_T)<<endl;
        //cout<<dr*(k-BINNUM/2)<<" "<<g[k]<<endl;
    }
}
