/*
Author: Yihao Liang
liangyihaosjtu@gmail.com
This code is for Event Chain Monte Carlo for pairwise interacting many body system
*/
#include "manage.hpp"
#include "basic.hpp"
#include <iostream>
#include "CellVetoList.hpp"
#include <cmath>
/*Used to manage systematic variables*/
double Lx=10,Ly=10,Lz=10;//0<=x<Lx...

vector<double>valence_of_type;
vector<double>valence_list;//record all the non-zero valence of beads
vector<CellVetoList*> Cell_Veto_Lists;//for all valences, you need to create a Cell-Veto List

int Create_Type() {//Create Type, return it's id
    Bead_Type new_type;
    new_type.index=Types.size();
    new_type.Interactions_with_Beads.clear();
    new_type.Interactions_with_Types.clear();
    new_type.X.clear();
    Types.push_back(new_type);
    valence_of_type.push_back(0);
    return new_type.index;
}

int Create_Bead(int type_id,double4 x){//Create Bead, return its secondary id. If failed, return -1.
    if(abs(valence_of_type[type_id])>EPSILON){
        cout<<"Error Creating Bead, should call Create_Charged_Bead"<<endl;
        exit(0);
    }
    if(type_id>=Types.size()) {
        cout<<"Error Creating Bead: type-id incorrect"<<endl;
        return -1;
    }
    
    if((x.x>Lx)||(x.x<0)){
        cout<<"Error Creating Bead: position out of boundary"<<endl;
        return -1;
    }

    if((x.y>Ly)||(x.y<0)){
        cout<<"Error Creating Bead: position out of boundary"<<endl;
        return -1;
    }

    if((x.z>Lz)||(x.z<0)){
        cout<<"Error Creating Bead: position out of boundary"<<endl;
        return -1;
    }
    x.w=0;//w field is reserved for charged particle
    Types[type_id].X.push_back(x);
    vector<int3>temp;
    temp.clear();
    Types[type_id].Interactions_with_Beads.push_back(temp);
    return Types[type_id].X.size()-1;
}

int Create_Charged_Type(double valence) {//Create Charged Type, return it's id
    if(abs(valence)<EPSILON){
        cout<<"Error: Charged type cannot have valence 0"<<endl;
        exit(0);
    }
    Bead_Type new_type;
    new_type.index=Types.size();
    new_type.Interactions_with_Beads.clear();
    new_type.Interactions_with_Types.clear();
    new_type.X.clear();
    Types.push_back(new_type);
    valence_of_type.push_back(valence);
    
    bool is_new_valence=true;
    for(int k=0;k<valence_list.size();k++)
        if(abs(valence_list[k]-valence)<EPSILON){
            is_new_valence=false;
            break;
        }
    if(is_new_valence)valence_list.push_back(valence);
    return new_type.index;
}

int Create_Charged_Bead(int type_id,double4 x){//Create Charged Bead, return its secondary id. If failed, return -1.
    if(type_id>=Types.size()) {
        cout<<"Error Creating Bead: type-id incorrect"<<endl;
        return -1;
    }
    
    if((x.x>Lx)||(x.x<0)){
        cout<<"Error Creating Bead: position out of boundary"<<endl;
        return -1;
    }

    if((x.y>Ly)||(x.y<0)){
        cout<<"Error Creating Bead: position out of boundary"<<endl;
        return -1;
    }

    if((x.z>Lz)||(x.z<0)){
        cout<<"Error Creating Bead: position out of boundary"<<endl;
        return -1;
    }
    x.w=valence_of_type[type_id];
    Types[type_id].X.push_back(x);
    vector<int3>temp;
    temp.clear();
    Types[type_id].Interactions_with_Beads.push_back(temp);
    return Types[type_id].X.size()-1;
}

void initialize_sys_charge(){
    CellVetoList*temp;
    double valence;
    for(int k=0;k<valence_list.size();k++){
        valence=valence_list[k];
        temp=new CellVetoList(Lx,Ly,Lz,valence,(&Types));
        Cell_Veto_Lists.push_back(temp);
    }
}
