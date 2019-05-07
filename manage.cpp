#include "manage.hpp"
#include "basic.hpp"
#include <iostream>
/*Used to manage systematic variables*/
double Lx=10,Ly=10,Lz=10;//0<=x<Lx...


int Create_Type() {//Create Type, return it's id
    Bead_Type new_type;
    new_type.index=Types.size();
    new_type.Interactions_with_Beads.clear();
    new_type.Interactions_with_Types.clear();
    new_type.X.clear();
    Types.push_back(new_type);
    return new_type.index;
}

int Create_Bead(int type_id,double4 x){//Create Bead, return its secondary id. If failed, return -1.
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

    Types[type_id].X.push_back(x);
    vector<int3>temp;
    temp.clear();
    Types[type_id].Interactions_with_Beads.push_back(temp);
    return Types[type_id].X.size()-1;
}
