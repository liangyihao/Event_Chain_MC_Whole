/*
Author: Yihao Liang
liangyihaosjtu@gmail.com
This code is for Event Chain Monte Carlo for pairwise interacting many body system
*/
#include "manage.hpp"
#include "basic.hpp"
#include <iostream>
#include "Input_File_Parser.hpp"
#include "dcd_writer.hpp"
#include "CellVetoList.hpp"
#include <cmath>
/*Used to manage systematic variables*/
double Lx=10,Ly=10,Lz=10;//0<=x<Lx...
double Bjerrum_Length=7.117;
vector<Instruction>Instruction_list;
int loop_times;

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

void Run(char*InputFileName){
    rand_init(1);
    Input_File_Parser(InputFileName);
    Hard_Repulsion_Checker();
    initialize_sys_charge();
    for(int l=0;l<Cell_Veto_Lists.size();l++){
        Cell_Veto_Lists[l]->check_print();//For debug
        Cell_Veto_Lists[l]->update_max_num_per_cell();
        Cell_Veto_Lists[l]->check_print();//For debug
    }
    Output_DCD_init(InputFileName);
    for(int l=0;l<loop_times;l++){
        for(int k=0;k<Instruction_list.size();k++){
            if(Instruction_list[k].Command==0){//Do ECMC
                Monte_Carlo(Instruction_list[k].Double_Para[0],Instruction_list[k].Int_Para[0]);
            }else if(Instruction_list[k].Command==1){//Do output
                if(l%Instruction_list[k].Int_Para[1]==0)cout<<l<<endl;
                if(l<Instruction_list[k].Int_Para[0])continue;
                if(l%Instruction_list[k].Int_Para[1]==0){
                    Output_DCD();
                    next_input_file_writer(InputFileName);
                }
            }else if(Instruction_list[k].Command==2){//Do reconstruct of CellVeto List
                if(l%Instruction_list[k].Int_Para[0]!=0)continue;
                for(int zz=0;zz<Cell_Veto_Lists.size();zz++){
                    Cell_Veto_Lists[zz]->Reconstruct_CellVeto_List(Instruction_list[k].Int_Para[1]);
                    Cell_Veto_Lists[zz]->update_max_num_per_cell();
                }            
            }else if(Instruction_list[k].Command==3){//Do refresh of CellVeto List
                if(l%Instruction_list[k].Int_Para[0]!=0)continue;
                for(int zz=0;zz<Cell_Veto_Lists.size();zz++)Cell_Veto_Lists[zz]->update_max_num_per_cell();
            }else if(Instruction_list[k].Command==4){//Do check of CellVeto List
                if(l%Instruction_list[k].Int_Para[0]!=0)continue;
                for(int zz=0;zz<Cell_Veto_Lists.size();zz++)Cell_Veto_Lists[zz]->check_print();
            }else if(Instruction_list[k].Command==5){//Do compress the system
                //Add code here
            }
        }
    }
    Output_DCD_Close();
    next_input_file_writer(InputFileName);
}
