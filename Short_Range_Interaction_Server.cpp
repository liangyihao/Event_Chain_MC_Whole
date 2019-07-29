#include "Short_Range_Interaction_Server.hpp"
vector<Parameter_List> Parameter_List_For_Short_Range_Interaction;
Short_Range_Interaction_Between_Types::Short_Range_Interaction_Between_Types(int type_id_1, int type_id_2, vector<double4>*X_type_1, vector<double4>*X_type_2, double(*Event_Time_Generator)(double4,double4,int,double*,double), double*Params, double SHORT_INTERACTION_RANGE,bool Using_CellList_1,bool Using_CellList_2){
    this->type_id_1=type_id_1;
    this->type_id_2=type_id_2;
    this->X_type_1=X_type_1;
    this->X_type_2=X_type_2;
    this->Event_Time_Generator=Event_Time_Generator;
    this->Params=Params;
    this->SHORT_INTERACTION_RANGE=SHORT_INTERACTION_RANGE;
    this->Using_CellList_1=Using_CellList_1;
    this->Using_CellList_2=Using_CellList_2;
    double Lx,Ly,Lz;
    Lx=Params[0];
    Ly=Params[1];
    Lz=Params[2];
    Cell_List_Pointer_1=nullptr;
    Cell_List_Pointer_2=nullptr;
    if((Lx/SHORT_INTERACTION_RANGE<5)||(Ly/SHORT_INTERACTION_RANGE<5)||(Lz/SHORT_INTERACTION_RANGE<5)){
        this->Using_CellList_1=false;
        this->Using_CellList_2=false;
        cout<<"Cell List closed for this interaction"<<endl;
    }

    //double L_dim_mean=(Lx+Ly+Lz)/3.0;
    if(X_type_1->size()<9){
        this->Using_CellList_1=false;
        cout<<"Cell List closed for this interaction, type "<<type_id_1<<endl;
    }

    if(X_type_2->size()<9){
        this->Using_CellList_2=false;
        cout<<"Cell List closed for this interaction, type "<<type_id_2<<endl;
    }


    if(this->Using_CellList_1){
        Cell_List_Pointer_1=new CellList(Lx,Ly,Lz,SHORT_INTERACTION_RANGE,type_id_1,X_type_1,Event_Time_Generator,Params);
    }
    if(this->Using_CellList_2){
        if(type_id_1==type_id_2){
            this->Using_CellList_2=this->Using_CellList_1;
            Cell_List_Pointer_2=Cell_List_Pointer_1;
        }else{
            Cell_List_Pointer_2=new CellList(Lx,Ly,Lz,SHORT_INTERACTION_RANGE,type_id_2,X_type_2,Event_Time_Generator,Params);
        }
    }
}


void Short_Range_Interaction_Between_Types::Get_Event(TwoBody_Event&Event, int2 Active_Bead, double4 X_Active_Bead, int axis){
    if(Active_Bead.x==type_id_1){
        if(Using_CellList_2){
            Cell_List_Pointer_2->Get_Event(Event,Active_Bead,X_Active_Bead,axis);
        }else{
            //Get Event Directly
            double t;
            for(int l=0;l<X_type_2->size();l++){
                if((type_id_1==type_id_2)&&(Active_Bead.y==l))continue;
                t=Event_Time_Generator(X_Active_Bead,(*X_type_2)[l],axis,Params,Event.Event_Time);
                if(t<Event.Event_Time){
                    Event.Event_Time=t;
                    Event.Target_Bead.x=type_id_2;
                    Event.Target_Bead.y=l;
                }
            }
        }
    }else if(Active_Bead.x==type_id_2){
        if(Using_CellList_1){
            Cell_List_Pointer_1->Get_Event(Event,Active_Bead,X_Active_Bead,axis);
        }else{
            //Get Event Directly
            double t;
            for(int l=0;l<X_type_1->size();l++){
                if((type_id_1==type_id_2)&&(Active_Bead.y==l))continue;
                t=Event_Time_Generator(X_Active_Bead,(*X_type_1)[l],axis,Params,Event.Event_Time);
                if(t<Event.Event_Time){
                    Event.Event_Time=t;
                    Event.Target_Bead.x=type_id_1;
                    Event.Target_Bead.y=l;
                }
            }
        }
    }
}

void Short_Range_Interaction_Between_Types::Update(int2 ids, double4 New_X){
    if(ids.x==type_id_1){
        if(Using_CellList_1)Cell_List_Pointer_1->Update(ids,New_X);
    }else{
        if(Using_CellList_2)Cell_List_Pointer_2->Update(ids,New_X);
    }
}
