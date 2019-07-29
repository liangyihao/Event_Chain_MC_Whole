//deal with short-range interactions
//range:[0,Lx],[0,Ly],[0,Lz]
#ifndef SHORT_RANGE_INTERACTION_SERVER
#define SHORT_RANGE_INTERACTION_SERVER

#include <vector>
#include <iostream>
#include "public.hpp"
#include "CellList.hpp"

using namespace std;
extern vector<Parameter_List> Parameter_List_For_Short_Range_Interaction;
class Short_Range_Interaction_Between_Types{
	private:
	double(*Event_Time_Generator)(double4,double4,int,double*,double);
	double*Params;//Params[0]:Lx, Params[1]:Ly, Params[2]:Lz
	double SHORT_INTERACTION_RANGE;
	int type_id_1,type_id_2;
	vector<double4>*X_type_1;
	vector<double4>*X_type_2;

	/*For each type, create Cell List seperately.
	Sometimes, if the number of beads in such type is small, 
	users can choose turn-off the corresponding cell list.
	*/
	CellList* Cell_List_Pointer_1;
	CellList* Cell_List_Pointer_2;
	bool Using_CellList_1,Using_CellList_2;
	public:
	Short_Range_Interaction_Between_Types(int type_id_1, int type_id_2, vector<double4>*X_type_1, vector<double4>*X_type_2, double(*Event_Time_Generator)(double4,double4,int,double*,double), double*Params, double SHORT_INTERACTION_RANGE,bool Using_CellList_1,bool Using_CellList_2);
	void Get_Event(TwoBody_Event&Event, int2 Active_Bead, double4 X_Active_Bead, int axis);
	void Update(int2 ids, double4 New_X);
	
	double (*get_Gen())(double4,double4,int,double*,double){
		return Event_Time_Generator; 
	}
	double get_cutoff_r(){
		return SHORT_INTERACTION_RANGE;
	}
	int2 Interacting_Types(){
		int2 temp;
		temp.x=type_id_1;
		temp.y=type_id_2;
		return temp;
	}
};
#endif
