/*
Author: Yihao Liang
liangyihaosjtu@gmail.com
This code is for Event Chain Monte Carlo for pairwise interacting many body system
*/
#ifndef PUBLIC
#define PUBLIC
#define EPSILON (1E-14)
#include <vector>
using namespace std;

typedef struct{
    double x,y,z,w;
}double4;

typedef struct{
    int x,y;
}int2;

typedef struct{
    int x,y,z;
}int3;

typedef struct{
    double Event_Time;
    int2 Target_Bead;
}TwoBody_Event;
/*
We classify beads into types.
The data structure "Bead_Type" has all information of current bead type
And these bead type informations will be organized into arrays.
When we refer a bead, it has two ids, the type-id and the bead-id in that type
The id of a bead is usually represented in an int2 data. id.x: type-id,   id.y: bead-id in that type
*/
typedef struct{
    int index;//the id of this bead type

    vector<double4>X;//size=N
    //the coordinates of beads in this type. X[i] stores information of i'th bead in current bead type.
    //X has 4 field: x,y,z,w. (x,y,z) is position, w is valence of charge.

    vector<vector<int3> >Interactions_with_Beads;//size=N
    //Interactions_with_Beads[i] stores information on who else has special interaction with bead i in current bead type
    //It is a list of special interactions for each bead
    //Interactions_with_Beads[i][..].x    the type-id of the bead who interact with bead i in the current bead type
    //Interactions_with_Beads[i][..].y    the secondary bead-id of the bead who interact with bead i in the current bead type
    //Interactions_with_Beads[i][..].z    interaction id in the Event_Time_Generator_List_For_Bonds list

    vector<int>Interactions_with_Types;//size=N
    //Interactions_with_Types[i] stores information of the type-type interaction the current bead type interact with
    //Interactions_with_Types[i] the short-range-interaction-id in the short range interaction server list
}Bead_Type;

class Parameter_List{
    public:
    double*data;
    //data[0]: Lx
    //data[1]: Ly
    //data[2]: Lz
    ///...
    //When calling function for Event Time, just pass the address of data[0]
    Parameter_List(double Lx,double Ly,double Lz) {
        data=new double [8];
        data[0]=Lx;
        data[1]=Ly;
        data[2]=Lz;
    }
};

typedef struct{
  int Command;//0:ECMC 1:Output
  vector<double>Double_Para;
  vector<int>Int_Para;
}Instruction;

#endif
