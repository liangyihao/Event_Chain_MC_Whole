#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "basic.hpp"
#include "manage.hpp"
#include "Short_Range_potentials_Basic.hpp"
#include "Short_Range_potentials_extension.hpp"
#include "Input_File_Parser.hpp"
using namespace std;
typedef struct{
  string name;
  bool ischarged;
  double diameter;
}Type_Definition;
vector<Type_Definition>Type_Definition_List;
// For each type, it's id in Type_Definition_List equals it's id in Types

typedef struct{
  string name;
  string Interaction_Name;
  vector<double>para;
}Bond_Definition;
vector<Bond_Definition>Bond_Definition_List;

typedef struct{
    int Bond_Def_Id;
    int2 bead1,bead2;
}Bond;
vector<Bond>Bonds_List;
vector<int2>IO_Order_of_Beads;//Record the order of beads in the input file, output in this order
vector<vector<int> >IO_ids;//IO_ids[typeid][beadid]show you the id in input file
int Input_File_Parser(const char* fileName)
{  
  ifstream ifs(fileName);
  if (!ifs.good()) {
    cout << "Cannot read file " << fileName << std::endl;
    exit(0);
  }

  bool beads_pos_readed = 0;

  string s;
  int line_num=0;
  Type_Definition Type_Definition_Temp;
  Bond_Definition Bond_Definition_Temp;
  Instruction Instruction_Temp;

  while (!ifs.eof()) {
    s.clear();
    getline(ifs,s);line_num++;//when reading each line, ignore comment
    cout<<endl<<"Line"<<line_num<<':'<<s<<endl;
    int pos=s.find("%");
    if(pos!=s.npos){
      int len=s.length();
      s.erase(pos,len-pos);
      cout<<"Refined line: "<<s<<endl;
    }

    istringstream iss(s);
    string word;

    while (iss >> word) {
      if (word=="Lx") {iss>>Lx;continue;}
      if (word=="Ly") {iss>>Ly;continue;}
      if (word=="Lz") {iss>>Lz;continue;}
      if (word=="Bjerrum-Length") {iss>>Bjerrum_Length;continue;}
      if (word=="Type") {//Command format:   Type TypeName
        iss>> word;
        Type_Definition_Temp.name=word;
        Type_Definition_Temp.ischarged=0;
        Type_Definition_Temp.diameter=0;
        iss>> word;
        if(word=="-Diameter")iss>>Type_Definition_Temp.diameter;
        Type_Definition_List.push_back(Type_Definition_Temp);
        Create_Type();
        continue;
      }
      if (word=="ChargedType"){//ChargedType Ca +1 -Diameter 6 %Define Charged Type Ca with valence +1
        iss>> word;
        Type_Definition_Temp.name=word;
        Type_Definition_Temp.ischarged=1;
        double valence;
        iss>>valence;
        Type_Definition_Temp.diameter=0;
        iss>> word;
        if(word=="-Diameter")iss>>Type_Definition_Temp.diameter;
        Type_Definition_List.push_back(Type_Definition_Temp);
        Create_Charged_Type(valence);
        continue;
      }
      if (word=="Interaction-Type-Type") {
        //Command format: Interaction-Type-Type TypeName TypeName InteractionName Parameters
        int TypeId1,TypeId2;
        bool Using_CellList1=true,Using_CellList2=true;
        iss>> word;//TypeName
        for(TypeId1=0;TypeId1<Type_Definition_List.size();TypeId1++){
          if(word==Type_Definition_List[TypeId1].name)break;
        }
        if(TypeId1==Type_Definition_List.size()){
          cout<<"Error at line "<<line_num<<" there are undifined type"<<endl;
          exit(0);
        }
        iss>> word;//TypeName
        for(TypeId2=0;TypeId2<Type_Definition_List.size();TypeId2++){
          if(word==Type_Definition_List[TypeId2].name)break;
        }
        if(TypeId2==Type_Definition_List.size()){
          cout<<"Error at line "<<line_num<<" there are undifined type"<<endl;
        }
        iss>>word;
        if(word=="Hard"){
          double d;
          iss>>d;
          Create_Hard_Sphere_Interaction_Between_Types(TypeId1,TypeId2,Using_CellList1,Using_CellList2,d);
        }else if(word=="LJ"){
          double sigma,epsilon,rcut;
          iss>>sigma>>epsilon>>rcut;
          Create_LJ_Interaction_Between_Types(TypeId1,TypeId2,Using_CellList1,Using_CellList2,sigma,epsilon,rcut); 
        }else if(word=="Gauss"){
          double sigma,epsilon,rcut;
          iss>>sigma>>epsilon>>rcut;
          Create_Gauss_Interaction_Between_Types(TypeId1,TypeId2,Using_CellList1,Using_CellList2,sigma,epsilon,rcut);
        }else {
          cout<<word<<" isn't supported in this version"<<endl;
          exit(0);
        }
        continue;
      }

      if (word=="Bond-Type") {//Format: Bond-Type SP1 Spring 10.0 6.0
        iss>>(Bond_Definition_Temp.name);
        iss>>(Bond_Definition_Temp.Interaction_Name);
        double temp;
        Bond_Definition_Temp.para.clear();
        while(!iss.eof()){
          iss>>temp;
          Bond_Definition_Temp.para.push_back(temp);
        }
        if((Bond_Definition_Temp.Interaction_Name=="Spring")&&(Bond_Definition_Temp.para.size()!=2)){
          cout<<"Error at line "<<line_num<<" parameters incorrect"<<endl;
          exit(0);
        }
        Bond_Definition_List.push_back(Bond_Definition_Temp);
        continue;
      }

      if (word.find("Positions")==0) {
        IO_Order_of_Beads.clear();
        IO_ids.resize(Types.size());for(int k=0;k<Types.size();k++)IO_ids[k].clear();
        while(1) {
          getline(ifs,s);line_num++;//when reading each line, ignore comment
          cout<<endl<<"Line"<<line_num<<':'<<s<<endl;
          int pos=s.find("%");
          if(pos!=s.npos){
            int len=s.length();
            s.erase(pos,len-pos);
            cout<<"Refined line: "<<s<<endl;
          }

          istringstream iss_pos(s);
          string word_pos;
          double4 X;
          iss_pos>>word_pos;
          if(word_pos.find("End-Positions")==0)break;
          iss_pos>>X.x>>X.y>>X.z;
          X.w=0;
          int typeid_temp;
          for(typeid_temp=0;typeid_temp<Type_Definition_List.size();typeid_temp++){//search it's type id
            if(Type_Definition_List[typeid_temp].name==word_pos)break;
          }
          if(typeid_temp==Type_Definition_List.size()){
            cout<<"Error at line "<<line_num<<" No such type"<<endl;
            exit(0);
          }

          int beadid_temp;
          if(Type_Definition_List[typeid_temp].ischarged==0)beadid_temp=Create_Bead(typeid_temp,X);else beadid_temp=Create_Charged_Bead(typeid_temp,X);
          int2 wholeid_temp;
          wholeid_temp.x=typeid_temp;
          wholeid_temp.y=beadid_temp;
          IO_Order_of_Beads.push_back(wholeid_temp);
          IO_ids[typeid_temp].push_back(IO_Order_of_Beads.size()-1);
        }
        continue;
      }

      if (word.find("Bonds")==0) {//Format: Bond-Name TypeName1 BeadId1 TypeName2 BeadId2
        while(1) {
          getline(ifs,s);line_num++;//when reading each line, ignore comment
          cout<<endl<<"Line"<<line_num<<':'<<s<<endl;
          int pos=s.find("%");
          if(pos!=s.npos){
            int len=s.length();
            s.erase(pos,len-pos);
            cout<<"Refined line: "<<s<<endl;
          }

          istringstream iss_bon(s);
          string word_bon,TypeName1,TypeName2;
          int2 bid1,bid2;
          iss_bon>>word_bon;
          if(word_bon.find("End-Bonds")==0)break;
          iss_bon>>TypeName1>>bid1.y>>TypeName2>>bid2.y;

          for(bid1.x=0;bid1.x<Type_Definition_List.size();bid1.x++){//search it's type id
            if(Type_Definition_List[bid1.x].name==TypeName1)break;
          }
          if(bid1.x==Type_Definition_List.size()){
            cout<<"Error at line "<<line_num<<" No such type"<<endl;
            exit(0);
          }
          if(bid1.y>=Types[bid1.x].X.size()){
            cout<<"Error at line "<<line_num<<" No such bead"<<endl;
            exit(0);
          }
          for(bid2.x=0;bid2.x<Type_Definition_List.size();bid2.x++){//search it's type id
            if(Type_Definition_List[bid2.x].name==TypeName2)break;
          }
          if(bid2.x==Type_Definition_List.size()){
            cout<<"Error at line "<<line_num<<" No such type"<<endl;
            exit(0);
          }
          if(bid2.y>=Types[bid2.x].X.size()){
            cout<<"Error at line "<<line_num<<" No such bead"<<endl;
            exit(0);
          }

          int bond_def_id;
          for(bond_def_id=0;bond_def_id<Bond_Definition_List.size();bond_def_id++){
            if(Bond_Definition_List[bond_def_id].name==word_bon)break;
          }
          if(bond_def_id==Bond_Definition_List.size()){
            cout<<"Error at line "<<line_num<<" Bond undefined"<<endl;
            exit(0);
          }
//          cout<<bid1.x<<' '<<bid1.y<<' '<<bid2.x<<' '<<bid2.y<<endl;
          if(Bond_Definition_List[bond_def_id].Interaction_Name=="Spring"){
            Create_Spring_Interaction_Between_Beads(bid1,bid2,Bond_Definition_List[bond_def_id].para[0],Bond_Definition_List[bond_def_id].para[1]);
          }
          Bond Bond_Temp;
          Bond_Temp.Bond_Def_Id=bond_def_id;
          Bond_Temp.bead1=bid1;
          Bond_Temp.bead2=bid2;
          Bonds_List.push_back(Bond_Temp);
        }
        continue;
      }

      if (word=="loop") {
        iss>>loop_times;
        while(1) {
          getline(ifs,s);line_num++;//when reading each line, ignore comment
          cout<<endl<<"Line"<<line_num<<':'<<s<<endl;
          int pos=s.find("%");
          if(pos!=s.npos){
            int len=s.length();
            s.erase(pos,len-pos);
            cout<<"Refined line: "<<s<<endl;
          }

          string command;
          istringstream com_iss(s);
          com_iss>>command;
          if(command=="End-loop")break;
          if(command=="")continue;
          if(command=="ECMC"){
            double L;
            int axis;
            string axis_string;
            Instruction_Temp.Command=0;
            Instruction_Temp.Double_Para.clear();
            Instruction_Temp.Int_Para.clear();
            com_iss>>axis_string;
            com_iss>>L;
            if(axis_string=="+x"){
              Instruction_Temp.Int_Para.push_back(+1);
            }else if(axis_string=="-x"){
              Instruction_Temp.Int_Para.push_back(-1);
            }else if(axis_string=="+y"){
              Instruction_Temp.Int_Para.push_back(+2);
            }else if(axis_string=="-y"){
              Instruction_Temp.Int_Para.push_back(-2);
            }else if(axis_string=="+z"){
              Instruction_Temp.Int_Para.push_back(+3);
            }else if(axis_string=="-z"){
              Instruction_Temp.Int_Para.push_back(-3);
            }else{
              cout<<"Error at line "<<line_num<<" axis name incorrect"<<endl;
              exit(0);
            }
            Instruction_Temp.Double_Para.push_back(L);
            Instruction_list.push_back(Instruction_Temp);
          }else if(command=="Out"){
            Instruction_Temp.Command=1;
            Instruction_Temp.Double_Para.clear();
            Instruction_Temp.Int_Para.clear();
            int MinIt,PerNum;
            com_iss>>MinIt>>PerNum;
            Instruction_Temp.Int_Para.push_back(MinIt);
            Instruction_Temp.Int_Para.push_back(PerNum);
            Instruction_list.push_back(Instruction_Temp);
          }else if(command=="Reconstruct-CellVeto-List"){//Do reconstruct of CellVeto List
            Instruction_Temp.Command=2;
            Instruction_Temp.Double_Para.clear();
            Instruction_Temp.Int_Para.clear();
            int PerNum,MaxNPerCell;
            com_iss>>PerNum>>MaxNPerCell;
            Instruction_Temp.Int_Para.push_back(PerNum);
            Instruction_Temp.Int_Para.push_back(MaxNPerCell);
            Instruction_list.push_back(Instruction_Temp);          
          }else if(command=="Refresh-CellVeto-maxnumpercell"){//Do refresh of CellVeto List
            Instruction_Temp.Command=3;
            Instruction_Temp.Double_Para.clear();
            Instruction_Temp.Int_Para.clear();
            int PerNum;
            com_iss>>PerNum;
            Instruction_Temp.Int_Para.push_back(PerNum);
            Instruction_list.push_back(Instruction_Temp);                    
          }else if(command=="Check-CellVeto-List"){//Do check of CellVeto List
            Instruction_Temp.Command=4;
            Instruction_Temp.Double_Para.clear();
            Instruction_Temp.Int_Para.clear();
            int PerNum;
            com_iss>>PerNum;
            Instruction_Temp.Int_Para.push_back(PerNum);
            Instruction_list.push_back(Instruction_Temp);                    
          }else{
            cout<<"Error at line "<<line_num<<" No such instruction"<<endl;
            exit(0);
          }
        }
        continue;
      }

      cout<<"Error at line "<<line_num<<", undefined command"<<endl;
    }
  }

  ifs.close();
  //
  char XML_FILE_NAME[256];
  int i;
  int last_dot_pos=-1;
  for(i=0;fileName[i]!=0;i++){
    if(fileName[i]=='.')last_dot_pos=i;
    XML_FILE_NAME[i]=fileName[i];
  }
  if(last_dot_pos==-1)last_dot_pos=i;
  XML_FILE_NAME[last_dot_pos]='.';
  XML_FILE_NAME[last_dot_pos+1]='x';
  XML_FILE_NAME[last_dot_pos+2]='m';
  XML_FILE_NAME[last_dot_pos+3]='l';
  XML_FILE_NAME[last_dot_pos+4]=0;
  //
  xml_write(XML_FILE_NAME,0);
  return 1;
}


void xml_write(const char* fileName, const int timestep)
{
  ofstream ofs;
  ofs.open(fileName);
  int N_Beads=0;
  for(int k=0;k<Types.size();k++)N_Beads+=Types[k].X.size();
  ofs << "<\?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
  ofs << "<hoomd_xml version=\"1.5\">\n";
  ofs << "<configuration time_step=\"" << timestep;
  ofs << "\" dimensions=\"3\" natoms=\"" << N_Beads << "\" >\n";
  ofs << "<box lx=\"" << Lx << "\" ly=\"" << Ly << "\" lz=\"" << Lz;
  ofs << "\" xy=\"0\" xz=\"0\" yz=\"0\"/>\n";



  ofs << "<position num=\"" << N_Beads << "\">\n";
  ofs.precision(15);//Output position
  /*
  for(int type_id=0;type_id<Types.size();type_id++){
    for(int bead_id=0;bead_id<Types[type_id].X.size();bead_id++){
      ofs<<Types[type_id].X[bead_id].x<<' '<<Types[type_id].X[bead_id].y<<' '<<Types[type_id].X[bead_id].z<<endl;
    }
  }*/
  int type_id,bead_id;
  for(int k=0;k<IO_Order_of_Beads.size();k++){
    type_id=IO_Order_of_Beads[k].x;
    bead_id=IO_Order_of_Beads[k].y;
    ofs<<Types[type_id].X[bead_id].x<<' '<<Types[type_id].X[bead_id].y<<' '<<Types[type_id].X[bead_id].z<<endl;
  }
  ofs << "</position>"<<endl;

  ofs <<"<diameter num=\""<<N_Beads<<"\">"<<endl;//Output type of each bead
  /*
  for(int type_id=0;type_id<Types.size();type_id++){
    for(int bead_id=0;bead_id<Types[type_id].X.size();bead_id++){
      ofs<<Type_Definition_List[type_id].diameter<<endl;
    }
  }*/
  for(int k=0;k<IO_Order_of_Beads.size();k++){
    type_id=IO_Order_of_Beads[k].x;
    bead_id=IO_Order_of_Beads[k].y;
    ofs<<Type_Definition_List[type_id].diameter<<endl;
  }
  ofs << "</diameter>"<<endl;

  ofs <<"<type>"<<endl;//Output type of each bead
  /*
  for(int type_id=0;type_id<Types.size();type_id++){
    for(int bead_id=0;bead_id<Types[type_id].X.size();bead_id++){
      ofs<<Type_Definition_List[type_id].name<<endl;
    }
  }*/
  for(int k=0;k<IO_Order_of_Beads.size();k++){
    type_id=IO_Order_of_Beads[k].x;
    bead_id=IO_Order_of_Beads[k].y;
    ofs<<Type_Definition_List[type_id].name<<endl;
  }
  ofs << "</type>"<<endl;
  
  ofs << "<charge>"<<endl;//Output Charge of each bead
  /*
  for(int type_id=0;type_id<Types.size();type_id++){
    for(int bead_id=0;bead_id<Types[type_id].X.size();bead_id++){
      ofs<<valence_of_type[type_id]<<endl;
    }
  }*/
  for(int k=0;k<IO_Order_of_Beads.size();k++){
    type_id=IO_Order_of_Beads[k].x;
    bead_id=IO_Order_of_Beads[k].y;
    ofs<<valence_of_type[type_id]<<endl;
  }
  ofs << "</charge>"<<endl;
  
  ofs << "<bond num=\""<<Bonds_List.size()<<"\">\n";
  for(int k=0;k<Bonds_List.size();k++) {
    int2 bid1,bid2;
    int gid1,gid2;
    bid1=Bonds_List[k].bead1;
    bid2=Bonds_List[k].bead2;
    /*
    gid1=0;
    for(int l=0;l<bid1.x;l++)gid1+=Types[l].X.size();
    gid1+=bid1.y;

    gid2=0;
    for(int l=0;l<bid2.x;l++)gid2+=Types[l].X.size();
    gid2+=bid2.y;
    */
    gid1=IO_ids[bid1.x][bid1.y];
    gid1=IO_ids[bid2.x][bid2.y];
    ofs<<Bond_Definition_List[Bonds_List[k].Bond_Def_Id].name<<' '<<gid1<<' '<<gid2<<endl;
  }
  ofs << "</bond>\n";

  ofs << "</configuration>\n";
  ofs << "</hoomd_xml>\n";

  ofs.close();
}

void Hard_Repulsion_Checker(){
  for(int interaction_id=0;interaction_id<Short_Range_Interaction_Between_Types_List.size();interaction_id++){
    if(Short_Range_Interaction_Between_Types_List[interaction_id]->get_Gen()!=Event_Time_Hard_Sphere)continue;

    int2 type_ids;
    int type_id1,type_id2;
    type_ids=Short_Range_Interaction_Between_Types_List[interaction_id]->Interacting_Types();
    type_id1=type_ids.x;
    type_id2=type_ids.y;

    double d;
    d=Short_Range_Interaction_Between_Types_List[interaction_id]->get_cutoff_r();
    //Now start to check
    int bead_id1,bead_id2;
    for(bead_id1=0;bead_id1<Types[type_id1].X.size();bead_id1++)
      for(bead_id2=0;bead_id2<Types[type_id2].X.size();bead_id2++){
        if((type_id1==type_id2)&&(bead_id1==bead_id2))continue;
        double dx,dy,dz,dr;
        dx=Types[type_id1].X[bead_id1].x-Types[type_id2].X[bead_id2].x;
        dy=Types[type_id1].X[bead_id1].y-Types[type_id2].X[bead_id2].y;
        dz=Types[type_id1].X[bead_id1].z-Types[type_id2].X[bead_id2].z;
        while(dx<-Lx/2)dx+=Lx;
        while(dx>+Lx/2)dx-=Lx;
        while(dy<-Ly/2)dy+=Ly;
        while(dy>+Ly/2)dy-=Ly;
        while(dz<-Lz/2)dz+=Lz;
        while(dz>+Lz/2)dz-=Lz;
        dr=sqrt(dx*dx+dy*dy+dz*dz);
        if(dr<d*(1-EPSILON)){
          cout<<"Error, overlap between"<<endl;
          cout<<Type_Definition_List[type_id1].name<<'-'<<bead_id1<<" and "<<Type_Definition_List[type_id2].name<<'-'<<bead_id2<<endl;
          cout<<"distance is "<<dr<<endl;
          cout<<"but it should be greater than "<<d<<endl;
          exit(0);
        }
      }
  }
}

void next_input_file_writer(const char* fileName){
  ifstream ifs(fileName);

  //
  char NEXT_FILE_NAME[256];
  int i;
  int last_dot_pos=-1;
  for(i=0;fileName[i]!=0;i++){
    if(fileName[i]=='.')last_dot_pos=i;
    NEXT_FILE_NAME[i]=fileName[i];
  }
  if(last_dot_pos==-1)last_dot_pos=i;
  NEXT_FILE_NAME[last_dot_pos]='_';
  NEXT_FILE_NAME[last_dot_pos+1]='n';
  NEXT_FILE_NAME[last_dot_pos+2]='e';
  NEXT_FILE_NAME[last_dot_pos+3]='x';
  NEXT_FILE_NAME[last_dot_pos+4]='t';
  NEXT_FILE_NAME[last_dot_pos+5]='.';
  NEXT_FILE_NAME[last_dot_pos+6]='e';
  NEXT_FILE_NAME[last_dot_pos+7]='c';
  NEXT_FILE_NAME[last_dot_pos+8]='m';
  NEXT_FILE_NAME[last_dot_pos+9]='c';
  NEXT_FILE_NAME[last_dot_pos+10]=0;
  ofstream ofs(NEXT_FILE_NAME);
  //

  string s_in,s_proc;
  int line_num=0;

  while (!ifs.eof()) {
    s_in.clear();
    getline(ifs,s_in);line_num++;//when reading each line, ignore comment
    s_proc=s_in;
    int pos=s_in.find("%");
    if(pos!=s_in.npos){
      int len=s_in.length();
      s_proc.erase(pos,len-pos);
    }

    istringstream iss(s_proc);
    string word;

    iss >> word;
    if (word.find("Positions")==0) {
        ofs<<s_in<<endl;
        int type_id,bead_id;
        type_id=0;bead_id=-1;
        while(1) {
          getline(ifs,s_in);line_num++;//when reading each line, ignore comment
          s_proc=s_in;
          int pos=s_in.find("%");
          if(pos!=s_in.npos){
            int len=s_in.length();
            s_proc.erase(pos,len-pos);
          }

          istringstream iss_pos(s_proc);
          string word_pos;
          iss_pos>>word_pos;
          if(word_pos.find("End-Positions")==0){
            ofs<<s_in<<endl;
            break;
          }
          bead_id++;
          if(bead_id==Types[type_id].X.size()){bead_id=0;type_id++;}
          ofs<<fixed;
          ofs<<setprecision(15)<<Type_Definition_List[type_id].name<<' '<<Types[type_id].X[bead_id].x<<' '<<Types[type_id].X[bead_id].y<<' '<<Types[type_id].X[bead_id].z<<endl;
        }
        continue;
    }else{
        ofs<<s_in<<endl;
    }
  }

  ifs.close();
  ofs.close();
  return;
}