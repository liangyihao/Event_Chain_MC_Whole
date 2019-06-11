/*
Author: Yihao Liang
liangyihaosjtu@gmail.com
This code is for Event Chain Monte Carlo for pairwise interacting many body system
*/
//Cell-Veto List
//deal with long-range interactions
//range:[0,Lx],[0,Ly],[0,Lz]
#include <cmath>
#include "CellVetoList.hpp"
#include "Get_Max_q.hpp"
#include "public.hpp"
using namespace std;

void CellVetoList::Get_Colomb_Event_Exceptional_Particles(int2 id_active_particle, int axis, double Bjerrum_Length, double&time, int2&id_next_active_bead){
//Get first coloumb event from exceptional particles
//Note: time stores the current min event time
//Note: id_next_active_bead stores the current next active bead
    	int2 ids;
	    double4 X1,X2;
    	double t;
	    X1=(*Types_pointer)[id_active_particle.x].X[id_active_particle.y];
        //go over exceptional particle, they are always computed directly in this call
	    for(int i=0;i<Exception_Particle_List.size();i++){
            ids=Exception_Particle_List[i];
            if((ids.x==id_active_particle.x)&&(ids.y==id_active_particle.y))continue;
	    	X2=(*Types_pointer)[ids.x].X[ids.y];
		    t=ES->Event_Time_Colomb(X1,X2,axis,Bjerrum_Length,time);
    		if(t<time){
	    		time=t;
		    	id_next_active_bead=ids;
		    }
    	}
}


void CellVetoList::Get_Colomb_Event_Cell_Veto(int2 id_active_particle, int axis, double Bjerrum_Length, double&time, int2&id_next_active_bead){
//Get first coloumb event from veto-grid
//Note: time stores the current min event time
//Note: id_next_active_bead stores the current next active bead
    	double4 X1,X1_t;
	    X1=(*Types_pointer)[id_active_particle.x].X[id_active_particle.y];
        X1_t=X1;
        int3 IWC;
        IWC=In_Which_Veto_Cell(X1);

        double Q_axis_tot,dL_axis,L_axis;
        int *IWC_axis;
        
        double*X1_axis;
        double*X1_t_axis;
        Frequency_Generator*FG_axis;
        int NC_axis;

        if(abs(axis)==1){
            Q_axis_tot=Qx_tot;
            dL_axis=dLx;
            L_axis=Lx;
            IWC_axis=(&(IWC.x));
            X1_axis=(&(X1.x));
            X1_t_axis=(&(X1_t.x));
            FG_axis=FG_x;
            NC_axis=NC_x;
        }

        if(abs(axis)==2){
            Q_axis_tot=Qy_tot;
            dL_axis=dLy;
            L_axis=Ly;
            IWC_axis=(&(IWC.y));
            X1_axis=(&(X1.y));
            X1_t_axis=(&(X1_t.y));
            FG_axis=FG_y;
            NC_axis=NC_y;
        }

        if(abs(axis)==3){
            Q_axis_tot=Qz_tot;
            dL_axis=dLz;
            L_axis=Lz;
            IWC_axis=(&(IWC.z));
            X1_axis=(&(X1.z));
            X1_t_axis=(&(X1_t.z));
            FG_axis=FG_z;
            NC_axis=NC_z;
        }
        
        if(((*IWC_axis)==0)&&(abs((*X1_t_axis)-L_axis)<EPSILON*L_axis)){//!!!need to make tighter condition
            (*X1_t_axis)-=L_axis;
        }
        int sign=axis/abs(axis);
        double clock=0;
        double next_bound_clock,next_bound_coordinate;
        if(sign==1){
            next_bound_coordinate=((*IWC_axis)+1)*dL_axis;
            next_bound_clock=((*IWC_axis)+1)*dL_axis-(*X1_t_axis);
        }else if(sign==-1){
            next_bound_coordinate=(*IWC_axis)*dL_axis;
            next_bound_clock=(*X1_t_axis)-(*IWC_axis)*dL_axis;
        }

    	do{//loop for each possible cell
            double Time=time;
            int2 Id_Next_Active_Bead=id_next_active_bead;
            int2 ids;
            double4 X2;
            int3 IWC2,IWCR;
            double temp;

	    	//step 1: go over all the charged particles in the nearest neighboring cells
            #ifdef DEBUG
                //cout<<"clock:"<<clock<<endl;
                //cout<<"axis:"<<axis<<endl;
                //cout<<"temp coordinate of active particle:"<<X1_t.x<<' '<<X1_t.y<<' '<<X1_t.z<<endl;
                //cout<<"In Cell: "<<IWC.x<<' '<<IWC.y<<' '<<IWC.z<<endl;
            #endif
			for(int i=IWC.x-1;i<=IWC.x+1;i++)
				for(int j=IWC.y-1;j<=IWC.y+1;j++)
					for(int k=IWC.z-1;k<=IWC.z+1;k++){
						IWC2.x=i;IWC2.y=j;IWC2.z=k;
                        
                        if(IWC2.x>=NC_x)IWC2.x-=NC_x;
                        if(IWC2.x<0)IWC2.x+=NC_x;

                        if(IWC2.y>=NC_y)IWC2.y-=NC_y;
                        if(IWC2.y<0)IWC2.y+=NC_y;

                        if(IWC2.z>=NC_z)IWC2.z-=NC_z;
                        if(IWC2.z<0)IWC2.z+=NC_z;

                        for(int l=0;l<min((int)(Veto_Cells[IWC2.x][IWC2.y][IWC2.z].particle_list.size()),Num_Particle_Per_Cell);l++){
							ids=Veto_Cells[IWC2.x][IWC2.y][IWC2.z].particle_list[l];
                            if((ids.x==id_active_particle.x)&&(ids.y==id_active_particle.y))continue;
							X2=(*Types_pointer)[ids.x].X[ids.y];
							temp=ES->Event_Time_Colomb(X1_t,X2,axis,Bjerrum_Length,next_bound_clock-clock);//need check
							if(temp+clock<min(Time,L_axis)){//need use another one
								Time=temp+clock;//need use another one
								Id_Next_Active_Bead=ids;//need use another one
                                //cout<<"From nearest neighbors, event time updated to:"<<Time<<endl;
							}
						}
                    }
		    //step 2: generate events by cell-veto
            do{
                temp=Exponential_Random(abs(Q_axis_tot*X1.w*valence*Bjerrum_Length*Num_Particle_Per_Cell));//pre-event
                if(clock+temp>Time){//check if group event later than other known event
                    if(Time<next_bound_clock){
                        time=Time;
                        id_next_active_bead=Id_Next_Active_Bead;
                        //cout<<"reject at event time"<<time<<endl;
                        return;
                    }else{
                        (*X1_t_axis)=next_bound_coordinate;
                        clock=next_bound_clock;
                        break;
                    }
                }
                if(clock+temp>L_axis){//check if group event later than one period
                    time=Time;
                    id_next_active_bead=Id_Next_Active_Bead;
                    //cout<<"reject at event time"<<time<<endl;
                    return;
                }
                if(clock+temp>next_bound_clock){//check if group event later than the cell boundary, if yes, break
                    (*X1_t_axis)=next_bound_coordinate;
                    clock=next_bound_clock;
                    break;
                }

                clock+=temp;
                (*X1_t_axis)+=temp*sign;

                //assign it to a cell
                IWCR=Int_To_Int3(FG_axis->Gen());
                double q_max_cell;
                if(abs(axis)==1)q_max_cell=qx_max[IWCR.x][IWCR.y][IWCR.z];
                if(abs(axis)==2)q_max_cell=qy_max[IWCR.x][IWCR.y][IWCR.z];
                if(abs(axis)==3)q_max_cell=qz_max[IWCR.x][IWCR.y][IWCR.z];
                q_max_cell*=abs(X1.w*valence);

                #ifdef DEBUG
                double q_min_cell;
                    if(abs(axis)==1)q_min_cell=qx_min[IWCR.x][IWCR.y][IWCR.z];
                    if(abs(axis)==2)q_min_cell=qy_min[IWCR.x][IWCR.y][IWCR.z];
                    if(abs(axis)==3)q_min_cell=qz_min[IWCR.x][IWCR.y][IWCR.z];
                    q_min_cell*=abs(X1.w*valence);
                #endif

                if(axis*X1.w*valence<0){
                    if(abs(axis)==1)IWCR.x*=-1;
                    if(abs(axis)==2)IWCR.y*=-1;
                    if(abs(axis)==3)IWCR.z*=-1;
                }

                IWC2.x=IWC.x+IWCR.x;
                IWC2.y=IWC.y+IWCR.y;
                IWC2.z=IWC.z+IWCR.z;
                if(IWC2.x>=NC_x)IWC2.x-=NC_x;
                if(IWC2.y>=NC_y)IWC2.y-=NC_y;
                if(IWC2.z>=NC_z)IWC2.z-=NC_z;
                if(IWC2.x< 0)IWC2.x+=NC_x;
                if(IWC2.y< 0)IWC2.y+=NC_y;
                if(IWC2.z< 0)IWC2.z+=NC_z;

                //assign it to a particle in that cell
                if(Veto_Cells[IWC2.x][IWC2.y][IWC2.z].particle_list.size()==0)continue;
                double p;
                p=Uniform_Random()*Num_Particle_Per_Cell;
                if(p>=Veto_Cells[IWC2.x][IWC2.y][IWC2.z].particle_list.size())continue;
                ids=Veto_Cells[IWC2.x][IWC2.y][IWC2.z].particle_list[(int)p];
                if((ids.x==id_active_particle.x)&&(ids.y==id_active_particle.y))continue;
                
                //if assignment succeed,
                X2=X2=(*Types_pointer)[ids.x].X[ids.y];
                double q;
                q=ES->D_Potential_Ramp(X1_t,X2,axis);

                #ifdef DEBUG
                    double q_safe;
                    q_safe=max(ES->D_Potential(X1_t,X2,axis),0.0);
                    if(abs(q-q_safe)>EPSILON){
                        cout<<"Warning q incorrect q="<<q<<" q(std)="<<q_safe<<endl;
                        cout<<"active particle position:"<<X1_t.x<<','<<X1_t.y<<','<<X1_t.z<<" charge"<<X1_t.w<<endl;
                        cout<<"target particle position:"<<X2.x<<','<<X2.y<<','<<X2.z<<" charge"<<X2.w<<endl;
                        cout<<"axis: ";
                        if(axis==+1)cout<<"+x"<<endl;
                        if(axis==-1)cout<<"-x"<<endl;
                        if(axis==+2)cout<<"+y"<<endl;
                        if(axis==-2)cout<<"-y"<<endl;
                        if(axis==+3)cout<<"+z"<<endl;
                        if(axis==-3)cout<<"-z"<<endl;
                        cout<<endl;
                    }
                    if((q-q_max_cell)*(q-q_min_cell)>0){
                        cout<<"Warning, q out of bound. q="<<q<<" not in ["<<q_min_cell<<","<<q_max_cell<<"]"<<endl;
                        int Ir,Jr,Kr;
                        Ir=IWC2.x-IWC.x;Jr=IWC2.y-IWC.y;Kr=IWC2.z-IWC.z;
                        if(axis==-1)Ir*=-1;
                        if(axis==-2)Jr*=-1;
                        if(axis==-3)Kr*=-1;
                        while(Ir<0)Ir+=NC_x;while(Ir>NC_x)Ir-=NC_x;
                        while(Jr<0)Jr+=NC_y;while(Jr>NC_y)Jr-=NC_y;
                        while(Kr<0)Kr+=NC_z;while(Kr>NC_z)Kr-=NC_z;

                        cout<<"The cell index is"<<Ir<<' '<<Jr<<' '<<Kr<<endl<<endl;
                    }
                #endif
                if(Uniform_Random()<q/q_max_cell){//check if it is a real event
                    //real event
                    time=clock;
                    id_next_active_bead=ids;
                    #ifdef DEBUG
                        //cout<<"reject at event time(from cell veto):"<<time<<endl;
                    #endif
                    return;
                }
            }while(1);

        //step 3: update
            (*IWC_axis)+=sign;

            if((*IWC_axis)>=NC_axis){
                (*X1_t_axis)-=L_axis;
                (*IWC_axis)-=NC_axis;
            }
            if((*IWC_axis)<0){
                (*X1_t_axis)+=L_axis;
                (*IWC_axis)+=NC_axis;
            }

            if(sign==1){
                next_bound_coordinate=((*IWC_axis)+1)*dL_axis;
                next_bound_clock+=dL_axis;
                //next_bound_clock=((*IWC_axis)+1)*dL_axis-(*X1_t_axis)+clock;
            }else if(sign==-1){
                next_bound_coordinate=(*IWC_axis)*dL_axis;
                next_bound_clock+=dL_axis;
                //next_bound_clock=(*X1_t_axis)-(*IWC_axis)*dL_axis+clock;
            }

	    }while(1);
}

void CellVetoList::Get_Colomb_Event(int2 id_active_particle, int axis, double Bjerrum_Length, double&time, int2&id_next_active_bead){
//Get first coloumb event
//Note: time stores the current min event time
//Note: id_next_active_bead stores the current next active bead
    Get_Colomb_Event_Exceptional_Particles(id_active_particle, axis, Bjerrum_Length, time, id_next_active_bead);
    Get_Colomb_Event_Cell_Veto(id_active_particle, axis, Bjerrum_Length, time, id_next_active_bead);
}
