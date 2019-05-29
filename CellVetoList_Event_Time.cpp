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
            if((ids.x==id_active_particle.x)&&(ids.y==id_active_particle.y))continue;
    		ids=Exception_Particle_List[i];
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


        int sign=axis/abs(axis);
        double clock=0;
        double next_bound_clock;
        if(sign==1){
            next_bound_clock=((*IWC_axis)+1)*dL_axis-(*X1_t_axis);
        }else if(sign==-1){
            next_bound_clock=(*X1_t_axis)-(*IWC_axis)*dL_axis;
        }

    	do{
            double Time=time;
            int2 Id_Next_Active_Bead=id_next_active_bead;
            int2 ids;
            double4 X2;
            int3 IWC2;
            double temp;

	    	//step 1: go over all the charged particles in the nearest neighboring cells
			for(int i=IWC.x-1;i<=IWC.x+1;i++)
				for(int j=IWC.y-1;j<=IWC.y+1;j++)
					for(int k=IWC.z-1;k<=IWC.z+1;k++)
						for(int l=0;l<Veto_Cells[IWC.x][IWC.y][IWC.z].particle_list.size();l++){
							ids=Veto_Cells[IWC.x][IWC.y][IWC.z].particle_list[l];
                            if((ids.x==id_active_particle.x)&&(ids.y==id_active_particle.y))continue;
							X2=(*Types_pointer)[ids.x].X[ids.y];
							temp=ES->Event_Time_Colomb(X1_t,X2,axis,Bjerrum_Length,next_bound_clock-clock);//need check
							if(temp+clock<Time){//need use another one
								Time=temp+clock;//need use another one
								Id_Next_Active_Bead=ids;//need use another one
							}
						}
		    //step 2: generate events by cell-veto
            do{
                temp=Exponential_Random(abs(Q_axis_tot*X1.w*valence*Bjerrum_Length));//pre-event
                if(clock+temp>next_bound_clock)break;//check if this event later than the cell boundary, if yes, break
                if(clock+temp>Time){//check if this event later than other known event
                    time=Time;
                    id_active_particle=Id_Next_Active_Bead;
                    return;
                }
                if(clock+temp>L_axis){//check if this event later than one period
                    time=Time;
                    id_active_particle=Id_Next_Active_Bead;
                    return;
                }
                clock+=temp;
                (*X1_t_axis)=(*X1_axis)+clock*sign;

                //assign it to a cell
                IWC2=Int_To_Int3(FG_axis->Gen());
                double q_max_cell;
                if(abs(axis)==1)q_max_cell=qx_max[IWC2.x][IWC2.y][IWC2.z];
                if(abs(axis)==2)q_max_cell=qy_max[IWC2.x][IWC2.y][IWC2.z];
                if(abs(axis)==3)q_max_cell=qz_max[IWC2.x][IWC2.y][IWC2.z];
                q_max_cell*=abs(X1.w*valence);

                if(axis*X1.w*valence<0){
                    if(abs(axis)==1)IWC2.x*=-1;
                    if(abs(axis)==2)IWC2.y*=-1;
                    if(abs(axis)==3)IWC2.z*=-1;
                }

                IWC2.x=IWC.x+IWC2.x;
                IWC2.y=IWC.y+IWC2.y;
                IWC2.z=IWC.z+IWC2.z;
                if(IWC2.x>=NC_x)IWC2.x-=NC_x;
                if(IWC2.y>=NC_y)IWC2.y-=NC_y;
                if(IWC2.z>=NC_z)IWC2.z-=NC_z;
                if(IWC2.x< 0)IWC2.x+=NC_x;
                if(IWC2.y< 0)IWC2.y+=NC_y;
                if(IWC2.z< 0)IWC2.z+=NC_z;

                //assign it to a particle in that cell
                double p;
                p=Uniform_Random()*Num_Particle_Per_Cell;
                if(p>=Veto_Cells[IWC2.x][IWC2.y][IWC2.z].particle_list.size())continue;
                ids=Veto_Cells[IWC2.x][IWC2.y][IWC2.z].particle_list[(int)p];
                if((ids.x==id_active_particle.x)&&(ids.y==id_active_particle.y))continue;
                
                //if assignment succeed,
                X2=X2=(*Types_pointer)[ids.x].X[ids.y];
                double q;
                q=max(ES->D_Potential(X1_t,X2,axis),0.0);
                if(Uniform_Random()<q/q_max_cell){//check if it is a real event
                    //real event
                    time=clock;
                    id_next_active_bead=ids;
                    return;
                }
            }while(1);

        //step 3: update
            (*IWC_axis)+=sign;
            if(sign==1){
                next_bound_clock=((*IWC_axis)+1)*dL_axis-(*X1_t_axis);
            }else if(sign==-1){
                next_bound_clock=(*X1_t_axis)-(*IWC_axis)*dL_axis;
            }

            if((*IWC_axis)>=NC_axis){
                (*X1_t_axis)-=L_axis;
                (*IWC_axis)-=NC_axis;
            }
            if((*IWC_axis)<0){
                (*X1_t_axis)+=L_axis;
                (*IWC_axis)+=NC_axis;
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
