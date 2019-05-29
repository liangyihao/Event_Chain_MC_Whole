/*
Author: Yihao Liang
liangyihaosjtu@gmail.com
This code is for Event Chain Monte Carlo for pairwise interacting many body system
*/
#include "Ewald_Sum_Factor.hpp"
#include "Random_Number.hpp"
#include <cmath>
#include <iostream>


Ewald_Sum::Ewald_Sum(double Lx,double Ly,double Lz,double accuracy){
    Reset(Lx,Ly,Lz,accuracy);
}

void Ewald_Sum::Reset(double Lx,double Ly,double Lz,double accuracy){
    this->Lx=Lx;
    this->Ly=Ly;
    this->Lz=Lz;
    this->accuracy=accuracy;
    double max_L,min_L;
    max_L=max(max(Lx,Ly),Lz);
    min_L=min(min(Lx,Ly),Lz);
    if(max_L>2*min_L){cout<<"Error, L_max>2*L_min"<<endl;exit(0);}
    alpha=sqrt(M_PI/(min_L*max_L))*1.3;//alpha=sqrt(M_PI/(min_L*max_L))*(1.0~2.5)
    cut_off_n=sqrt(-log(accuracy))/(alpha*min_L)+1;
    cut_off_m=sqrt(-log(accuracy))*(alpha*max_L/M_PI)+1;
    cout<<"Cutoff in Real space: "<<cut_off_n<<endl;
    cout<<"Cutoff in Fourier space: "<<cut_off_m<<endl;

    if(cut_off_m>=20){
        cout<<"Ewald sum init error 1: Too accurate"<<endl;
        exit(0);
    }

    if(cut_off_n>=20){
        cout<<"Ewald sum init error 2: Too accurate"<<endl;
        exit(0);
    }

    double qx,qy,qz,q2;
    for(int mx=0;mx<=cut_off_m;mx++)
        for(int my=0;my<=cut_off_m;my++){
            for(int mz=0;mz<=cut_off_m;mz++){
                qx=2*M_PI/Lx*mx;
                qy=2*M_PI/Ly*my;
                qz=2*M_PI/Lz*mz;
                q2=qx*qx+qy*qy+qz*qz;
                B[mx][my][mz]=4.0*M_PI/(Lx*Ly*Lz)*exp(-q2/(4.0*alpha*alpha))/q2;
                C[mx][my][mz]=pow(2,(mx!=0)+(my!=0)+(mz!=0))*B[mx][my][mz];
                //cout<<C[mx][my][mz]<<endl;
        }
    }
    B[0][0][0]=0;
    C[0][0][0]=0;
}

double Ewald_Sum::Event_Time_Colomb(double4 X1,double4 X2,int axis_index,double Bjerrum_Length, double Max_Event_T){
	//Move X1 along axis_index direction, get reject event time of spring interaction with X2, with parameters Params.
    //If the event time is less than Max_Event_T and less than period in the moving direction, return real event time.
    //If the event time is greater than Max_Event_T, return 2*(period of moving direction)
    //If the event time is greater than one period of moving direction, return 2*(period of moving direction)
	int sign=axis_index/abs(axis_index);
	double4 Y_active,Y_target;
	Y_active.x=X1.x-X2.x;
	Y_active.y=X1.y-X2.y;
	Y_active.z=X1.z-X2.z;
	Y_active.w=X1.w;
	Y_target.x=0;
	Y_target.y=0;
	Y_target.z=0;
	Y_target.w=X2.w;

	while(Y_active.x<-Lx/2)Y_active.x+=Lx;
	while(Y_active.x>+Lx/2)Y_active.x-=Lx;

	while(Y_active.y<-Ly/2)Y_active.y+=Ly;
	while(Y_active.y>+Ly/2)Y_active.y-=Ly;

	while(Y_active.z<-Lz/2)Y_active.z+=Lz;
	while(Y_active.z>+Lz/2)Y_active.z-=Lz;


	double*x_dir_pointer;
	double x_dir0;
	double L;
	double d2_min;
	if(abs(axis_index)==1){
		x_dir_pointer=&(Y_active.x);
		L=Lx;
		d2_min=(Y_active.y)*(Y_active.y)+(Y_active.z)*(Y_active.z);
        if(abs(d2_min)<(1E-30)*L)Y_active.y=(1E-30)*L;
	}
	if(abs(axis_index)==2){
		x_dir_pointer=&(Y_active.y);
		L=Ly;
		d2_min=(Y_active.z)*(Y_active.z)+(Y_active.x)*(Y_active.x);
        if(abs(d2_min)<(1E-30)*L)Y_active.z=(1E-30)*L;
	}
	if(abs(axis_index)==3){
		x_dir_pointer=&(Y_active.z);
		L=Lz;
		d2_min=(Y_active.x)*(Y_active.x)+(Y_active.y)*(Y_active.y);
        if(abs(d2_min)<(1E-30)*L)Y_active.x=(1E-30)*L;
	}

	(*x_dir_pointer)=(*x_dir_pointer)*sign;
	x_dir0=(*x_dir_pointer);

    double S1,S2;
	if(x_dir0<0){
		S1=-x_dir0;
		S2=-x_dir0+L/2.0;
	}else{
		S1=-x_dir0+L/2.0;
		S2=-x_dir0+L;
	}

	/*
	From s=0 to s=L, we have 3 parts:
	[0,S1],[S1,S2],[S2,L]
	*/
    double q;
    q=-log(1-Uniform_Random());
    double s_left,s_right;
    double Q_left,Q_right;
    double Prod_Charge=X1.w*X2.w;
    
    //double temp1,temp2;
    double Potential_Left,Potential_Right;

    if((Prod_Charge>0)&&(x_dir0<0)){
	    	//search [0,S1],[S2,L]
	    	s_left=0;
	    	Q_left=0;
            
            s_right=S1;
            (*x_dir_pointer)=x_dir0+s_left;
            Potential_Left=Potential(Y_active, Y_target);
            (*x_dir_pointer)=x_dir0+s_right;
            Potential_Right=Potential(Y_active, Y_target);
            Q_right=Q_left+Bjerrum_Length*(Potential_Right-Potential_Left);

            if(q>Q_right){
            	s_left=S2;
                Q_left=Q_right;

            	s_right=L;
                (*x_dir_pointer)=x_dir0+s_left;
                Potential_Left=Potential(Y_active, Y_target);
                (*x_dir_pointer)=x_dir0+s_right;
                Potential_Right=Potential(Y_active, Y_target);
                Q_right=Q_left+Bjerrum_Length*(Potential_Right-Potential_Left);
            }
    }else if((Prod_Charge>0)&&(x_dir0>=0)){
	    	//only search [S1,S2]
	    	s_left=S1;
	    	Q_left=0;

            s_right=S2;
            (*x_dir_pointer)=x_dir0+s_left;
            Potential_Left=Potential(Y_active, Y_target);
            (*x_dir_pointer)=x_dir0+s_right;
            Potential_Right=Potential(Y_active, Y_target);
            Q_right=Q_left+Bjerrum_Length*(Potential_Right-Potential_Left);
    }else if((Prod_Charge<0)&&(x_dir0<=0)){
            //only search [S1,S2]
	    	s_left=S1;
	    	Q_left=0;

            s_right=S2;
            (*x_dir_pointer)=x_dir0+s_left;
            Potential_Left=Potential(Y_active, Y_target);
            (*x_dir_pointer)=x_dir0+s_right;
            Potential_Right=Potential(Y_active, Y_target);
            Q_right=Q_left+Bjerrum_Length*(Potential_Right-Potential_Left);
    }else if((Prod_Charge<0)&&(x_dir0>0)){
            //search [0,S1],[S2,L]
	    	s_left=0;
	    	Q_left=0;
            
            s_right=S1;
            (*x_dir_pointer)=x_dir0+s_left;
            Potential_Left=Potential(Y_active, Y_target);
            (*x_dir_pointer)=x_dir0+s_right;
            Potential_Right=Potential(Y_active, Y_target);
            Q_right=Q_left+Bjerrum_Length*(Potential_Right-Potential_Left);

            if(q>Q_right){
            	s_left=S2;
                Q_left=Q_right;

            	s_right=L;
                (*x_dir_pointer)=x_dir0+s_left;
                Potential_Left=Potential(Y_active, Y_target);
                (*x_dir_pointer)=x_dir0+s_right;
                Potential_Right=Potential(Y_active, Y_target);
                Q_right=Q_left+Bjerrum_Length*(Potential_Right-Potential_Left);
            }
    }
    if(q>Q_right)return 2*L;
        
    //Find solution
    double s_mid,Q_mid;
    double Potential_Mid;

    double s_left0,Q_left0,Potential_Left0;
    s_left0=s_left;
    Q_left0=Q_left;
    Potential_Left0=Potential_Left;

    while((abs(Q_right-Q_left)>(accuracy*Bjerrum_Length/L))&&(abs(s_right-s_left)>EPSILON*L)) {
        if(s_left>Max_Event_T)return 2*L;//This trick is for saving time
        s_mid=(s_left+s_right)/2;
        (*x_dir_pointer)=x_dir0+s_mid;
        Potential_Mid=Potential(Y_active, Y_target);
        Q_mid=Q_left0+Bjerrum_Length*(Potential_Mid-Potential_Left0);

        if(Q_mid>q){
            Q_right=Q_mid;
            s_right=s_mid;
            Potential_Right=Potential_Mid;
        }else {
            Q_left=Q_mid;
            s_left=s_mid;
            Potential_Left=Potential_Mid;
        }
    }

    return s_mid;

}

double Ewald_Sum::Potential(double4 X,double4 Y){
    double res_R=0,res_k=0;
    double dx0,dy0,dz0;
    double c1,c2;
    c1=X.w;
    c2=Y.w;
    dx0=X.x-Y.x;
    dy0=X.y-Y.y;
    dz0=X.z-Y.z;
    double dx,dy,dz,dr;
    for(int nx=-cut_off_n;nx<=cut_off_n;nx++){
        dx=dx0+nx*Lx;
        for(int ny=-cut_off_n;ny<=cut_off_n;ny++){
            dy=dy0+ny*Ly;
            for(int nz=-cut_off_n;nz<=cut_off_n;nz++){
                dz=dz0+nz*Lz;
                dr=sqrt(dx*dx+dy*dy+dz*dz);
                res_R+=erfc(alpha*dr)/dr;
            }
        }
    }
    double Qx,Qy,Qz;
    double qx0,qy0,qz0;
    qx0=2.0*M_PI*dx0/Lx;
    qy0=2.0*M_PI*dy0/Ly;
    qz0=2.0*M_PI*dz0/Lz;
    for(int mx=0;mx<=cut_off_m;mx++){
        Qx=qx0*mx;
        for(int my=0;my<=cut_off_m;my++){
            Qy=qy0*my;
            for(int mz=0;mz<=cut_off_m;mz++){
                Qz=qz0*mz;
                res_k+=C[mx][my][mz]*cos(Qx)*cos(Qy)*cos(Qz);
            }
        }
    }
    return c1*c2*(res_R+res_k);
}


double Ewald_Sum::Self_E(double4 X){
    double res=0;
    double dx,dy,dz,dr;
    for(int nx=-cut_off_n;nx<=cut_off_n;nx++){
        dx=nx*Lx;
        for(int ny=-cut_off_n;ny<=cut_off_n;ny++){
            dy=ny*Ly;
            for(int nz=-cut_off_n;nz<=cut_off_n;nz++){
                dz=nz*Lz;
                dr=sqrt(dx*dx+dy*dy+dz*dz);
                if(dr<1E-20)continue;
                res+=erfc(alpha*dr)/dr;
            }
        }
    }


    for(int mx=0;mx<=cut_off_m;mx++)
        for(int my=0;my<=cut_off_m;my++){
            for(int mz=0;mz<=cut_off_m;mz++){
                res+=C[mx][my][mz];
            }
        }
    res=res/2-alpha/sqrt(M_PI);
    return res*X.w*X.w;
}

double Ewald_Sum::D_Potential(double4 X,double4 Y,int axis){
    double res_R=0,res_k=0;
    double dx0,dy0,dz0;
    double c1,c2;
    c1=X.w;
    c2=Y.w;
    dx0=X.x-Y.x;
    dy0=X.y-Y.y;
    dz0=X.z-Y.z;
    double dx,dy,dz,dr;
    double qx,qy,qz;
    
    double*dx_axis;
    double*q_axis;

    double(*fx)(double);
    double(*fy)(double);
    double(*fz)(double);

    fx=cos;
    fy=cos;
    fz=cos;

    if(abs(axis)==1){
        if(axis==-1)dx0=-dx0;
        dx_axis=&dx;
        q_axis=&qx;    
        fx=sin;
    }
    if(abs(axis)==2){
        if(axis==-2)dy0=-dy0;
        dx_axis=&dy;
        q_axis=&qy;
        fy=sin;
    }
    if(abs(axis)==3){
        if(axis==-3)dz0=-dz0;
        dx_axis=&dz;
        q_axis=&qz;
        fz=sin;
    }


    for(int nx=-cut_off_n;nx<=cut_off_n;nx++){
        dx=dx0+nx*Lx;
        for(int ny=-cut_off_n;ny<=cut_off_n;ny++){
            dy=dy0+ny*Ly;
            for(int nz=-cut_off_n;nz<cut_off_n;nz++){
                dz=dz0+nz*Lz;
                dr=sqrt(dx*dx+dy*dy+dz*dz);
                res_R+=(erfc(alpha*dr)/dr+2*alpha*exp(-alpha*alpha*dr*dr)/sqrt(M_PI))*(*dx_axis)/(dr*dr);
            }
        }
    }
    double qx0,qy0,qz0;
    qx0=2.0*M_PI/Lx;
    qy0=2.0*M_PI/Ly;
    qz0=2.0*M_PI/Lz;
    
    for(int mx=0;mx<=cut_off_m;mx++){
        qx=qx0*mx;
        for(int my=0;my<=cut_off_m;my++){
            qy=qy0*my;
            for(int mz=0;mz<=cut_off_m;mz++){
                qz=qz0*mz;
                res_k+=(*q_axis)*C[mx][my][mz]*fx(qx*dx0)*fy(qy*dy0)*fz(qz*dz0);
            }
        }
    }
    return -(res_R+res_k)*c1*c2;
}


double Ewald_Sum::D_Potential_Check(double4 X,double4 Y,int axis){
        double4 XX;
        double ep=1E-9;
        double epsilon;
        XX=X;
        int sign;
        sign=axis/abs(axis);
        if(abs(axis)==1)epsilon=sign*(ep)*Lx;
        if(abs(axis)==2)epsilon=sign*(ep)*Ly;
        if(abs(axis)==3)epsilon=sign*(ep)*Lz;

        if(abs(axis)==1)XX.x+=epsilon;
        if(abs(axis)==2)XX.y+=epsilon;
        if(abs(axis)==3)XX.z+=epsilon;

        double P1,P2;
        P1=Potential(X,Y);
        P2=Potential(XX,Y);
        return (P2-P1)/abs(epsilon);
}
