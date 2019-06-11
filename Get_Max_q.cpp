/*
Author: Yihao Liang
liangyihaosjtu@gmail.com
This code is for Event Chain Monte Carlo for pairwise interacting many body system
*/
#include "public.hpp"
#include "Ewald_Sum_Factor.hpp"
#include <limits>

double4 get_max_q_xface(double xb,double yb1,double yb2,double zb1,double zb2, Ewald_Sum&ES,int axis);
double4 get_max_q_yface(double yb,double xb1,double xb2,double zb1,double zb2, Ewald_Sum&ES,int axis);
double4 get_max_q_zface(double zb,double xb1,double xb2,double yb1,double yb2, Ewald_Sum&ES,int axis);

double get_max_q(double xb1, double xb2, double yb1, double yb2, double zb1, double zb2, Ewald_Sum&ES,int axis){
	double4 M;
	double q_max=0;
	M=get_max_q_xface(xb1, yb1, yb2, zb1, zb2, ES, axis);if(q_max<M.w)q_max=M.w;
	M=get_max_q_xface(xb2, yb1, yb2, zb1, zb2, ES, axis);if(q_max<M.w)q_max=M.w;
	M=get_max_q_yface(yb1, xb1, xb2, zb1, zb2, ES, axis);if(q_max<M.w)q_max=M.w;
	M=get_max_q_yface(yb2, xb1, xb2, zb1, zb2, ES, axis);if(q_max<M.w)q_max=M.w;
	M=get_max_q_zface(zb1, xb1, xb2, yb1, yb2, ES, axis);if(q_max<M.w)q_max=M.w;
	M=get_max_q_zface(zb2, xb1, xb2, yb1, yb2, ES, axis);if(q_max<M.w)q_max=M.w;
	return q_max;
}

double4 get_max_q_xface(double xb,double yb1,double yb2,double zb1,double zb2, Ewald_Sum&ES,int axis){
			//search face x==xb, yb1<y<yb2, zb1<z<zb2
			//const int NB=50;
			const int NB=20;
			double hy,hz;
			
			double4 X1;
			X1.w=1;X1.x=0;X1.y=0;X1.z=0;//fixed

			double4 X2,G;
			X2.w=1;
			X2.x=xb;

			hy=(yb2-yb1)/NB;
			hz=(zb2-zb1)/NB;

			double q,max_q=-1,my0,mz0;

			for(X2.y=yb1;X2.y<=yb2+EPSILON*(yb2-yb1);X2.y+=hy)
				for(X2.z=zb1;X2.z<=zb2+EPSILON*(zb2-zb1);X2.z+=hz) {

					if(X2.y>yb2)X2.y=yb2;
					if(X2.z>zb2)X2.z=zb2;

					q=ES.D_Potential(X1,X2,axis);
					if(q<0)q=0;
					if(q>=max_q){max_q=q;my0=X2.y;mz0=X2.z;}
				}

			double my,mz;
			my=my0;
			mz=mz0;

			for(X2.y=my0-hy;X2.y<=my0+hy;X2.y+=(hy/NB))
				for(X2.z=mz0-hz;X2.z<=mz0+hz;X2.z+=(hz/NB)) {
					if(X2.y<yb1)continue;
					if(X2.z<zb1)continue;

					if(X2.y>yb2)continue;
					if(X2.z>zb2)continue;

					q=ES.D_Potential(X1,X2,axis);
					if(q<0)q=0;
					if(q>=max_q){max_q=q;my=X2.y;mz=X2.z;}
				}

			X2.w=max_q;
			X2.y=my;
			X2.z=mz;
			return X2;
}


double4 get_max_q_yface(double yb,double xb1,double xb2,double zb1,double zb2, Ewald_Sum&ES,int axis){
			//search face xb1<x<xb2, y=yb, zb1<z<zb2
			//const int NB=50;
			const int NB=20;
			double hx,hz;
			
			double4 X1;
			X1.w=1;//fixed
			X1.x=0;//fixed
			X1.y=0;//fixed
			X1.z=0;//fixed

			double4 X2,G;
			X2.w=1;
			X2.y=yb;

			hx=(xb2-xb1)/NB;
			hz=(zb2-zb1)/NB;

			double q,max_q=-1,mx0,mz0;

			for(X2.x=xb1;X2.x<=xb2+EPSILON*(xb2-xb1);X2.x+=hx)
				for(X2.z=zb1;X2.z<=zb2+EPSILON*(zb2-zb1);X2.z+=hz) {

					//if(X2.x>xb2)X2.x=xb2;
					if(X2.x>xb2)X2.x=xb2;
					if(X2.z>zb2)X2.z=zb2;

					q=ES.D_Potential(X1,X2,axis);
					if(q<0)q=0;
					if(q>=max_q){max_q=q;mx0=X2.x;mz0=X2.z;}
				}

			double mx,mz;
			mx=mx0;
			mz=mz0;

			for(X2.x=mx0-hx;X2.x<=mx0+hx;X2.x+=(hx/NB))
				for(X2.z=mz0-hz;X2.z<=mz0+hz;X2.z+=(hz/NB)) {
					if(X2.x<xb1)continue;
					if(X2.z<zb1)continue;

					if(X2.x>xb2)continue;
					if(X2.z>zb2)continue;

					q=ES.D_Potential(X1,X2,axis);
					if(q<0)q=0;
					if(q>=max_q){max_q=q;mx=X2.x;mz=X2.z;}
				}

			X2.w=max_q;
			X2.x=mx;
			X2.z=mz;
			return X2;
}


double4 get_max_q_zface(double zb,double xb1,double xb2,double yb1,double yb2, Ewald_Sum&ES,int axis){
			//search face xb1<x<xb2, yb1<y<yb2, z=zb
			//const int NB=50;
			const int NB=20;
			double hx,hy;
			
			double4 X1;
			X1.w=1;//fixed
			X1.x=0;//fixed
			X1.y=0;//fixed
			X1.z=0;//fixed

			double4 X2,G;
			X2.w=1;
			X2.z=zb;

			hx=(xb2-xb1)/NB;
			hy=(yb2-yb1)/NB;

			double q,max_q=-1,mx0,my0;

			for(X2.x=xb1;X2.x<=xb2+EPSILON*(xb2-xb1);X2.x+=hx)
				for(X2.y=yb1;X2.y<=yb2+EPSILON*(yb2-yb1);X2.y+=hy) {

					if(X2.x>xb2)X2.x=xb2;
					if(X2.y>yb2)X2.y=yb2;

					q=ES.D_Potential(X1,X2,axis);
					if(q<0)q=0;
					if(q>=max_q){max_q=q;mx0=X2.x;my0=X2.y;}
				}

			double mx,my;
			mx=mx0;
			my=my0;

			for(X2.x=mx0-hx;X2.x<=mx0+hx;X2.x+=(hx/NB))
				for(X2.y=my0-hy;X2.y<=my0+hy;X2.y+=(hy/NB)) {
					if(X2.y<yb1)continue;
					if(X2.x<xb1)continue;

					if(X2.y>yb2)continue;
					if(X2.x>xb2)continue;

					q=ES.D_Potential(X1,X2,axis);
					if(q<0)q=0;
					if(q>=max_q){max_q=q;mx=X2.x;my=X2.y;}
				}

			X2.w=max_q;
			X2.x=mx;
			X2.y=my;
			return X2;
}



double get_min_q(double xb1, double xb2, double yb1, double yb2, double zb1, double zb2, Ewald_Sum&ES,int axis){
	double4 M;
	double q_min;
	q_min=numeric_limits<double>::max();

	//search face x==xb, yb1<y<yb2, zb1<z<zb2
	const int NB=10;
	double hx,hy,hz;

	double4 X1;
	X1.w=1;X1.x=0;X1.y=0;X1.z=0;//fixed

	double4 X2,G;
	X2.w=1;

	hx=(xb2-xb1)/NB;
	hy=(yb2-yb1)/NB;
	hz=(zb2-zb1)/NB;

	double q,mx0,my0,mz0;

	for(X2.x=xb1;X2.x<=xb2+EPSILON*(xb2-xb1);X2.x+=hx)
		for(X2.y=yb1;X2.y<=yb2+EPSILON*(yb2-yb1);X2.y+=hy)
			for(X2.z=zb1;X2.z<=zb2+EPSILON*(zb2-zb1);X2.z+=hz) {
				if(X2.x>xb2)X2.x=xb2;
				if(X2.y>yb2)X2.y=yb2;
				if(X2.z>zb2)X2.z=zb2;

				q=ES.D_Potential(X1,X2,axis);
				if(q<0)return 0;
				if(q<q_min){
					q_min=q;
					mx0=X2.x;
					my0=X2.y;
					mz0=X2.z;
				}
			}

			double mx,my,mz;
			mx=mx0;
			my=my0;
			mz=mz0;

	for(X2.x=mx0-hx;X2.x<=mx0+hx;X2.x+=(hx/NB))		
		for(X2.y=my0-hy;X2.y<=my0+hy;X2.y+=(hy/NB))
			for(X2.z=mz0-hz;X2.z<=mz0+hz;X2.z+=(hz/NB)) {
					if(X2.x<xb1)continue;
					if(X2.y<yb1)continue;
					if(X2.z<zb1)continue;

					if(X2.x>xb2)continue;
					if(X2.y>yb2)continue;
					if(X2.z>zb2)continue;

					q=ES.D_Potential(X1,X2,axis);
					if(q<0)return 0;
					if(q<q_min){
						q_min=q;
						mx=X2.x;
						my=X2.y;
						mz=X2.z;
					}
				}

			return q_min;
}
