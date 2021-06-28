#include "../include/iodkf.h"
#include "../include/arrays.h"
#include "../include/Legendre.h"
#include "../include/AccelHarmonic.h"
#include <math.h>
#include <stdio.h>

//--------------------------------------------------------------------------
//
// AccelHarmonic.m
//
// Purpose:
//   Computes the acceleration due to the harmonic gravity field of the 
//   central body
//
// Inputs:
//   r           Satellite position vector in the inertial system
//   E           Transformation matrix to body-fixed system
//   n_max       Maximum degree
//   m_max       Maximum order (m_max<=n_max; m_max=0 for zonals, only)
//
// Output:
//   a           Acceleration (a=d^2r/dt^2)
//
// Last modified:   2015/08/12   M. Mahooti
// 
//--------------------------------------------------------------------------

double *AccelHarmonic(double *r, double **E,int nE, int n_max, int m_max){
	extern double **Cnm,**Snm;
	double r_ref,gm,*r_bf,d,latgc,**pnm,**dpnm,dUdlatgc,dUdlon,dUdr,q1,q2,q3,
			b1,b2,b3,*a_bf, r2xy,lon;
	int n,m;
	
	a_bf=vector(3);
	r_ref=6378.1363e3;
	gm=398600.4415e9;
	
	//Body-fixed position
	r_bf=mat_x_vec(E,3,3,r,3);
	//Auxiliary quantities
	d=norma(r_bf,3);
	latgc=asin(r_bf[2]/d);
	lon=atan2(r_bf[1],r_bf[0]);
	
	

	Legendre(n_max,m_max,latgc,&pnm,&dpnm);
	
	dUdr=0;
	dUdlatgc=0;
	dUdlon=0;
	q3=0;
	q2=q3;
	q1=q2;
	
	for(n=0;n<=n_max;++n){
		b1=(-gm/(d*d))*pow((r_ref/d),n)*(n+1);
		b2=(-gm/d)*pow((r_ref/d),n);
		b3=(-gm/d)*pow((r_ref/d),n);
		for(m=0;m<=m_max;++m){
			q1=q1+pnm[n][m]*(Cnm[n][m]*cos(m*lon)+Snm[n][m]*sin(m*lon));
			q2=q2+dpnm[n][m]*(Cnm[n][m]*cos(m*lon)+Snm[n][m]*sin(m*lon));
			q3=q3+m*pnm[n][m]*(Snm[n][m]*cos(m*lon)-Cnm[n][m]*sin(m*lon));
		}
		dUdr=dUdr+q1*b1;
		dUdlatgc=dUdlatgc+q2*b2;
		dUdlon=dUdlon+q3*b3;
		q3=0;
		q2=q3;
		q1=q2;
	}
	
	//Body-fixed acceleration
	
	r2xy=r_bf[0]*r_bf[0]+r_bf[1]*r_bf[1];
	
	a_bf[0]=(1/d*dUdr - r_bf[2]/(d*d*sqrt(r2xy))*dUdlatgc)*r_bf[0]-(1/r2xy*dUdlon)*r_bf[1];
	a_bf[1]=(1/d*dUdr - r_bf[2]/(d*d*sqrt(r2xy))*dUdlatgc)*r_bf[1]+(1/r2xy*dUdlon)*r_bf[0];
	a_bf[2]=1/d*dUdr*r_bf[2]+sqrt(r2xy)/(d*d)*dUdlatgc;
	
	//Inertial acceleration
	
	freeArray(pnm,n_max+2,m_max+2);
	freeArray(dpnm,n_max+2,m_max+2);
	
	return(mat_x_vec(trasp(E,nE),nE,nE,a_bf,3));
	
	
	
}