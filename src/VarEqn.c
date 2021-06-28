#include <stdio.h> 
#include "../include/globales.h" 
#include "../include/arrays.h" 
#include "../include/iodkf.h" 
#include "../include/IERS.h" 
#include "../include/timediff.h" 
#include "../include/PrecMatrix.h" 
#include "../include/NutMatrix.h"
#include "../include/Mjday_TDB.h" 
#include "../include/JPL_Eph_DE430.h" 
#include "../include/PoleMatrix.h"
#include "../include/GHAMatrix.h"  
#include "../include/Accel.h" 
#include "../include/AccelHarmonic.h"
#include "../include/AccelPointMass.h" 
#include "../include/G_AccelHarmonic.h"
#include "../include/VarEqn.h"

//------------------------------------------------------------------------------
//
// VarEqn.m
//
// Purpose:
//   Computes the variational equations, i.e. the derivative of the state vector
//   and the state transition matrix
//
// Input:
//   x           Time since epoch in [s]
//   yPhi        (6+36)-dim vector comprising the state vector (y) and the
//               state transition matrix (Phi) in column wise storage order
//
// Output:
//   yPhip       Derivative of yPhi
// 
// Last modified:   2015/08/12   M. Mahooti
//
//------------------------------------------------------------------------------

void VarEqn(double x, double *yPhi, double **yPhip){
	
	extern Param AuxParam;
	double x_pole,y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC,UT1_TAI,
		UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC, Mjd_UT1, Mjd_TT, **P, **N, **T, **E, *r, *v, 
		**Phi, *a, **G, **dfdy, **Phip;
	int i,j;
	double *aux=vector(3);
	r=vector(3);
	v=vector(3);
	Phi=array(6,6);
	
	IERS(AuxParam.Mjd_UTC,'l', &x_pole,&y_pole, &UT1_UTC, &LOD, &dpsi, &deps, &dx_pole, &dy_pole, &TAI_UTC); 

	timediff(UT1_UTC, TAI_UTC, &UT1_TAI, &UTC_GPS, &UT1_GPS, &TT_UTC, &GPS_UTC); 
	
	Mjd_UT1 = AuxParam.Mjd_TT + (UT1_UTC-TT_UTC)/86400;

	// Transformation matrix
	P = PrecMatrix(MJD_J2000, AuxParam.Mjd_TT + x/86400); 
	N = NutMatrix(AuxParam.Mjd_TT + x/86400); 
	T=prod(N,3,3,P,3,3);
	E = prod(prod(PoleMatrix(x_pole,y_pole),3,3,GHAMatrix(Mjd_UT1),3,3),3,3,T,3,3);

	// State vector components
	for(i=0;i<3;i++){
		r[i] = yPhi[i];
		v[i] = yPhi[i+3];
	}

	// State transition matrix
	for(j=0;i<6;++j){
		for(i=0;i<6;++i){
			Phi[i][j] = yPhi[6*(j+1)+1+i];
		}
	}

	// Acceleration and gradient
	a = AccelHarmonic ( r, E, 3, 20, 20 );
	G = G_AccelHarmonic ( r, E, 3, 20, 20 );

	// Time derivative of state transition matrix
	dfdy = array(6,6);
	
	for(i=0;i<3;++i){
		for(j=0;j<3;++j){
			
			dfdy[i][j]=0.0;
			dfdy[i+3][j]=G[i][j];
			if(i==j)
				dfdy[i][j+3]=1;
			else
				dfdy[i][j+3]=0;
			dfdy[i+3][j+3]=0.0;
		}
		
	}

	Phip = prod(dfdy,6,6,Phi,6,6);

	// Derivative of combined state vector and state transition matrix
	for(i=0;i<3;++i){
		aux[i]=v[i];
		aux[i+3]=a[i];
	}
	
	for(i=0;i<6;++i){
		for(j=0;j<6;++j){
			aux[5*j+i]=Phip[i][j];
		}
	}
	
	*yPhip=aux;

}