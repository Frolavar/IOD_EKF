#include "../include/G_AccelHarmonic.h"
#include "../include/Legendre.h"
#include "../include/AccelHarmonic.h"
#include "../include/arrays.h"
#include <stdio.h>
//--------------------------------------------------------------------------
//
// G_AccelHarmonic.m
//
// Purpose:
//   Computes the gradient of the Earth's harmonic gravity field 
//
// Inputs:
//   r           Satellite position vector in the true-of-date system
//   U           Transformation matrix to body-fixed system
//   n           Gravity model degree
//   m 			Gravity model order
//
// Output:
//   G    		Gradient (G=da/dr) in the true-of-date system
//
// Last modified:   2015/08/12   M. Mahooti
//
//--------------------------------------------------------------------------
double **G_AccelHarmonic(double *r, double **U, int nU, int n_max,int m_max){
	int i,j;
	double **G,*dr,*da,*da1,*da2,d; 
	
	d=1.0;
	G=array(3,3);
	dr=vector(3);
	
	// Gradient
	for(i=0;i<3;++i){
		    // Set offset in i-th component of the position vector
		for(j=0;j<3;++j)
			dr[j]=0.0;
		dr[i]=d/2;
		
		// Acceleration difference
		da1=AccelHarmonic(sumV(r,3,dr,3),U,nU, n_max, m_max);
		da2=AccelHarmonic(sumV(r,3,esc_x_vec(-1,dr,3),3),U,nU, n_max, m_max );
		da=sumV(da1,3,esc_x_vec(-1,da2,3),3);
		
		// Derivative with respect to i-th axis
		for(j=0;j<3;++j)
			G[j][i]=da[j]/d;
		freeVector(da,3);
		freeVector(da1,3);
		freeVector(da2,3);
	}
	freeVector(dr,3);

	return G;
}
