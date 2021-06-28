#include "../include/Cheb3D.h"
#include "../include/arrays.h"
#include <stdio.h>
#include <stdlib.h>

//--------------------------------------------------------------------------
//
// Chebyshev approximation of 3-dimensional vectors
//
// Inputs:
//     N       Number of coefficients
//     Ta      Begin interval
//     Tb      End interval
//     Cx      Coefficients of Chebyshev polyomial (x-coordinate)
//     Cy      Coefficients of Chebyshev polyomial (y-coordinate)
//     Cz      Coefficients of Chebyshev polyomial (z-coordinate)
//
// Last modified:   2018/01/27   M. Mahooti
// 
//--------------------------------------------------------------------------

double *Cheb3D(double t,int N,double Ta,double Tb,double *Cx,double *Cy,double *Cz){
	double *f1,*f2,*old_f1,*aux,tau;
	int i,j;
	
	// Check validity
	if( (t<Ta) || (Tb<t) ){
		printf("ERROR: Time out of range in Cheb3D::Value\n");
		exit(EXIT_FAILURE);
	}

	// Clenshaw algorithm
	tau = (2*t-Ta-Tb)/(Tb-Ta);  

	f1 = vector(3);
	f2 = vector(3);
	old_f1=vector(3);
	aux=vector(3);
	
	for(i=N-1;i>0;--i){
		for(j=0;j<3;++j){
			old_f1[j]=f1[j];
		}
		aux[0]=Cx[i];
		aux[1]=Cy[i];
		aux[2]=Cz[i];
		f1=sumV(sumV(esc_x_vec(2*tau,f1,3),3,esc_x_vec(-1,f2,3),3),3,aux,3);
	}

	aux[0]=Cx[0];
	aux[1]=Cy[0];
	aux[2]=Cz[0];

	return (sumV(sumV(esc_x_vec(tau,f1,3),3,esc_x_vec(-1,f2,3),3),3,aux,3));
}