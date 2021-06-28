#include <math.h>
#include <stdio.h>
#include "../include/eccanom.h"

//--------------------------------------------------------------------------
//
// Purpose:
//   Computes the eccentric anomaly for elliptic orbits
//
// Inputs:
//   M         Mean anomaly in [rad]
//   e         Eccentricity of the orbit [0,1]
// 
// Output:
//             Eccentric anomaly in [rad]
//
// Last modified:   2015/08/12   M. Mahooti
// 
//--------------------------------------------------------------------------

double EccAnom(double M, double e){
	double E, f, eps=2.22044604925031*1e-16,M_pi=3.14159265358979;
	int i, maxit;
	
	maxit = 15;
	i = 1;


	M = fmod(M, 2.0*M_pi);

	if (e<0.8)
		E = M; 
	else
		E = M_pi;

	f = E - e*sin(E) - M;
	E = E - f / ( 1.0 - e*cos(E) );

	while (fabs(f) > 1e2*eps){
		f = E - e*sin(E) - M;
		E = E - f / ( 1.0 - e*cos(E) );
		i = i+1;
		if (i==maxit)
			printf("convergence problems in EccAnom");
    }  
	return E;
}