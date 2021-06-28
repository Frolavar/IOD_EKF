#include "../include/AzElPa.h"
#include "../include/iodkf.h"
#include "../include/arrays.h"
#include <math.h>
#include <stdio.h>
//--------------------------------------------------------------------------
//
// Purpose:
//  Computes azimuth, elevation and partials from local tangent coordinates
//
// Input:
//   s      Topocentric local tangent coordinates (East-North-Zenith frame)
// 
// Outputs:
//   A      Azimuth [rad]
//   E      Elevation [rad]
//   dAds   Partials of azimuth w.r.t. s
//   dEds   Partials of elevation w.r.t. s
//
// Last modified:   2015/08/12   M. Mahooti
//
//--------------------------------------------------------------------------

void AzElPa(double *s,int n, double *Az, double *El,double **dAds, int n1,double **dEds, int n2){

	double rho;
	int i;

	rho = sqrt(s[0]*s[0]+s[1]*s[1]);

	// Angles
	*Az = atan2(s[0],s[1]);

	if (*Az<0.0){ 
		*Az = *Az+pi2;
	}
	
	*El = atan( s[2] / rho );

	// Partials
	(*dAds)[0]=s[1]/(rho*rho);
	(*dAds)[1]=-s[0]/(rho*rho);
	(*dAds)[2]=0.0;
	(*dEds)[0]=-s[0]*s[2]/rho;
	(*dEds)[1]=-s[1]*s[2]/rho;
	(*dEds)[2]=rho;
	*dEds=esc_x_vec(1/dot(s,n,s,n),*dEds,n2);
}