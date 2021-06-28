#include "../include/LTC.h"
#include "../include/arrays.h"
#include "../include/R_y.h"
#include "../include/R_z.h"

//--------------------------------------------------------------------------
// 
// LTC.m
//
// Purpose:
//   Transformation from Greenwich meridian system to 
//   local tangent coordinates
//
// Inputs:
//   lon      -Geodetic East longitude [rad]
//   lat      -Geodetic latitude [rad]
//   
// Output:
//   M        -Rotation matrix from the Earth equator and Greenwich meridian
//             to the local tangent (East-North-Zenith) coordinate system
//
// Last modified:   2015/08/12   M. Mahooti
//
//--------------------------------------------------------------------------

double **LTC(double lon,double lat){

	double **M,Aux;
	int j;

	M = prod(R_y(-1.0*lat),3,3,R_z(lon),3,3);
	for(j=0;j<3;++j){
		Aux=M[0][j]; 
		M[0][j]=M[1][j]; 
		M[1][j]=M[2][j]; 
		M[2][j]= Aux;
	}

	return M;
}