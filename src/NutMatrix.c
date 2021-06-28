#include "../include/R_x.h"
#include "../include/R_z.h"
#include "../include/NutAngles.h"
#include "../include/MeanObliquity.h"
#include "../include/NutMatrix.h"
#include "../include/arrays.h"

//--------------------------------------------------------------------------
//
// NutMatrix.m
//
// Purpose:
//   Transformation from mean to true equator and equinox
//
// Input:
//   Mjd_TT    Modified Julian Date (Terrestrial Time)
//
// Output:
//   NutMat    Nutation matrix
//
// Last modified:   2015/08/12   M. Mahooti
// 
//--------------------------------------------------------------------------

double **NutMatrix(double Mjd_TT){
	double eps,dpsi,deps;

	// Mean obliquity of the ecliptic
	eps = MeanObliquity (Mjd_TT);

	// Nutation in longitude and obliquity
	NutAngles(&dpsi,&deps,Mjd_TT);

	// Transformation from mean to true equator and equinox
	return (prod(prod(R_x(-eps-deps),3,3,R_z(-dpsi),3,3),3,3,R_x(+eps),3,3));
}