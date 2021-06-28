#include "../include/GHAMatrix.h"
#include "../include/gast.h"
#include "../include/R_z.h"

//--------------------------------------------------------------------------
//
// GHAMatrix.m
//
// Purpose:
//   Transformation from true equator and equinox to Earth equator and 
//   Greenwich meridian system 
//
// Input:
//   Mjd_UT1   Modified Julian Date UT1
// 
// Output:
//   GHAmat    Greenwich Hour Angle matrix
//
// Last modified:   2015/08/12   M. Mahooti
// 
//--------------------------------------------------------------------------
double **GHAMatrix(double Mjd_UT1){
	return R_z( gast(Mjd_UT1) );
}