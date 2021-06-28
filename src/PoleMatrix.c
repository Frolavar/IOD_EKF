#include "../include/R_x.h"
#include "../include/R_y.h"
#include "../include/PoleMatrix.h"
#include "../include/arrays.h"

//--------------------------------------------------------------------------
//
// PoleMatrix.m
//
// Purpose:
//   Transformation from pseudo Earth-fixed to Earth-fixed coordinates
//   for a given date
//
// Input:
//   Pole coordinte(xp,yp)
//
// Output:
//   PoleMat   Pole matrix
//
// Last modified:   2015/08/12   M. Mahooti
// 
//--------------------------------------------------------------------------

double **PoleMatrix(double xp,double yp){
	return (prod(R_y(-xp),3,3,R_x(-yp),3,3));
}