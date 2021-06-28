#include "../include/EqnEquinox.h"
#include "../include/NutAngles.h"
#include "../include/MeanObliquity.h"
#include <math.h>

//--------------------------------------------------------------------------
//
// EqnEquinox.m
//
// Purpose:
//   Computation of the equation of the equinoxes
//
// Input:
//   Mjd_TT    Modified Julian Date (Terrestrial Time)
// 
// Output:
//    EqE      Equation of the equinoxes
//
// Notes:
//   The equation of the equinoxes dpsi*cos(eps) is the right ascension of
//   the mean equinox referred to the true equator and equinox and is equal
//   to the difference between apparent and mean sidereal time.
//
// Last modified:   2015/08/12   M. Mahooti
// 
//--------------------------------------------------------------------------

double EqnEquinox(double Mjd_TT){
	double dpsi,deps;
	
	// Nutation in longitude and obliquity
	NutAngles(&dpsi,&deps,Mjd_TT);

	// Equation of the equinoxes
	return(dpsi * cos ( MeanObliquity(Mjd_TT) ));
}