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

//--------------------------------------------------------------------------
//
// Accel.m
//
// Purpose:
//   Computes the acceleration of an Earth orbiting satellite due to 
//    - the Earth's harmonic gravity field, 
//    - the gravitational perturbations of the Sun and Moon
//    - the solar radiation pressure and
//    - the atmospheric drag
//
// Inputs:
//   Mjd_TT      Terrestrial Time (Modified Julian Date)
//   Y           Satellite state vector in the ICRF/EME2000 system
//
// Output:
//   dY		    Acceleration (a=d^2r/dt^2) in the ICRF/EME2000 system
//
// Last modified:   2015/08/12   M. Mahooti
// 
//--------------------------------------------------------------------------

void Accel(double x, double *Y, double **dY){
	extern Param AuxParam; 
	double x_pole,y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC, UT1_TAI, 
		UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC, Mjd_UT1, Mjd_TT, **P, **N, **T, **E, MJD_TDB, *a, *r_Mercury, 
		*r_Venus, *r_Earth, *r_Mars, *r_Jupiter, *r_Saturn, *r_Uranus, *r_Neptune, *r_Pluto, *r_Moon, *r_Sun; 
	int i;
	
	IERS(AuxParam.Mjd_UTC + x/86400.0, 'l', &x_pole,&y_pole, &UT1_UTC, &LOD, &dpsi, &deps, &dx_pole, &dy_pole, &TAI_UTC); 
	
	timediff(UT1_UTC, TAI_UTC, &UT1_TAI, &UTC_GPS, &UT1_GPS, &TT_UTC, &GPS_UTC); 
	
	Mjd_UT1 = AuxParam.Mjd_UTC + x/86400.0 + UT1_UTC/86400.0; 
	Mjd_TT = AuxParam.Mjd_UTC + x/86400.0 + TT_UTC/86400.0; 
	
	P = PrecMatrix(MJD_J2000, Mjd_TT); 
	N = NutMatrix(Mjd_TT); 
	T=prod(N,3,3,P,3,3);
	E = prod(prod(PoleMatrix(x_pole,y_pole),3,3,GHAMatrix(Mjd_UT1),3,3),3,3,T,3,3);
	
	MJD_TDB = Mjday_TDB(Mjd_TT); 
	JPL_Eph_DE430(&r_Mercury, &r_Venus, &r_Earth, &r_Mars, &r_Jupiter, &r_Saturn, &r_Uranus, &r_Neptune, &r_Pluto, &r_Moon, &r_Sun,MJD_TDB);
	
	//Acceleration due to harmonic gravity field
	a=AccelHarmonic(Y,E,3,AuxParam.n,AuxParam.m);
	if(AuxParam.sun)
		a=sumV(a,3,AccelPointMass(Y,3,r_Sun,3,GM_Sun),3);
	if(AuxParam.moon)
		a=sumV(a,3,AccelPointMass(Y,3,r_Moon,3,GM_Moon),3);
	if(AuxParam.planets){
		a=sumV(a,3,AccelPointMass(Y,3,r_Mercury,3,GM_Mercury),3);
		a=sumV(a,3,AccelPointMass(Y,3,r_Venus,3,GM_Venus),3);
		a=sumV(a,3,AccelPointMass(Y,3,r_Mars,3,GM_Mars),3);
		a=sumV(a,3,AccelPointMass(Y,3,r_Jupiter,3,GM_Jupiter),3);
		a=sumV(a,3,AccelPointMass(Y,3,r_Saturn,3,GM_Saturn),3);
		a=sumV(a,3,AccelPointMass(Y,3,r_Uranus,3,GM_Uranus),3);
		a=sumV(a,3,AccelPointMass(Y,3,r_Neptune,3,GM_Neptune),3);
		a=sumV(a,3,AccelPointMass(Y,3,r_Pluto,3,GM_Pluto),3);
	}
	double *aux=vector(6);
	for(i=0;i<3;++i)
		aux[i]=Y[3+i];
	for(i=0;i<3;++i)
		aux[i+3]=a[i];
	*dY=aux;
}