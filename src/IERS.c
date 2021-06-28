#include "../include/IERS.h"
#include "../include/iodkf.h"
#include "../include/arrays.h"
#include <math.h>
#include <stdio.h>

//--------------------------------------------------------------------------
//
// IERS: Management of IERS time and polar motion data
//  
// Last modified:   2018/02/01   M. Mahooti
// 
//--------------------------------------------------------------------------

void IERS(double Mjd_UTC,char interp, double *x_pole,double *y_pole, double *UT1_UTC,double *LOD,double *dpsi,double *deps, double *dx_pole,double *dy_pole,double *TAI_UTC){
	extern double **eopdata;
	double mjd,*vv,*preeop,*nexteop,fixf,mfme, *eop;
	int i,j;
	if(interp=='l'){
		//linear interpolation
		mjd=(floor(Mjd_UTC));
		vv=vector(21413);
		for(i=0;i<21413;++i){
			vv[i]=eopdata[3][i];
		}
		
		i=find2(vv,21413,mjd);
		freeVector(vv,21413);
		preeop=vector(13);
		nexteop=vector(13);
		
		for(j=0;j<13;++j){
			preeop[j]=eopdata[j][i];
			nexteop[j]=eopdata[j][i+1];
		}
		mfme=1440*(Mjd_UTC-floor(Mjd_UTC));
		fixf=mfme/1440.0;
		
		// Setting of IERS Earth rotation parameters
		// (UT1-UTC [s], TAI-UTC [s], x ["], y ["])
		*x_pole  = preeop[4]+(nexteop[4]-preeop[4])*fixf;
		*y_pole  = preeop[5]+(nexteop[5]-preeop[5])*fixf;
		*UT1_UTC = preeop[6]+(nexteop[6]-preeop[6])*fixf;
		*LOD     = preeop[7]+(nexteop[7]-preeop[7])*fixf;
		*dpsi    = preeop[8]+(nexteop[8]-preeop[8])*fixf;
		*deps    = preeop[9]+(nexteop[9]-preeop[9])*fixf;
		*dx_pole = preeop[10]+(nexteop[10]-preeop[10])*fixf;
		*dy_pole = preeop[11]+(nexteop[11]-preeop[11])*fixf;
		*TAI_UTC = preeop[12];
		
		*x_pole  = *x_pole/(Arcs);  // Pole coordinate [rad]
		*y_pole  = *y_pole/(Arcs);  // Pole coordinate [rad]
		*dpsi    = *dpsi/(Arcs);
		*deps    = *deps/(Arcs);
		*dx_pole = *dx_pole/(Arcs); // Pole coordinate [rad]
		*dy_pole = *dy_pole/(Arcs); // Pole coordinate [rad]
		
		freeVector(preeop,13);
		freeVector(nexteop,13);
	}else if(interp=='n'){
		mjd = (floor(Mjd_UTC));
		vv=vector(21413);
		for(i=0;i<21413;++i){
			vv[i]=eopdata[3][i];
		}
		i=find2(vv,21413,mjd);
		freeVector(vv,21413);
		
		nexteop=vector(13);
		
		for(j=0;j<13;++j){
			eop[j]=eopdata[j][i];
		}
		
		// Setting of IERS Earth rotation parameters
		// (UT1-UTC [s], TAI-UTC [s], x ["], y ["])
		*x_pole  = eop[4]/(Arcs);  // Pole coordinate [rad]
		*y_pole  = eop[5]/(Arcs);  // Pole coordinate [rad]
		*UT1_UTC = eop[6];             // UT1-UTC time difference [s]
		*LOD     = eop[7];             // Length of day [s]
		*dpsi    = eop[8]/(Arcs);
		*deps    = eop[9]/(Arcs);
		*dx_pole = eop[10]/(Arcs); // Pole coordinate [rad]
		*dy_pole = eop[11]/(Arcs); // Pole coordinate [rad]
		*TAI_UTC = eop[12];            // TAI-UTC time difference [s]
		freeVector(eop,13);
	}
}