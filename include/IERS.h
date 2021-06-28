#ifndef IERS_h_
#define IERS_h_

void IERS(double Mjd_UTC,char interp, double *x_pole,double *y_pole, double *UT1_UTC,double *LOD,double *dpsi,double *deps, double *dx_pole,double *dy_pole,double *TAI_UTC);

#endif