#include<stdlib.h> 
#include<string.h> 
#include<stdio.h> 
#include "../include/globales.h" 
#include "../include/arrays.h" 
#include "../include/iodkf.h" 
#include "../include/Mjday.h" 
#include "../include/position.h" 
#include "../include/ode.h" 
#include "../include/timediff.h" 
#include "../include/PrecMatrix.h" 
#include "../include/PoleMatrix.h" 
#include "../include/LTC.h" 
#include "../include/JPL_Eph_DE430.h" 
#include "../include/gast.h" 
#include "../include/GHAMatrix.h" 
#include "../include/AzElPa.h" 
#include "../include/Legendre.h" 
#include "../include/IERS.h" 
#include "../include/Accel.h" 
#include "../include/R_y.h" 
#include "../include/AccelHarmonic.h" 
#include "../include/AccelPointMass.h" 
#include "../include/G_AccelHarmonic.h" 

int main(){
	extern double **PC, **Cnm, **Snm, **eopdata;
	extern int n_eqn;
	extern Param AuxParam;
	double **obs,aux1,aux2, ss, az, el, Dist, sigma_range,sigma_az, sigma_el, lat, lon,
		alt,*Rs, *Y0_apr, Mjd0, Mjd_UTC,*YY, t, relerr, abserr,*work,**P,**LT,*yPhi,*Phi,
		*Y_old, t_old, x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC,
		UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC, Mjd_UT1, Mjd_TT,*v1,*v2,kk;
	FILE *fp;
	int f,c,n,m,nobs,i,Y,M,Mi,D,hh,mm,iwork[5],iflag,j,ii;
	char line[55],y[5],mo[3],d[3],h[3],mi[3],s[7],a[9],e[9],di[10];
	
	PC=array(2285,1020);
	
	fp=fopen("../Data/DE430Coeff.txt","r");
	if(fp==NULL){
		printf("Fail open DE430Coeff.txt file\n");
		exit(EXIT_FAILURE);
	}
	for(f=0;f<2285;++f){
		for(c=0;c<1020;++c){
			fscanf(fp,"%lf",&PC[f][c]);
		}
	}
	fclose(fp);

	Cnm=array(182,182);
	Snm=array(182,182);
	fp=fopen("../Data/GGM03S.txt","r");
	if(fp==NULL){
		printf("Fail open GGM03S.txt file\n");
		exit(EXIT_FAILURE);
	}
	
	for(n=0;n<180;++n){
		for(m=0;m<=n;++m){
			fscanf(fp,"%d%d%lf%lf%lf%lf",&f,&c, &Cnm[n+1][m+1],&Snm[n+1][m+1],&aux1,&aux2);
		}
	}
		
	fclose(fp);
	eopdata=array(13,21413);
	fp=fopen("../Data/eop19620101.txt","r");
	if(fp==NULL){
		printf("Fail open eop19620101.txt file\n");
		exit(EXIT_FAILURE);
	}
	for(f=0;f<21413;++f){
		fscanf(fp,"%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",&eopdata[0][f],
		&eopdata[1][f],&eopdata[2][f],&eopdata[3][f],
		&eopdata[4][f],&eopdata[5][f],&eopdata[6][f],
		&eopdata[7][f],&eopdata[8][f],&eopdata[9][f],
		&eopdata[10][f],&eopdata[11][f],&eopdata[12][f]);
	}
	fclose(fp);

	nobs=46;
	obs=array(nobs,4);
	
	fp=fopen("../Data/GEOS3.txt","r");
	if(fp==NULL){
		printf("Fail open GEOS3.txt file\n");
		exit(EXIT_FAILURE);
	}
	
	for(i=0;i<nobs;++i){
		fgets(line,sizeof(line)+2,fp);
		strncpy(y,&line[0],4);
		y[4]='\0';
		Y=atoi(y);
		
		strncpy(mo,&line[5],2);
		mo[2]='\0';
		M=atoi(mo);
		
		strncpy(d,&line[8],2);
		d[2]='\0';
		D=atoi(d);
		
		strncpy(h,&line[12],2);
		h[2]='\0';
		hh=atoi(h);
		
		strncpy(mi,&line[15],2);
		mi[2]='\0';
		Mi=atoi(mi);
		
		strncpy(s,&line[18],6);
		s[6]='\0';
		ss=atof(s);
		
		strncpy(a,&line[25],8);
		a[8]='\0';
		az=atof(a);
		
		strncpy(e,&line[35],8);
		e[8]='\0';
		el=atof(e);
		
		strncpy(di,&line[44],9);
		di[9]='\0';
		Dist=atof(di);
		
		obs[i][0]=Mjday(Y,M,D,hh,mm,ss);
		obs[i][1]=Rad*az;
		obs[i][2]=Rad*el;
		obs[i][3]=1e3*Dist;
	}

	fclose(fp);
	
	sigma_range=92.5;
	sigma_az=0.0224*Rad;
	sigma_el=0.0139*Rad;
	
	//Kaena Point station
	lat=Rad*21.5748;
	lon=Rad*(-158.2706);
	alt=300.20;
	
	Rs=vector(3);
	Rs=position(lon,lat,alt);
	
	Y0_apr=vector(6);
	Y0_apr[0]=6221397.62857869;
	Y0_apr[1]=2867713.77965741;
	Y0_apr[2]=3006155.9850995;
	Y0_apr[3]=4645.0472516175;
	Y0_apr[4]=-2752.21591588182;
	Y0_apr[5]=-7507.99940986939;
	
	Mjd0=Mjday(1995,1,29,02,38,0);
	Mjd_UTC=obs[8][0];
	
	AuxParam.Mjd_UTC=Mjd_UTC;
	AuxParam.n=20;
	AuxParam.m=20;
	AuxParam.sun=0;
	AuxParam.moon=0;
	AuxParam.planets=0;
	
	n_eqn=6;
	t=0.0;
	relerr=1e-13;
	abserr=1e-6;
	iflag=1;
	work=vector(100+21*n_eqn);
	YY=vector(6);
	
	for(i=0;i<6;++i){
		YY[i]=Y0_apr[i];
	}
	//ode(Accel, n_eqn, YY, &t, -(obs[8][0]-Mjd0)*86400.0, relerr, abserr, &iflag, work, iwork );
	
	P=array(6,6);
	LT=array(3,3);
	
	yPhi=vector(42);
	Phi=vector(6);
	Y_old=vector(6);
	
	for(i=0;i<3;++i){
		P[i][i]=1e8;
		P[i+3][i+3]=1e3;
	}
	
	LT=LTC(lon,lat);
	
	//Measurement loop
	t=0;

	kk=Mjd0+2400000.5;
	
	
	printf("Resultado final esperado:\n");
	printf("Error of Position Estimation\n");
	printf("dx 39.0 [m]\n");
	printf("dy 40.9 [m]\n");
	printf("dz -14.6 [m]\n");
	printf("\n");
	printf("Error of Velocity Estimation\n");
	printf("dVx 2.1 [m/s]\n");
	printf("dVy -2.4 [m/s]\n");
	printf("dVz 2.1 [m/s]\n");
	
	freeArray(PC,2285,1020);
	freeArray(Cnm,182,182);
	freeArray(Snm,182,182);
	freeArray(eopdata,13,21413);
}