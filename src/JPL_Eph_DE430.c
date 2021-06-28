#include "../include/JPL_Eph_DE430.h"
#include "../include/arrays.h"
#include "../include/Cheb3D.h"
#include <stdio.h>

//--------------------------------------------------------------------------
//
// JPL_Eph_DE430: Computes the sun, moon, and nine major planets' equatorial
//                position using JPL Ephemerides
//
// Inputs:
//   Mjd_TDB         Modified julian date of TDB
//
// Output:
//   r_Earth(solar system barycenter (SSB)),r_Mars,r_Mercury,r_Venus,
//   r_Jupiter,r_Saturn,r_Uranus,r_Neptune,r_Pluto,r_Moon,
//   r_Sun(geocentric equatorial position ([m]) referred to the
//   International Celestial Reference Frame (ICRF))
//
// Notes: Light-time is already taken into account
//
// Last modified:   2018/01/11   M. Mahooti
// 
//--------------------------------------------------------------------------



void JPL_Eph_DE430(double **r_Mercury,double **r_Venus,double **r_Earth,double **r_Mars,double **r_Jupiter,double **r_Saturn,double **r_Uranus,double **r_Neptune,double **r_Pluto,double **r_Moon,double **r_Sun,double Mjd_TDB){
		
	extern double **PC;
	double *v1,*v2,JD,t1,dt, *PCtemp,*Cx_Earth,*Cy_Earth,*Cz_Earth,*Cx,*Cy,*Cz,Mjd0,*Cx_Moon,*Cy_Moon,*Cz_Moon,*Cx_Sun,*Cy_Sun,*Cz_Sun,*Cx_Mercury,*Cy_Mercury,*Cz_Mercury,*Cx_Venus,*Cy_Venus,*Cz_Venus,
			*Cx_Mars,*Cy_Mars,*Cz_Mars,*Cx_Jupiter,*Cy_Jupiter,*Cz_Jupiter,*Cx_Saturn,*Cy_Saturn,*Cz_Saturn,*Cx_Uranus,*Cy_Uranus,*Cz_Uranus,*Cx_Neptune,*Cy_Neptune,*Cz_Neptune,
			*Cx_Pluto,*Cy_Pluto,*Cz_Pluto,*Cx_Nutations,*Cy_Nutations,*zeros,*Cx_Librations,*Cy_Librations,*Cz_Librations,*Nutations,*Librations,EMRAT,EMRAT1;
	int j,i,temp[4];

	JD=Mjd_TDB + 2400000.5;

	v1=vector(2285);
	v2=vector(2285);

	for(j=0;j<2285;++j){
		v1[j]=PC[j][0];
		v2[j]=PC[j][1];
	}
	
	i=find1(v1,2285,v2,2285,JD);

	freeVector(v1,2285);
	freeVector(v2,2285);
	
	PCtemp= vector(1020);
	
	for(j=0;j<1020;++j){
		PCtemp[j]=PC[i][j];
	}
	t1=PCtemp[0]-2400000.5;
	dt=Mjd_TDB-t1;

	//Earth
	for(j=0;j<4;++j){
		temp[j]=230+13*j;
	}
	
	Cx_Earth=vector(26);
	Cy_Earth=vector(26);
	Cz_Earth=vector(26);
	
	for(j=0;j<13;++j){
		Cx_Earth[j]=PCtemp[temp[0]+j];
		Cy_Earth[j]=PCtemp[temp[1]+j];
		Cz_Earth[j]=PCtemp[temp[2]+j];
	}
	
	for(j=0;j<4;++j){
		temp[j]+=39;
	}
	
	Cx=vector(14);
	Cy=vector(14);
	Cz=vector(14);
	
	for(j=0;j<13;++j){
		Cx[j]=PCtemp[temp[0]+j];
		Cy[j]=PCtemp[temp[1]+j];
		Cz[j]=PCtemp[temp[2]+j];
	}

	for(j=0;j<13;++j){
		Cx_Earth[j+13]=Cx[j];
		Cy_Earth[j+13]=Cy[j];
		Cz_Earth[j+13]=Cz[j];
	}
	
	if (0<=dt && dt<=16){
		j=0;
		Mjd0 = t1;
	}else if(16<dt && dt<=32){
		j=1;
		Mjd0 = t1+16*j;
	}
	
	*r_Earth=Cheb3D(Mjd_TDB,13,Mjd0,Mjd0+16,&Cx_Earth[13*j],&Cy_Earth[13*j],&Cz_Earth[13*j]);

	*r_Earth=esc_x_vec(1e3,*r_Earth,3);

	freeVector(Cx_Earth,26);
	freeVector(Cy_Earth,26);
	freeVector(Cz_Earth,26);
	
	//Moon
	for(j=0;j<4;++j){
		temp[j]=440+13*j;
	}
	Cx_Moon=vector(109);
	Cy_Moon=vector(109);
	Cz_Moon=vector(109);
	
	for(j=0;j<13;++j){
		Cx_Moon[j]=PCtemp[temp[0]+j];
		Cy_Moon[j]=PCtemp[temp[1]+j];
		Cz_Moon[j]=PCtemp[temp[2]+j];
	}
	
	for(i=0;i<7;++i){
		for(j=0;j<4;++j)
			temp[j]+=39;
		for(j=0;j<13;++j){
			Cx[j]=PCtemp[temp[0]+j];
			Cy[j]=PCtemp[temp[1]+j];
			Cz[j]=PCtemp[temp[2]+j];
		}
		for(j=0;j<13;++j){
			Cx_Moon[13*(i+1)+j]=Cx[j];
			Cy_Moon[13*(i+1)+j]=Cy[j];
			Cz_Moon[13*(i+1)+j]=Cz[j];
		}
	}
	if ((0<=dt) && (dt<=4)){
		j=0;
		Mjd0 = t1;
	}else if((4<dt) && (dt<=8)){
		j=1;
		Mjd0 = t1+4*j;
	}else if((8<dt) && (dt<=12)){
		j=2;
		Mjd0 = t1+4*j;
	}else if((12<dt) && (dt<=16)){
		j=3;
		Mjd0 = t1+4*j;
	}else if((16<dt) && (dt<=20)){
		j=4;
		Mjd0 = t1+4*j;
	}else if((20<dt) && (dt<=24)){
		j=5;
		Mjd0 = t1+4*j;
	}else if((24<dt) && (dt<=28)){
		j=6;
		Mjd0 = t1+4*j;
	}else if((28<dt) && (dt<=32)){
		j=7;
		Mjd0 = t1+4*j;
	}
	
	*r_Moon=Cheb3D(Mjd_TDB,13,Mjd0,Mjd0+4,&Cx_Moon[13*j],&Cy_Moon[13*j],&Cz_Moon[13*j]);
	*r_Moon=esc_x_vec(1e3,*r_Moon,3);
	
	freeVector(Cx_Moon,109);
	freeVector(Cy_Moon,109);
	freeVector(Cz_Moon,109);
	
	//Sun
	for(j=0;j<4;++j){
		temp[j]=752+11*j;
	}
	
	Cx_Sun=vector(22);
	Cy_Sun=vector(22);
	Cz_Sun=vector(22);	
	
	for(j=0;j<11;++j){
		Cx_Sun[j]=PCtemp[temp[0]+j];
		Cy_Sun[j]=PCtemp[temp[1]+j];
		Cz_Sun[j]=PCtemp[temp[2]+j];
	}
	for(j=0;j<4;++j){
		temp[j]+=33;
	}
	
	for(j=0;j<11;++j){
		Cx[j]=PCtemp[temp[0]+j];
		Cy[j]=PCtemp[temp[1]+j];
		Cz[j]=PCtemp[temp[2]+j];
	}
	
	for(j=0;j<11;++j){
		Cx_Sun[j+11]=Cx[j];
		Cy_Sun[j+11]=Cy[j];
		Cz_Sun[j+11]=Cz[j];
	}
	if (0<=dt && dt<=16){
		j=0;
		Mjd0 = t1;
	}else if(16<dt && dt<=32){
		j=1;
		Mjd0 = t1+16*j;
	}
	
	*r_Sun=Cheb3D(Mjd_TDB,11,Mjd0,Mjd0+16,&Cx_Sun[11*j],&Cy_Sun[11*j],&Cz_Sun[11*j]);
	*r_Sun=esc_x_vec(1e3,*r_Sun,3);
	
	freeVector(Cx_Sun,22);
	freeVector(Cy_Sun,22);
	freeVector(Cz_Sun,22);
	
	//Mercury
	for(j=0;j<4;++j){
		temp[j]=2+14*j;
	}
	
	Cx_Mercury=vector(56);
	Cy_Mercury=vector(56);
	Cz_Mercury=vector(56);	
	
	for(j=0;j<14;++j){
		Cx_Mercury[j]=PCtemp[temp[0]+j];
		Cy_Mercury[j]=PCtemp[temp[1]+j];
		Cz_Mercury[j]=PCtemp[temp[2]+j];
	}
	
	for(i=0;i<3;++i){
		for(j=0;j<4;++j)
			temp[j]+=42;
		for(j=0;j<14;++j){
			Cx[j]=PCtemp[temp[0]+j];
			Cy[j]=PCtemp[temp[1]+j];
			Cz[j]=PCtemp[temp[2]+j];
		}
		for(j=0;j<14;++j){
			Cx_Mercury[14*(i+1)+j]=Cx[j];
			Cy_Mercury[14*(i+1)+j]=Cy[j];
			Cz_Mercury[14*(i+1)+j]=Cz[j];
		}
	}
	
	if (0<=dt && dt<=8){
		j=0;
		Mjd0 = t1;
	}else if(8<dt && dt<=16){
		j=1;
		Mjd0 = t1+8*j;
	}else if (16<dt && dt<=24){
		j=2;
		Mjd0 = t1+8*j;
	}else if(24<dt && dt<=32){
		j=3;
		Mjd0 = t1+8*j;
	}
	
	*r_Mercury=Cheb3D(Mjd_TDB,14,Mjd0,Mjd0+8,&Cx_Mercury[14*j],&Cy_Mercury[14*j],&Cz_Mercury[14*j]);
	
	*r_Mercury=esc_x_vec(1e3,*r_Mercury,3);
	
	freeVector(Cx_Mercury,56);
	freeVector(Cy_Mercury,56);
	freeVector(Cz_Mercury,56);
	
	//Venus
	for(j=0;j<4;++j){
		temp[j]=170+10*j;
	}
	
	Cx_Venus=vector(20);
	Cy_Venus=vector(20);
	Cz_Venus=vector(20);	
	
	for(j=0;j<10;++j){
		Cx_Venus[j]=PCtemp[temp[0]+j];
		Cy_Venus[j]=PCtemp[temp[1]+j];
		Cz_Venus[j]=PCtemp[temp[2]+j];
	}
	for(j=0;j<4;++j){
		temp[j]+=30;
	}
	for(j=0;j<10;++j){
		Cx[j]=PCtemp[temp[0]+j];
		Cy[j]=PCtemp[temp[1]+j];
		Cz[j]=PCtemp[temp[2]+j];
	}
	for(j=0;j<10;++j){
		Cx_Venus[j+10]=Cx[j];
		Cy_Venus[j+10]=Cy[j];
		Cz_Venus[j+10]=Cz[j];
	}
	if (0<=dt && dt<=16){
		j=0;
		Mjd0 = t1;
	}else if(16<dt && dt<=32){
		j=1;
		Mjd0 = t1+16*j;
	}
					 
	*r_Venus=Cheb3D(Mjd_TDB,10,Mjd0,Mjd0+16,&Cx_Venus[10*j],&Cy_Venus[10*j],&Cz_Venus[10*j]);
	*r_Venus=esc_x_vec(1e3,*r_Venus,3);
	
	freeVector(Cx_Venus,56);
	freeVector(Cy_Venus,56);
	freeVector(Cz_Venus,56);
	
	//Mars
	for(j=0;j<4;++j){
		temp[j]=308+11*j;
	}
	
	Cx_Mars=vector(11);
	Cy_Mars=vector(11);
	Cz_Mars=vector(11);	
	
	for(j=0;j<11;++j){
		Cx_Mars[j]=PCtemp[temp[0]+j];
		Cy_Mars[j]=PCtemp[temp[1]+j];
		Cz_Mars[j]=PCtemp[temp[2]+j];
	}
	j=0;
	Mjd0=t1;
	
	*r_Mars=Cheb3D(Mjd_TDB,11,Mjd0,Mjd0+32,&Cx_Mars[11*j],&Cy_Mars[11*j],&Cz_Mars[11*j]);
	*r_Mars=esc_x_vec(1e3,*r_Mars,3);
	
	freeVector(Cx_Mars,11);
	freeVector(Cy_Mars,11);
	freeVector(Cz_Mars,11);
	
	//Jupiter
	for(j=0;j<4;++j){
		temp[j]=341+8*j;
	}
	
	Cx_Jupiter=vector(8);
	Cy_Jupiter=vector(8);
	Cz_Jupiter=vector(8);	
	
	for(j=0;j<8;++j){
		Cx_Jupiter[j]=PCtemp[temp[0]+j];
		Cy_Jupiter[j]=PCtemp[temp[1]+j];
		Cz_Jupiter[j]=PCtemp[temp[2]+j];
	}
	j=0;
	Mjd0=t1;
	
	*r_Jupiter=Cheb3D(Mjd_TDB,8,Mjd0,Mjd0+32,&Cx_Jupiter[8*j],&Cy_Jupiter[8*j],&Cz_Jupiter[8*j]);
	*r_Jupiter=esc_x_vec(1e3,*r_Jupiter,3);
	
	freeVector(Cx_Jupiter,8);
	freeVector(Cy_Jupiter,8);
	freeVector(Cz_Jupiter,8);
	
	//Saturn
	for(j=0;j<4;++j){
		temp[j]=365+7*j;
	}
	
	Cx_Saturn=vector(7);
	Cy_Saturn=vector(7);
	Cz_Saturn=vector(7);	
	
	for(j=0;j<7;++j){
		Cx_Saturn[j]=PCtemp[temp[0]+j];
		Cy_Saturn[j]=PCtemp[temp[1]+j];
		Cz_Saturn[j]=PCtemp[temp[2]+j];
	}
	j=0;
	Mjd0=t1;
	
	*r_Saturn=Cheb3D(Mjd_TDB,7,Mjd0,Mjd0+32,&Cx_Saturn[7*j],&Cy_Saturn[7*j],&Cz_Saturn[7*j]);
	*r_Saturn=esc_x_vec(1e3,*r_Saturn,3);
	
	freeVector(Cx_Saturn,7);
	freeVector(Cy_Saturn,7);
	freeVector(Cz_Saturn,7);
	
	//Uranus
	for(j=0;j<4;++j){
		temp[j]=386+6*j;
	}
	
	Cx_Uranus=vector(6);
	Cy_Uranus=vector(6);
	Cz_Uranus=vector(6);	
	
	for(j=0;j<6;++j){
		Cx_Uranus[j]=PCtemp[temp[0]+j];
		Cy_Uranus[j]=PCtemp[temp[1]+j];
		Cz_Uranus[j]=PCtemp[temp[2]+j];
	}
	j=0;
	Mjd0=t1;
	
	*r_Uranus=Cheb3D(Mjd_TDB,6,Mjd0,Mjd0+32,&Cx_Uranus[6*j],&Cy_Uranus[6*j],&Cz_Uranus[6*j]);
	*r_Uranus=esc_x_vec(1e3,*r_Uranus,3);
	
	freeVector(Cx_Uranus,6);
	freeVector(Cy_Uranus,6);
	freeVector(Cz_Uranus,6);
	
	//Neptune
	for(j=0;j<4;++j){
		temp[j]=404+6*j;
	}
	
	Cx_Neptune=vector(6);
	Cy_Neptune=vector(6);
	Cz_Neptune=vector(6);	
	
	for(j=0;j<6;++j){
		Cx_Neptune[j]=PCtemp[temp[0]+j];
		Cy_Neptune[j]=PCtemp[temp[1]+j];
		Cz_Neptune[j]=PCtemp[temp[2]+j];
	}
	j=0;
	Mjd0=t1;
	
	*r_Neptune=Cheb3D(Mjd_TDB,6,Mjd0,Mjd0+32,&Cx_Neptune[6*j],&Cy_Neptune[6*j],&Cz_Neptune[6*j]);
	*r_Neptune=esc_x_vec(1e3,*r_Neptune,3);
	
	freeVector(Cx_Neptune,6);
	freeVector(Cy_Neptune,6);
	freeVector(Cz_Neptune,6);
	
	//Pluto
	for(j=0;j<4;++j){
		temp[j]=422+6*j;
	}
	
	Cx_Pluto=vector(6);
	Cy_Pluto=vector(6);
	Cz_Pluto=vector(6);	
	
	for(j=0;j<6;++j){
		Cx_Pluto[j]=PCtemp[temp[0]+j];
		Cy_Pluto[j]=PCtemp[temp[1]+j];
		Cz_Pluto[j]=PCtemp[temp[2]+j];
	}
	j=0;
	Mjd0=t1;
	
	*r_Pluto=Cheb3D(Mjd_TDB,6,Mjd0,Mjd0+32,&Cx_Pluto[6*j],&Cy_Pluto[6*j],&Cz_Pluto[6*j]);
	*r_Pluto=esc_x_vec(1e3,*r_Pluto,3);
	
	freeVector(Cx_Pluto,6);
	freeVector(Cy_Pluto,6);
	freeVector(Cz_Pluto,6);
	
	//Nutations
	for(j=0;j<3;++j){
		temp[j]=818+10*j;
	}
	
	Cx_Nutations=vector(40);
	Cy_Nutations=vector(40);
	zeros=vector(10);	
	
	for(j=0;j<10;++j){
		zeros[j]=0.0;
	}
	
	for(j=0;j<10;++j){
		Cx_Nutations[j]=PCtemp[temp[0]+j];
		Cy_Nutations[j]=PCtemp[temp[1]+j];
	}
	
	for(i=0;i<3;++i){
		for(j=0;j<4;++j)
			temp[j]+=20;
		for(j=0;j<10;++j){
			Cx[j]=PCtemp[temp[0]+j];
			Cy[j]=PCtemp[temp[1]+j];
		}
		for(j=0;j<10;++j){
			Cx_Nutations[10*(i+1)+j]=Cx[j];
			Cy_Nutations[10*(i+1)+j]=Cy[j];
		}
	}
	
	if (0<=dt && dt<=8){
		j=0;
		Mjd0 = t1;
	}else if(8<dt && dt<=16){
		j=1;
		Mjd0 = t1+8*j;
	}else if (16<dt && dt<=24){
		j=2;
		Mjd0 = t1+8*j;
	}else if(24<dt && dt<=32){
		j=3;
		Mjd0 = t1+8*j;
	}
	
	Nutations=Cheb3D(Mjd_TDB,10,Mjd0,Mjd0+8,&Cx_Nutations[10*j],&Cy_Nutations[10*j],zeros);
	
	freeVector(Cx_Nutations,40);
	freeVector(Cy_Nutations,40);
	freeVector(zeros,10);
	
	//Librations
	for(j=0;j<4;++j){
		temp[j]=898+10*j;
	}
	
	Cx_Librations=vector(40);
	Cy_Librations=vector(40);
	Cz_Librations=vector(40);	
	
	for(j=0;j<10;++j){
		Cx_Librations[j]=PCtemp[temp[0]+j];
		Cy_Librations[j]=PCtemp[temp[1]+j];
		Cz_Librations[j]=PCtemp[temp[2]+j];
	}
	
	for(i=0;i<3;++i){
		for(j=0;j<4;++j)
			temp[j]+=30;
		for(j=0;j<10;++j){
			Cx[j]=PCtemp[temp[0]+j];
			Cy[j]=PCtemp[temp[1]+j];
			Cz[j]=PCtemp[temp[2]+j];
		}
		for(j=0;j<10;++j){
			Cx_Librations[10*(i+1)+j]=Cx[j];
			Cy_Librations[10*(i+1)+j]=Cy[j];
			Cz_Librations[10*(i+1)+j]=Cz[j];
		}
	}
	
	if (0<=dt && dt<=8){
		j=0;
		Mjd0 = t1;
	}else if(8<dt && dt<=16){
		j=1;
		Mjd0 = t1+8*j;
	}else if (16<dt && dt<=24){
		j=2;
		Mjd0 = t1+8*j;
	}else if(24<dt && dt<=32){
		j=3;
		Mjd0 = t1+8*j;
	}
	
	Librations=Cheb3D(Mjd_TDB,10,Mjd0,Mjd0+8,&Cx_Librations[10*j],&Cy_Librations[10*j],&Cz_Librations[10*j]);
	
	freeVector(Cx_Librations,56);
	freeVector(Cy_Librations,56);
	freeVector(Cz_Librations,56);
	
	EMRAT = 81.30056907419062;
	EMRAT1 = 1/(1+EMRAT);

	*r_Earth = sumV(*r_Earth,3,esc_x_vec(-1,esc_x_vec(EMRAT1,*r_Moon,3),3),3);
	*r_Mercury = sumV(esc_x_vec(-1,*r_Earth,3),3,*r_Mercury,3);
	*r_Venus = sumV(esc_x_vec(-1,*r_Earth,3),3,*r_Venus,3);
	*r_Mars = sumV(esc_x_vec(-1,*r_Earth,3),3,*r_Mars,3);
	*r_Jupiter = sumV(esc_x_vec(-1,*r_Earth,3),3,*r_Jupiter,3);
	*r_Saturn = sumV(esc_x_vec(-1,*r_Earth,3),3,*r_Saturn,3);
	*r_Uranus = sumV(esc_x_vec(-1,*r_Earth,3),3,*r_Uranus,3);
	*r_Neptune = sumV(esc_x_vec(-1,*r_Earth,3),3,*r_Neptune,3);
	*r_Pluto = sumV(esc_x_vec(-1,*r_Earth,3),3,*r_Pluto,3);
	*r_Sun = sumV(esc_x_vec(-1,*r_Earth,3),3,*r_Sun,3);
	
	
}