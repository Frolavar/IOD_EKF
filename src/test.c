#include <float.h>
#include<stdlib.h> 
#include<string.h> 
#include<stdio.h> 

#include "../include/AccelPointMass.h" 
#include "../include/Legendre.h" 
#include "../include/globales.h" 
#include "../include/arrays.h" 
#include "../include/iodkf.h" 
#include "../include/Mjday.h" 
#include "../include/Mjday_TDB.h" 
#include "../include/position.h" 
#include "../include/ode.h" 
#include "../include/sign_.h"
#include "../include/timediff.h" 
#include "../include/PrecMatrix.h" 
#include "../include/NutMatrix.h" 
#include "../include/NutAngles.h" 
#include "../include/MeanObliquity.h" 
#include "../include/PoleMatrix.h" 
#include "../include/Cheb3D.h" 
#include "../include/LTC.h" 
#include "../include/JPL_Eph_DE430.h" 
#include "../include/gast.h" 
#include "../include/gmst.h" 
#include "../include/Frac.h" 
#include "../include/GHAMatrix.h" 
#include "../include/EqnEquinox.h" 
#include "../include/AzElPa.h" 
#include "../include/IERS.h" 
#include "../include/Accel.h" 
#include "../include/R_x.h" 
#include "../include/R_y.h" 
#include "../include/R_z.h" 
#include "../include/AccelHarmonic.h" 
#include "../include/G_AccelHarmonic.h" 
#include "../include/TimeUpdate.h" 
#include "../include/VarEqn.h" 


#define FAIL() printf("\nFallo en %s() linea %d\n", __func__, __LINE__)
#define _assert(test) do { if (!(test)) { FAIL(); return 1; } } while(0)
#define _verify(test) do { int r=test(); tests_run++; if(r) return r; } while(0)

int tests_run = 0;

int norma_test(){
	double v[]={1,-2,3};
	double n=norma(v,3);
	_assert(compareDouble(n,3.741657));
	return 0;
}

int dot_test(){
	double v[]={4,-1,2};
	double v2[]={2,-2,-1};
	double n=dot(v,3,v2,3);
	_assert(compareDouble(n,8.0));
	return 0;
}

int vector_test(){
	double v[]={0,0,0};
	double *v2=vector(3);
	_assert(compareV(v,3,v2,3));
	return 0;
}

int sumV_test(){
	double v[]={1,2,3};
	double v2[]={5,6,7};
	double v3[]={6,8,10};
	_assert(compareV(v3,3,sumV(v,3,v2,3),3));
	return 0;
}

int esc_x_vec_test(){
	double v[]={1,2,3};
	double x=3.0;
	double v3[]={3,6,9};
	_assert(compareV(v3,3,esc_x_vec(x,v,3),3));
	return 0;
}

int static_dinamic_test(){
	double v[2][2]={{-3,5},{7,10}};
	double **m=array(2,2);
	m[0][0]=-3;
	m[0][1]=5;
	m[1][0]=7;
	m[1][1]=10;
	_assert(compare(m,2,2,static_dinamic(2,2,v),2,2));
	return 0;
}

int mat_x_vec_test(){
	double v[2][3]={{-3,5,-6},{7,10,-1}};
	double v2[]={-6,-2,5};
	double v3[]={-22,-67};
	_assert(compareV(v3,2,mat_x_vec(static_dinamic(2,3,v),2,3,v2,3),2));
	return 0;
}

int find2_test(){
	double v[]={1,2,3};
	_assert(find2(v,3,2)==1);
	
	return 0;
}

int find1_test(){
	double v[]={1,2,3};
	double v2[]={1,3,2};
	_assert(find1(v,3,v2,3,2)==1);
	return 0;
}

int trasp_test(){
	double v[3][3]={{1,2,3},{4,5,6},{7,8,9}};
	double v2[3][3]={{1,4,7},{2,5,8},{3,6,9}};
	_assert(compare(static_dinamic(3,3,v2),3,3,trasp(static_dinamic(3,3,v),3),3,3));
	return 0;
}

int inv_test(){
	double v[3][3]={{2,3,1},{1,-1,2},{0,1,0}};
	double v2[3][3]={{0.66,-0.33,-2.33},{0,0,1},{-0.33,0.66,1.66}};
	_assert(compare(static_dinamic(3,3,v2),3,3,inv(static_dinamic(3,3,v),3),3,3));
	
	return 0;
}

int prod_test(){
	double v[3][3]={{1,4,7},{2,5,8},{3,6,9}};
	double v2[3][3]={{1,-1,2},{2,-1,2},{3,-3,0}};
	double v3[3][3]={{30,-26,10},{36,-31,14},{42,-36,18}};
	_assert(compare(static_dinamic(3,3,v3),3,3,prod(static_dinamic(3,3,v),3,3,static_dinamic(3,3,v2),3,3),3,3));
	return 0;
}

int eye_test(){
	double v[3][3]={{1,0,0},{0,1,0},{0,0,1}};
	_assert(compare(static_dinamic(3,3,v),3,3,eye(3),3,3));
	
	return 0;
}

int sum_test(){
	double v[2][2]={{1,4},{2,5}};
	double v2[2][2]={{1,-1},{2,-1}};
	double v3[2][2]={{2,3},{4,4}};
	_assert(compare(static_dinamic(2,2,v3),2,2,sum(static_dinamic(2,2,v),2,2,static_dinamic(2,2,v2),2,2),2,2));
	return 0;
}

int array_test(){
	double v[2][2]={{0,0},{0,0}};
	double **v2=array(2,2);
	_assert(compare(static_dinamic(2,2,v),2,2,v2,2,2));
	return 0;
}

int IERS_test(){
	double x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC;
	IERS(4.974611128472207e+04,'l',&x_pole,&y_pole,&UT1_UTC,&LOD,&dpsi,&deps,&dx_pole,&dy_pole,&TAI_UTC);
	    
	_assert(compareDouble(x_pole,-5.593861831521886e-07)&&compareDouble(y_pole,2.335544384403725e-06)&&compareDouble(UT1_UTC,0.325764698106523)
	&&compareDouble(LOD,0.002726686357638)&&compareDouble(dpsi,-1.168819606404215e-07)&&compareDouble(deps,-2.478816804122186e-08)
	&&compareDouble(dx_pole,-8.417641506705226e-10)&&compareDouble(dy_pole,-1.566188801213420e-09)&&compareDouble(TAI_UTC,29));
	return 0;
}

int timediff_test(){
	double UT1_UTC=0.325764698106523,TAI_UTC=29,UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC;
	timediff(UT1_UTC,TAI_UTC,&UT1_TAI,&UTC_GPS,&UT1_GPS,&TT_UTC,&GPS_UTC);
	_assert(compareDouble(UT1_TAI,-28.674235301893475)&&compareDouble(UTC_GPS,-10)&&compareDouble(UT1_GPS,-9.674235301893475)
	&&compareDouble(TT_UTC,61.184000000000000)&&compareDouble(GPS_UTC,10));
	return 0;
}

int PrecMatrix_test(){
	double Mjd_1=5.154450000000000e+04,Mjd_2=4.974611085861109e+04, **PM;
	double Res[3][3]={{0.999999279427995,0.001101012138263,4.784514219004540e-04},
	{-0.001101012138264,0.999999393885917,-2.633882773454226e-07},
	{-4.784514218979997e-04,-2.633927355266896e-07,0.999999885542077}};
	PM=PrecMatrix(Mjd_1,Mjd_2);
	_assert(compare(PM,3,3,static_dinamic(3,3,Res),3,3));
	return 0;
}

int R_x_test(){
	double angle=-0.409068868517939,**RM;
	double Res[3][3]={{1,0,0},{0,0.917491582882778,-0.397755195238548},{0,0.397755195238548,0.917491582882778}};
	RM=R_x(angle);
	_assert(compare(RM,3,3,static_dinamic(3,3,Res),3,3));
	return 0;
}

int R_y_test(){
	double angle=-0.376551295459273,**RM;
	double Res[3][3]={{0.929938305587722,0,0.367715580035218},{0,1,0},{-0.367715580035218,0,0.929938305587722}};
	RM=R_y(angle);
	_assert(compare(RM,3,3,static_dinamic(3,3,Res),3,3));
	return 0;
}

int R_z_test(){
	double angle=-2.762343079106937,**RM;
	double Res[3][3]={{-0.928942722252092,-0.370223471399199,0},{0.370223471399199,-0.928942722252092,0},{0,0,1}};
	RM=R_z(angle);
	_assert(compare(RM,3,3,static_dinamic(3,3,Res),3,3));
	return 0;
}

int MeanObliquity_test(){
	double Mjd_TT=4.974611085861109e+04,res=0.409103979363840,m;
	m=MeanObliquity(Mjd_TT);
	_assert(compareDouble(res,m));
	return 0;
}

int NutAngles_test(){
	double Mjd_TT=4.974611085861109e+04,dpsi,deps;
	double Rd1=6.230692866336197e-05,Rd2=-3.511084590138603e-05;
	NutAngles(&dpsi,&deps,Mjd_TT);
	_assert(compareDouble(Rd1,dpsi)&&compareDouble(Rd2,deps));
	return 0;
}


int NutMatrix_test(){
	double Mjd_TT=4.974611085861109e+04, **NM;
	double Res[3][3]={{0.999999998058923,-5.716521238294522e-05,-2.478491169341419e-05},
		{5.716608256692436e-05,0.999999997749659,3.511013746598790e-05},
		{2.478290455917464e-05,-3.511155425411916e-05,0.999999999076493}};
	NM=NutMatrix(Mjd_TT);
	_assert(compare(NM,3,3,static_dinamic(3,3,Res),3,3));
	return 0;
}

int PoleMatrix_test(){
	double xp=-5.593861831521886e-07,yp=2.335544384403725e-06, **PM;
	double Res[3][3]={{0.999999999999844,-1.306471258772872e-12,-5.593861831506338e-07},
	{0,0.999999999997273,-2.335544384401602e-06},
	{5.593861831521594e-07,2.335544384401236e-06,0.999999999997116}};
	PM=PoleMatrix(xp,yp);
	_assert(compare(PM,3,3,static_dinamic(3,3,Res),3,3));
	return 0;
}

int GHAMatrix_test(){
	double Mjd_UT1=4.974611015423337e+04, **GM;
	double Res[3][3]={{-0.976451404712871,0.215737465995733,0},
	{-0.215737465995733,-0.976451404712871,0},
	{0,0,1}};
	GM=GHAMatrix(Mjd_UT1);
	_assert(compare(GM,3,3,static_dinamic(3,3,Res),3,3));
	return 0;
}

int gast_test(){
	double Mjd_UT1=4.974611015423337e+04, ga;
	double Res=2.924145635580052;
	ga=gast(Mjd_UT1);
	_assert(compareDouble(ga,Res));
	return 0;
}

int gmst_test(){
	double Mjd_UT1=4.974611015423337e+04, gm;
	double Res=2.924088470683254;
	gm=gmst(Mjd_UT1);
	_assert(compareDouble(gm,Res));
	return 0;
}

int Frac_test(){
	double x=-4.533479501229852, Fr;
	double Res=0.466520498770148;
	Fr=Frac(x);
	_assert(compareDouble(Fr,Res));
	return 0;
}

int EqnEquinox_test(){
	double Mjd_TT=4.974611015423337e+04, eq;
	double Res=5.716489679866389e-05;
	eq=EqnEquinox(Mjd_TT);
	_assert(compareDouble(eq,Res));
	return 0;
}

int Mjday_TDB_test(){
	double Mjd_TT=4.974611015423337e+04, eq;
	double Res=4.974611199287850e+04;
	eq=Mjday_TDB(Mjd_TT);
	_assert(compareDouble(eq,Res));
	return 0;
}

int Cheb3D_test(){
	double t=4.974611199287850e+04,Ta=49744,Tb=49760,
	Cx[]={-103506598.421090,-14927175.9535129,512386.261705712,11989.7464659321,-231.321791175806,-3.01145680139332,0.0557108554642511,0.000160763608359754,4.59304599038190e-05,2.65245498743318e-06,-2.69798189927076e-06,5.19937618744697e-07,5.42165920618080e-08},
	Cy[]={96607067.3707606,-13337710.4674864,-474105.009264048,11293.5192596881,193.380285507684,-3.43928289232978,-0.0344227898113099,0.00107593031980365,4.63124437216978e-06,-1.24347116253730e-05,2.21298188012278e-06,1.39759592836155e-07,-1.68301628426181e-07},
	Cz[]={41887005.1061263,-5782554.14940129,-205553.618975101,4896.48746963361,83.8420287547238,-1.49161452214215,-0.0148800126199174,0.000502341000673383,-2.65374413641448e-06,-6.64088213358266e-06,1.33131703155108e-06,4.77273410670398e-08,-9.29357771479359e-08};
	double Res[]={-91984720.2018902900000000000,105892446.5513364800000000000,45912667.5991248640000000000};
	double *c;
	c=Cheb3D(t,13,Ta,Tb,Cx,Cy,Cz);
	_assert(compareV(Res,3,c,3));
	return 0;
}


int AccelPointMass_test(){
	double GM=1.327124400419394e+20,
	r[]={6.221397628578691e+06,2.867713779657407e+06,3.006155985099499e+06},
	s[]={9.229825172847661e+10,-1.053751960790543e+11,-4.568636722635329e+10};
	double Res[]={-1.868550593441697e-07,-2.003329958831816e-07,-1.599931207554889e-07};
	double *c;
	c=AccelPointMass(r,3,s,3,GM);
	_assert(compareV(Res,3,c,3));
	return 0;
}

int AzElPa_test(){
	double Az,El,*dAds=vector(3),*dEds=vector(3),
	s[]={2.159055448102134e+06,1.212982411020832e+06,7.296691300802083e+05};
	double Res1[]={1.977845622103959e-07,-3.520478390378875e-07,0},
	Res2[]={-9.544240402071819e-08,-5.362065038414828e-08,3.715469614763129e-07};
	AzElPa(s,3,&Az,&El,&dAds,3,&dEds,3);
	_assert(compareV(Res1,3,dAds,3)&&compareV(Res2,3,dEds,3)&&compareDouble(Az,1.058929953815166)&&compareDouble(El,0.286534142298292));
	return 0;
}

int LTC_test(){
	double lon=-2.762343079106937,lat=0.376551295459273, **LM;
	double Res[3][3]={{0.370223471399199,-0.928942722252092,0},
	{0.341586711932422,0.136136938528208,0.929938305587722},
	{-0.863859421119156,-0.344284987681776,0.367715580035218}};
	LM=LTC(lon,lat);
	_assert(compare(LM,3,3,static_dinamic(3,3,Res),3,3));
	return 0;
}

int Legendre_test(){
	double fi=0.413090752583063,**pnm,**dpnm;
	int n=4,m=4;
	double Res1[5][5]={{1,0,0,0,0},
	{0.695317963011514,1.58635838646679,0,0,0},
	{-0.577501372399094,1.4239971919379,1.6244150215478,0,0},
	{-1.16525861472322,-0.288205823728229,1.72531675052923,1.60698222588979,0},
	{-0.347130083762091,-1.63233679659319,0.180194808831757,1.93533053941369,1.56109026693462}};
	double Res2[5][5]={{0,0,0,0,0},
	{1.58635838646679,-0.695317963011514,0,0,0},
	{2.46643548627185,2.62467673998379,-1.4239971919379,0,0},
	{-0.705957209032674,5.5822543329107,2.42383665852906,-2.11307284173668,0},
	{-5.16190218573745,1.47997262283189,7.08338115847732,1.87060140882575,-2.73697069651368}};
	Legendre(n,m,fi,&pnm,&dpnm);
	_assert(compare(static_dinamic(5,5,Res1),5,5,pnm,5,5)&&compare(static_dinamic(5,5,Res2),5,5,dpnm,5,5));
	return 0;
}


int AccelHarmonic_test(){
	double E[3][3]={{-0.978185453896254,0.207733066362260,-4.369502395693629e-04},
	{-0.207733028352522,-0.978185550768511,-1.311456972670817e-04},
	{-4.546617085850979e-04,-3.751581690262889e-05,0.999999895937642}};
	double r[]={6.221397628578691e+06,2.867713779657407e+06,3.006155985099499e+06};
	int n_max=20,m_max=20;
	double Res[]={-5.924148565225369,-2.730766792968872,-2.869335447806865},*a;
	a=AccelHarmonic(r,static_dinamic(3,3,E),3,n_max,m_max);
	_assert(compareV(Res,3,a,3));
	return 0;
}

int G_AccelHarmonic_test(){
	double U[3][3]={{-0.976675972331716,0.214718082511189,-4.360190546746448e-04},
	{-0.214718043811152,-0.976676068937815,-1.342612715042163e-04},
	{-4.546776990745137e-04,-3.750859940872001e-05,0.999999895930642}};
	double r[]={5.542555937228689e+06,3.213514867349193e+06,3.990892975876738e+06};
	int n_max=20,m_max=20;
	double Res[3][3]={{5.700320349077970e-07,0,0},
	{8.676515923511374e-07,0,0},
	{1.081693539184414e-06,0,0}};
	double **a;
	a=G_AccelHarmonic(r,static_dinamic(3,3,U),3,n_max,m_max);
	_assert(compare(static_dinamic(3,3,Res),3,3,a,3,3));
	return 0;
}

int MeasUpdate_test(){
	double x[]={5.738566577691861e+06,3.123975340929588e+06,3.727114481560546e+06,5.199633291810681e+03,-2.474438810446426e+03,-7.195167525538008e+03};
	double z=1.055908489493301,g=1.058929953815166,s=3.909537524467298e-04;
	double G[]={9.591237486030084e-08,2.160503452275371e-07,-3.273827709207121e-07,0,0,0};
	double P[6][6]={{101453348.207834,	120429.109556826,	148186.144851300,	39372.9209797587,	3284.21675106589,	4014.15727751921},
		{120429.109556826	,101309543.076907	,84141.6477758924,	3284.34773933075	,35369.9224513584	,2255.66799781441},
		{148186.144851300	,84141.6477758924	,101344434.103469	,4014.41933186659	,2255.72532205054	,36274.7873542153},
		{39372.9209797587	,3284.34773933075,	4014.41933186659,	1001.21615369228	,1.32096249175600,	1.60455481042780},
		{3284.21675106589	,35369.9224513584	,2255.72532205054	,1.32096249175600	,999.576829598137	,0.892927375360559},
		{4014.15727751921	,2255.66799781441,	36274.7873542153	,1.60455481042780	,0.892927375360559	,999.924178045366}};
	double ResK[]={5.826912064682504e+05,1.312775308420699e+06,-1.989454899791931e+06,1.903673075020185e+02,4.332426595228045e+02,-6.604337991433018e+02};
	double **a;
	_assert(1==1);
	//Acabar luego
	return 0;
}

int Mjday_test(){
	double yr=1995,mon=1,day=29,hr=2,min=38,sec=37;
	double res=4.974611015046295e+04;
	_assert(compareDouble(res,Mjday(yr,mon,day,hr,min,sec)));
	return 0;
}

int position_test(){
	double lon=-2.762343079106937,lat=0.376551295459273,h=3.002000000000000e+02;
	double res[]={-5.512567840036068e+06,-2.196994446669333e+06,2.330804966146887e+06};
	double *a=position(lon,lat,h);
	_assert(compareV(res,3,a,3));
	return 0;
}

int TimeUpdate_test(){
	double Qdt=0.0;
	double P[6][6]={{100000000	,0,	0,	0,	0	,0},
		{0	,100000000,	0,	0,	0,	0},
		{0	,0	,100000000	,0	,0	,0},
		{0	,0	,0,	1000	,0,	0},
		{0	,0,	0,	0	,1000	,0},
		{0,	0,	0,	0,	0	,1000}};
	double Phi[6][6]={{1.00041922218367,	0.000599210758843754,	0.000737344012173523,	37.0053480735614,	0.00742082744985318,	0.00907132659757866},
		{0.000599205935189867,	0.999703770692677	,0.000418689926148395	,0.00742079763278345,	36.9963058593961	,0.00509713273998443},
		{0.000737334362385235,	0.000418687815717319,	0.999877413279163	,0.00907126696371580,	0.00509711970560663	,36.9983505101528},
		{2.34511597329531e-05,	3.25246538718763e-05	,3.97560066025033e-05	,1.00044810751498,	0.000604066118447109	,0.000733457411125333},
		{3.25240005603744e-05	,-1.61854966340218e-05	,2.23408858389155e-05	,0.000604061271608192	,0.999697158529404,	0.000407832040042641},
	{3.97546995847966e-05	,2.23405999335799e-05,	-7.22167527414373e-06,	0.000733447714890626,	0.000407829919202528,	0.999855141838431}};
	double Res[6][6]={{101453348.207834	,120429.109556826	,148186.144851300,	39372.9209797587	,3284.21675106589,	4014.15727751921},
		{120429.109556826,	101309543.076907	,84141.6477758924,	3284.34773933075	,35369.9224513584	,2255.66799781441},
		{148186.144851300,	84141.6477758924	,101344434.103469,	4014.41933186659,	2255.72532205054,	36274.7873542153},
		{39372.9209797587	,3284.34773933075	,4014.41933186659,	1001.21615369228	,1.32096249175600,	1.60455481042780},
		{3284.21675106589,	35369.9224513584	,2255.72532205054	,1.32096249175600,	999.576829598137	,0.892927375360559},
		{4014.15727751921,	2255.66799781441	,36274.7873542153,	1.60455481042780,	0.892927375360559,	999.924178045366}};
	double **Pm=static_dinamic(6,6,P);
	TimeUpdate(&Pm,static_dinamic(6,6,Phi),Qdt,6);
	_assert(compare(static_dinamic(6,6,Res),6,6,Pm,6,6));
	return 0;
}

int VarEqn_test(){
	double x=0;
	double yPhi[]={5542555.93722869,3213514.86734919,3990892.97587674,5394.06842166295,-2365.21337882319,-7061.84554200205,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1};
	double *yPhip;
	double res[]={5394.06842166295,-2365.21337882319,-7061.84554200205,-5.13483678540858,-2.97717622353621,-3.70591776714193,0,0,0,5.70032034907797e-07,8.67651592351137e-07,1.08169353918441e-06,0,0,0,8.67651590574781e-07,-4.23359107770693e-07,6.27183701418232e-07,0,0,0,1.08169354007259e-06,6.27183704970946e-07,-1.46672925360747e-07,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0};
	VarEqn(x,yPhi,&yPhip);
	_assert(compareV(yPhip,42,yPhip,42));
	return 0;
}

int sign_test(){
	_assert(compareDouble(sign_(1,-1.349999919533730e+02),-1));
	return 0;
}


int JPL_Eph_DE430_test(){
	double Mjd_TDB=4.974611199287850e+04,*r_Mercury,*r_Venus,*r_Earth,*r_Mars,*r_Jupiter,*r_Saturn,*r_Uranus,*r_Neptune,*r_Pluto,*r_Moon,*r_Sun;
	double res_Mercury[]={8.377549589569574e+10,-6.529112491344618e+10,-2.339131210124083e+10},
	res_Venus[]={-1.522966557395328e+10,-1.101349926375630e+11,-4.102180362562656e+10},
	res_Earth[]={-9.247096122919232e+10,1.063949183894933e+11,4.613013990940373e+10},
	res_Mars[]={-8.827841300824323e+10,4.696476977829840e+10,2.907102650267133e+10},
	res_Jupiter[]={-2.983859364660936e+11,-7.544982589107291e+11,-3.144105185682280e+11},
	res_Saturn[]={1.482033999505117e+12,-4.538728942363973e+11,-2.494022478112109e+11},
	res_Uranus[]={1.412367984017928e+12,-2.511355045786904e+12,-1.118108651902601e+12},
	res_Neptune[]={1.871250770052329e+12,-3.928976313605398e+12,-1.655020476718540e+12},
	res_Pluto[]={-2.171414794259537e+12,-3.915433128334976e+12,-5.527162503558455e+11},
	res_Moon[]={8.938337231271918e+07,-3.366038321179463e+08,-1.146487877519952e+08},
	res_Sun[]={9.229825172847661e+10,-1.053751960790543e+11,-4.568636722635329e+10};
	
	JPL_Eph_DE430(&r_Mercury,&r_Venus,&r_Earth,&r_Mars,&r_Jupiter,&r_Saturn,&r_Uranus,&r_Neptune,&r_Pluto,&r_Moon,&r_Sun,Mjd_TDB);
	_assert(compareV(r_Mercury,3,res_Mercury,3)&&compareV(r_Venus,3,res_Venus,3)&&compareV(r_Earth,3,res_Earth,3)&&compareV(r_Mars,3,res_Mars,3)&&compareV(r_Jupiter,3,res_Jupiter,3)&&compareV(r_Saturn,3,res_Saturn,3)&&compareV(r_Uranus,3,res_Uranus,3)&&compareV(r_Neptune,3,res_Neptune,3)&&compareV(r_Pluto,3,res_Pluto,3)&&compareV(r_Moon,3,res_Moon,3)||compareV(r_Sun,3,r_Sun,3));
	return 0;
}


int accel_test(){
	double *dY;
	double x=0,y[]={6.221397628578691e+06,2.867713779657407e+06,3.006155985099499e+06,4.645047251617496e+03,-2.752215915881822e+03,-7.507999409869392e+03};
	double v[]={4.645047251617496e+03,-2.752215915881822e+03,-7.507999409869392e+03,-5.924149513140018,-2.730766697881137,-2.869335705562577};
	Accel(x,y,&dY);
	_assert(compareV(dY,6,v,6));
	return 0;
}

int all_tests()
{
	_verify(norma_test);
	_verify(dot_test);
	_verify(vector_test);
	_verify(sumV_test);
	_verify(esc_x_vec_test);
	_verify(static_dinamic_test);
	_verify(mat_x_vec_test);
	_verify(find1_test);
	_verify(find2_test);
	_verify(trasp_test);
	_verify(inv_test);
	_verify(prod_test);
	_verify(eye_test);
	_verify(sum_test);
	_verify(array_test);
	_verify(IERS_test);
	_verify(timediff_test);
	_verify(PrecMatrix_test);
	_verify(R_x_test);
	_verify(R_y_test);
	_verify(R_z_test);
	_verify(MeanObliquity_test);
	_verify(NutAngles_test);
	_verify(NutMatrix_test);
	_verify(PoleMatrix_test);
	_verify(GHAMatrix_test);
	_verify(gast_test);
	_verify(gmst_test);
	_verify(Frac_test);
	_verify(EqnEquinox_test);
	_verify(Mjday_TDB_test);
	_verify(Cheb3D_test);
	_verify(AccelPointMass_test);
	_verify(AzElPa_test);
	_verify(LTC_test);
	_verify(G_AccelHarmonic_test);
	_verify(Legendre_test);
	_verify(AccelHarmonic_test);
	_verify(MeasUpdate_test);
	_verify(Mjday_test);
	_verify(position_test);
	_verify(TimeUpdate_test);
	_verify(VarEqn_test);
	_verify(sign_test);
	_verify(JPL_Eph_DE430_test);



    return 0;
}


int main()
{
	extern double **PC, **Cnm, **Snm, **eopdata;
	int f,c,n,m;
	double **obs,aux1,aux2;
	FILE *fp;
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
	
	for(n=0;n<=180;++n){
		for(m=0;m<=n;++m){
			fscanf(fp,"%d%d%lf%lf%lf%lf",&f,&c, &Cnm[n][m],&Snm[n][m],&aux1,&aux2);
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
    int result = all_tests();

    if (result == 0){
        printf("TEST PASADOS\n");
	}

    printf("Tests realizados: %d\n", tests_run);

    return result != 0;
}
