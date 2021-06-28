#ifndef GLOBALES_h_
#define GLOBALES_h_

typedef struct{
	double Mjd_UTC, Mjd_TT;
	int n,m,sun,moon,planets;
} Param;

double **PC, **Cnm, **Snm, **eopdata;
int n_eqn;
Param AuxParam;

#endif