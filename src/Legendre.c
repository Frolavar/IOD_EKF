#include "../include/Legendre.h" 
#include <math.h>
#include "../include/arrays.h"
#include <stdio.h>

void Legendre(int n, int m, double fi, double*** pnm, double*** dpnm){
	
	//voy a dejar la fila y columna 0 libres
	double** aux_pnm = array(n + 1, m + 1);
	double** aux_dpnm = array(n + 1, m + 1);

	aux_pnm[0][0] = 1;
	aux_dpnm[0][0] = 0;

	aux_pnm[1][1] = sqrt(3) * cos(fi);
	aux_dpnm[1][1] = -sqrt(3) * sin(fi);

	/* Coeficientes diagonales */

	for (int i = 1; i < n; i++)
	{
		aux_pnm[i+1][i+1] = sqrt((2.0 * (i+1) + 1) / (2.0 * (i+1))) * cos(fi) * aux_pnm[i][i];
	}

	for (int i = 1; i < n; i++)
	{
		aux_dpnm[i + 1][i + 1] = sqrt((2.0 * (i+1) + 1) / (2.0 * (i+1))) * ((cos(fi) * aux_dpnm[i][i]) - (sin(fi) * aux_pnm[i][i]));
	}
	/* Primer paso coeficientes horizontales*/
	for (int i = 0; i < n; i++)
	{
		aux_pnm[i + 1] [i] = sqrt(2.0 * (i+1) + 1) * sin(fi) * aux_pnm[i][i];
	}
	for (int i = 0; i < n; i++)
	{
		aux_dpnm[i + 1][i] = sqrt(2.0 * (i+1) + 1) * ((cos(fi) * aux_pnm[i][i]) + (sin(fi) * aux_dpnm[i][i]));
	}
	
	/* Segundo paso coeficientes horizontales*/
	int j = -1, k = 1;
	while(1)
	{
		for (int i = k; i < n; i++)
		{
			aux_pnm[i + 1][j + 1] = sqrt((2.0 * (i+1) + 1) / ((1.0 * (i+1) - (j+1)) * (1.0 * (i+1) + (j+1)))) * ((sqrt(2.0 * (i+1) - 1) * sin(fi) * aux_pnm[i][j + 1])- (sqrt(((1.0 * (i+1) + (j+1) - 1.0) * (1.0 * (i+1) - (j+1) - 1.0)) / (2.0 * (i+1) - 3)) * aux_pnm[i - 1][j + 1]));
		}
		j = j + 1;
		k = k + 1;
		if (j >= m)
		{
			break;
		}	
	}
	j = -1;
	k = 1;
	
	while (1)
	{
		for (int i = k; i < n; i++)
		{
			aux_dpnm[i + 1][j + 1] = sqrt((2.0 * (i+1) + 1) / ((1.0 * (i+1) - 1.0 * (j+1)) * (1.0 * (i+1) + (j+1)))) * ((sqrt(2.0 * (i+1) - 1) * sin(fi) * aux_dpnm[i][j + 1])
				+ (sqrt(2.0 * (i+1) - 1) * cos(fi) * aux_pnm[i][j + 1]) - (sqrt(((1.0 * (i+1) +(j+1) - 1.0) * (1.0 * (i+1) - (j+1) - 1.0)) / (2.0 * (i+1) - 3)) * aux_dpnm[i - 1][j + 1]));
		}
		j = j + 1;
		k = k + 1;
		if (j >= m)
		{
			break;
		}
	}
	


	
	*pnm=array(n+1,m+1);
	*dpnm=array(n+1,m+1);
	
	*pnm=aux_pnm;
	*dpnm=aux_dpnm;
	
	//freeArray(aux_pnm, n + 2, m + 2);
	//freeArray(aux_dpnm, n + 2, m + 2);

}