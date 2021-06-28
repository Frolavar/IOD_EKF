#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include "../include/arrays.h"

// matrices

double *mat_x_vec(double **m,int nf,int nc,double *v,int n){
	double *res;
	int i,j;
	
	res=vector(n);
	for(i=0;i<nf;i++){
		for(j=0;j<nc;j++){
			res[i]=res[i]+m[i][j]*v[j];
		}
	}
	return res;
}
double *esc_x_vec(double k,double *v, int n){
	double *result;
	int i;
	
	result=vector(n);
	for(i=0;i<n;++i){
		result[i]=k*v[i];
	}
	return result;
}

double *sumV(double *v1, int n1, double *v2, int n2){
	double *result;
	int i;
	
	if(n1!=n2){
		printf("sumV: n1 != n2\n");
		exit(EXIT_FAILURE);
	}
	
	result=vector(n1);
	for(i=0;i<n1;++i){
		result[i]=v1[i]+v2[i];
	}
	return result;
}

int find2(double *v,int n,double JD){
	int i=0;
	while(i<n){
		if(fabs(v[i]-JD)<pow(10,-10)){
			break;
		}
		++i;
	}
	return i;
}

int find1(double *v1,int n1,double *v2,int n2,double JD){
	
	int i1=0,i2=0, ind1=0,ind2=0;
	
	if(n1 != n2){
        printf("find1: n1 != n2\n");
        exit(EXIT_FAILURE);
    }
	while(i1<n1){
		if(v1[i1]>JD){
			ind1=i1;
			break;
		}
		++i1;
	}
	
	while(i2<n2){
		if(v2[i2]>JD){
			ind2=i2;
			break;
		}
		++i2;
	}
	return (ind1<ind2 ? ind1:ind2);			
}

double dot(double *v1, int n1, double *v2, int n2)
{
    double sum;
    int i;
    
    if(n1 != n2){
        printf("dot: n1 != n2\n");
        exit(EXIT_FAILURE);
    }

    sum = 0.0;
    for(i = 0; i < n1; i++)
        sum += v1[i]*v2[i];
    
    return(sum);
}

double norma(double *v, int n)
{
    double sum;
    int i;
    
    sum = 0.0;
    for(i = 0; i < n; i++)
        sum += v[i]*v[i];
    
    return (sqrt(sum));
}


double *vector(int n)
{
    double *v;
 
    v = (double *) calloc(n, sizeof(double));
    if (v == NULL){
        printf("vector: memory not allocated\n");
        exit(EXIT_FAILURE);
    }

    return v;
}

void freeVector(double *v, int n)
{
    free(v);
}

void printVector(double *v, int n)
{
    int i;

    for(i = 0; i < n; i++)
        printf("%5.19lf ", v[i]);
     printf("\n");
}


double **trasp(double **m, int n)
{
    double **result;
    int i, j;
    
    result = array(n, n);
    
    for(i = 0; i < n; i++)
        for(j = 0; j < n; j++)
            result[i][j] = m[j][i];

    return result;
}

// Using Gauss Jordan Method
double **inv(double **m, int n)
{
    double **mat, ratio, **result;
    int i, j, k;
    
    mat = array(n, 2*n);
    result = array(n, n);
    
    for(i = 0; i < n; i++)
        for(j = 0; j < n; j++)
            mat[i][j] = m[i][j];
    
    for(i = 0; i < n; i++){
         for(j = 0; j < n; j++){
              if(i == j){
                  mat[i][j+n] = 1;
              } else {
                  mat[i][j+n] = 0;
              }
         }
    }
    
    /* Applying Gauss Jordan Elimination */
    
    for(i = 0; i < n; i++){
         if(mat[i][i] == 0.0){
             printf("inv: error\n");
             exit(EXIT_FAILURE);
         }
         for(j = 0; j < n; j++){
              if(i != j){
                   ratio = mat[j][i]/mat[i][i];
                   for(k = 0; k < 2*n; k++){
                       mat[j][k] = mat[j][k] - ratio*mat[i][k];
                   }
              }
         }
    }
    
    /* Row Operation to Make Principal Diagonal to 1 */
    
    for(i = 0; i < n; i++){
         for(j = n; j < 2*n;j++){
             mat[i][j] = mat[i][j]/mat[i][i];
         }
    }
    
    for(i = 0; i < n; i++)
        for(j = 0; j < n; j++)
            result[i][j] = mat[i][j+n];
 
    freeArray(mat, n, 2*n);
    
    return result;
}

double **prod(double **mat1, int nf1, int nc1,
              double **mat2, int nf2, int nc2)
{
    double **result;
    int i, j, k;
    
    if(nf1 != nc2){
        printf("prod: nf1 != nc2\n");
        exit(EXIT_FAILURE);
    }
    result = array(nf1, nc2);
    for(i = 0; i < nf1; i++){
        for(j = 0; j < nc2; j++){
            result[i][j] = 0;//inicializar
            for(k = 0; k < nc1; k++){
                result[i][j] = result[i][j] + mat1[i][k]*mat2[k][j];
            }
        }
     }
    
    return result;
}

double **eye(int n)
{
    double **result;
    int i;

    result = array(n, n);
    for (i = 0; i < n; ++i)
        result[i][i] = 1.0;
 
    return result;
}

double **sum(double **mat1, int nf1, int nc1,
             double **mat2, int nf2, int nc2)
{
    double **result;
    int i, j;
    
    if((nf1 != nf2) || (nc1 != nc2)){
        printf("sum: different number of rows or columns\n");
        exit(EXIT_FAILURE);
    }
    result = array(nf1, nc2);
    for (i = 0; i < nf1; ++i)
        for (j = 0; j < nc1; ++j) {
            result[i][j] = mat1[i][j] + mat2[i][j];
        }
    
    return result;
}

int compareV(double *v1, int n1, double *v2, int n2)
{
    int i, j;
    
    if(n1 != n2){
        printf("vec: different dimension\n");
        return 0;
    }
    
    for (i = 0; i < n1; ++i)
        if(fabs(v1[i]-v2[i]) > pow(10, -1)){
                return 0;
        }
 
    return 1;
}

int compare(double **mat1, int nf1, int nc1,
            double **mat2, int nf2, int nc2)
{
    int i, j;
    
    if((nf1 != nf2) || (nc1 != nc2)){
        printf("sum: different number of rows or columns\n");
        return 0;
    }
    
    for (i = 0; i < nf1; ++i)
        for (j = 0; j < nc1; ++j) {
            if(fabs(mat1[i][j]-mat2[i][j]) > pow(10, -2))
                return 0;
        }
 
    return 1;
}

double **array(int nf, int nc)
{
    double **m;
    int i;

    m = (double **) calloc(nf, sizeof(double *));
    if (m == NULL){
        printf("array: memory not allocated\n");
        exit(EXIT_FAILURE);
    }
    for(i = 0; i < nf; i++){
        m[i] = (double *) calloc(nc, sizeof(double));
        if (m[i] == NULL){
            printf("array: memory not allocated\n");
            exit(EXIT_FAILURE);
        }
    }
    
    return m;
}

void freeArray(double **mat, int nf, int nc)
{
    int i;

    for(i = 0; i < nf; i++)
        free(mat[i]);
    free(mat);
}

void printArray(double **mat, int nf, int nc)
{
    int i, j;

    for(i = 0; i < nf; i++){
        for(j = 0; j < nc; j++){
            printf("%5.19lf ", mat[i][j]);
        }
    printf("\n");
    }
}


int compareDouble(double d1,double d2)
{
    if(fabs(d1-d2) > pow(10, -2)){
		return 0;
	}
    return 1;
}

double **static_dinamic(int nf, int nc,double m[nf][nc]){
	int i,j;
	double **x=array(nf,nc);
	for(i=0;i<nf;++i){
		for(j=0;j<nc;++j){
			x[i][j]=m[i][j];
		}
	}
	return x;
}
