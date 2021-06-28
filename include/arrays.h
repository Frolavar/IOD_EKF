#ifndef ARRAY_h_
#define ARRAY_h_

// vectores

double norma(double *v, int n);
double dot(double *v1, int n1, double *v2, int n2);
double *vector(int n);
void freeVector(double *v, int n);
void printVector(double *v, int n);
double *sumV(double *v1, int n1, double *v2, int n2);
double *esc_x_vec(double k,double *v, int n);

// matrices

double **static_dinamic(int nf, int nc,double m[nf][nc]);
double *mat_x_vec(double **m,int nf,int nc,double *v,int n);
int find2(double *v1,int n1,double JD);
int find1(double *v1,int n1,double *v2,int n2,double JD);
double **trasp(double **m, int n);
double **inv(double **m, int n);
double **prod(double **mat1, int nf1, int nc1,double **mat2, int nf2, int nc2);
double **eye(int n);
double **sum(double **mat1, int nf1, int nc1,double **mat2, int nf2, int nc2);
int compareV(double *v1, int n1, double *v2, int n2);
int compare(double **mat1, int nf1, int nc1,double **mat2, int nf2, int nc2);
double **array(int nf, int nc);
void freeArray(double **mat, int nf, int nc);
void printArray(double **mat, int nf, int nc);

//Comparar double
int compareDouble(double d1,double d2);

#endif
