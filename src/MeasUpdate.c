#include "../include/MeasUpdate.h"
#include "../include/arrays.h"

void MeasUpdate(double ***K, double **x, double *z,double *g, double *s,double **G,double ***P,int n){
	
	double **Inv_W,**Menos;
	int i,j;

	Inv_W = array(n,n);
	for(i=0;i<n;++i){
		Inv_W[i][i] = s[i]*s[i];    // Inverse weight (measurement covariance)
	}

	// Kalman gain
	*K = prod(prod(*P,n,n,inv(G,n),n,n),n,n,inv(sum(Inv_W,n,n,prod(prod(G,n,n,*P,n,n),n,n,inv(G,n),n,n),n,n),n),n,n);

	// State update
	*x = sumV(*x,n,mat_x_vec(*K,n,n,sumV(z,n,esc_x_vec(-1,g,n),n),n),n);

	// Covariance update

	Menos=prod(*K,n,n,G,n,n);

	for(i=0;i<n;++i){
		for(j=0;j<n;++j){
			Menos[i][j] = Menos[i][j]*-1;
		}
	}

	*P=prod(sum(eye(n),n,n,Menos,n,n),n,n,*P,n,n); 
}