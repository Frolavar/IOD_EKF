#include "../include/AccelPointMass.h"
#include "../include/arrays.h"
#include <math.h>


double *AccelPointMass(double *r, int nr, double *s, int ns, double GM) { 
	double *d; 
	// Relative position vector of satellite w.r.t. point mass 
	d = sumV(r, nr, esc_x_vec(-1.0, s, ns), ns); 
	// Acceleration 
	return(sumV(esc_x_vec(-GM/pow(norma(d,3),3.0), d, ns), ns, esc_x_vec(-GM/pow(norma(s,3),3.0), s, ns), ns)); 
} 