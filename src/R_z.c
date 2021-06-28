#include "../include/R_z.h"
#include "../include/arrays.h"
#include <math.h>

//--------------------------------------------------------------------------
//  input:
//    angle       - angle of rotation [rad]
//
//  output:
//    rotmat      - vector result
//--------------------------------------------------------------------------

double **R_z(double angle){
	double C,S,**rotmat;
	C = cos(angle);
	S = sin(angle);
	rotmat = array(3,3);

	rotmat[0][0] = C;  
	rotmat[0][1] = S;  
	rotmat[0][2] = 0.0;
	rotmat[1][0] = -1.0*S;  
	rotmat[1][1] = C;  
	rotmat[1][2] = 0.0;
	rotmat[2][0] = 0.0;  
	rotmat[2][1] = 0.0;  
	rotmat[2][2] = 1.0;
	return rotmat;
}