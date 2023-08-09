/*
 *  levelsetElementaryFunctions.h
 *
 *  Created on: September 18, 2014
 *      Author: Reid
 */
#ifndef LEVELSETELEMENTARYFUNCTIONS_H_
#define LEVELSETELEMENTARYFUNCTIONS_H_

#include "definitions.h"

// contains elementary methods that operate on level sets such as gradient, divergence, heaviside, etc


// gaussian filter with sigma = 1;
ArrayXd gaussFilter(const ArrayXd & vec, const int & xd, const int & yd, const int & zd) {
	ArrayXd outg = vec;
	double corn = 0.020586282804687;
	double edge = 0.033941042344737;
	double face = 0.055959318463501;
	double cent = 0.092261318644656;
	for (int z = 1; z < zd-1; z++) {
		for (int y = 1; y < yd-1; y++) {
			for (int x = 1; x < xd-1; x++) {
													// center
				outg(z*yd*xd + y*xd + x) = cent*vec(z*yd*xd + y*xd + x) + 
													// 6 faces
													face*(vec(z*yd*xd + y*xd + x+1)+vec(z*yd*xd + y*xd + x-1)+
															vec(z*yd*xd + (y+1)*xd + x)+vec(z*yd*xd + (y-1)*xd + x)+
															vec((z+1)*yd*xd + y*xd + x)+vec((z-1)*yd*xd + y*xd + x)) +
													// 12 edges: z+1 edges
													edge*(vec((z+1)*yd*xd + y*xd + x+1)+vec((z+1)*yd*xd + y*xd + x-1)+
															vec((z+1)*yd*xd + (y+1)*xd + x)+vec((z+1)*yd*xd + (y-1)*xd + x)+
															// z-1 edges
															vec((z-1)*yd*xd + y*xd + x+1)+vec((z-1)*yd*xd + y*xd + x-1)+
															vec((z-1)*yd*xd + (y+1)*xd + x)+vec((z-1)*yd*xd + (y-1)*xd + x)+
															// z edges
															vec(z*yd*xd + (y-1)*xd + x+1)+vec(z*yd*xd + (y-1)*xd + x-1)+
															vec(z*yd*xd + (y+1)*xd + x+1)+vec(z*yd*xd + (y+1)*xd + x-1)) +
													// 8 corners: z+1 corners
													corn*(vec((z+1)*yd*xd + (y+1)*xd + x+1)+vec((z+1)*yd*xd + (y+1)*xd + x-1)+
															vec((z+1)*yd*xd + (y-1)*xd + x+1)+vec((z+1)*yd*xd + (y-1)*xd + x-1)+
															// z-1 corners
															vec((z-1)*yd*xd + (y+1)*xd + x+1)+vec((z-1)*yd*xd + (y+1)*xd + x-1)+
															vec((z-1)*yd*xd + (y-1)*xd + x+1)+vec((z-1)*yd*xd + (y-1)*xd + x-1));
			}
		}
	}
	return outg;
}



ArrayXd gradx(const ArrayXd & vec, const int & xd, const int & yd, const int & zd) {
	ArrayXd outGx(vec.size());
	for (int z = 0; z < zd; z++) {
		for (int y = 0; y < yd; y++) {
			// first order for x = 0 face
			outGx(z*yd*xd + y*xd) = vec(z*yd*xd + y*xd + 1) - vec(z*yd*xd + y*xd);
			// second order for inside
			for (int x = 1; x < xd-1; x++) {
				outGx(z*yd*xd + y*xd + x) = (vec(z*yd*xd + y*xd + x + 1) - vec(z*yd*xd + y*xd + x - 1 ))/2.0;
			}
			// first order for x = xd face
			outGx(z*yd*xd + y*xd + xd-1) = vec(z*yd*xd + y*xd + xd-1) - vec(z*yd*xd + y*xd + xd-2);
		}
	}
	return outGx;
}


ArrayXd gradz(const ArrayXd & vec, const int & xd, const int & yd, const int & zd) {
	// get outGz
	ArrayXd outGz(vec.size());
	for (int y = 0; y < yd; y++) {
		for (int x = 0; x < xd; x++) {
			// first order for z = 0 face
			outGz(y*xd + x) = vec(yd*xd + y*xd + x) - vec(y*xd + x);
			// second order for inside
			for (int z = 1; z < zd-1; z++) {
				outGz(z*yd*xd + y*xd + x) = (vec((z+1)*yd*xd + y*xd + x) - vec((z-1)*yd*xd + y*xd + x ))/2.0;
			}
			// first order for z = zd face
			outGz((zd-1)*yd*xd + y*xd + x) = vec((zd-1)*yd*xd + y*xd + x) - vec((zd-2)*yd*xd + y*xd + x);
		}
	}
	return outGz;
}


ArrayXd grady(const ArrayXd & vec, const int & xd, const int & yd, const int & zd) {
	ArrayXd outGy(vec.size());
	for (int z = 0; z < zd; z++) {
		for (int x = 0; x < xd; x++) {
			// first order for y = 0 face
			outGy(z*yd*xd + x) = vec(z*yd*xd + xd + x) - vec(z*yd*xd + x);
			// second order for inside
			for (int y = 1; y < yd-1; y++) {
				outGy(z*yd*xd + y*xd + x) = (vec(z*yd*xd + (y+1)*xd + x) - vec(z*yd*xd + (y-1)*xd + x ))/2.0;
			}
			// first order for y = yd face
			outGy(z*yd*xd + (yd-1)*xd + x) = vec(z*yd*xd + (yd-1)*xd + x) - vec(z*yd*xd + (yd-2)*xd + x);
		}
	}
	return outGy;
}


// finds the gradient given a level set vec, dimensions in xd, yd, and zd, and outputs result into outGx, outGy, and outGz
void gradient(const ArrayXd & vec, const int & xd, const int & yd, const int & zd, ArrayXd & outGx, ArrayXd & outGy, ArrayXd & outGz) {
	outGx = gradx(vec, xd, yd, zd);
	outGy = grady(vec, xd, yd, zd);
	outGz = gradz(vec, xd, yd, zd);
}

// finds magnitude of gradient
ArrayXd gradmag(const ArrayXd & vec, const int & xd, const int & yd, const int & zd) {
	ArrayXd outgmag(vec.size());
	ArrayXd gx, gy, gz;
	gradient(vec, xd, yd, zd, gx, gy, gz);
	outgmag = sqrt(gx*gx + gy*gy + gz*gz);
	return outgmag;
}

ArrayXd divergence(const ArrayXd & v1, const ArrayXd & v2, const ArrayXd & v3, const int & xd, const int & yd, const int & zd) {
	return gradx(v1,xd,yd,zd) + grady(v2,xd,yd,zd) + gradz(v3,xd,yd,zd);
}

// smooth heaviside function
ArrayXd sHeavi(const ArrayXd & vec, double eps = 1.5) {
	ArrayXd h(vec.size());
	for (int i = 0; i < vec.size(); i++) {
		if (vec(i) > eps) {
			h(i) = 1;
		}
		else if (vec(i) < -eps) {
			h(i) = 0;
		}
		else {
			h(i) = 0.5*(1. + vec(i)/eps + sin(M_PI*vec(i)/eps)/M_PI);
		}
	}
	return h;
}

// smooth dirac function
ArrayXd sDirac(const ArrayXd & vec, double eps = 1.5) {
	ArrayXd d(vec.size());
	for (int i = 0; i < vec.size(); i++) {
		if (fabs(vec(i)) > eps) {
			d(i) = 0;
		}
		else {
			d(i) = (1.+cos(M_PI*vec(i)/eps))/(2.*eps);
		}
	}
	return d;
}



ArrayXd dp(const ArrayXd & vec) {
	ArrayXd dp(vec.size());
	for (int i = 0; i < vec.size(); i++) {
		if (vec(i) >= 1) {
			dp(i) = (vec(i)-1)/vec(i);
		}
		else if (vec(i) == 0) {
			dp(i) = 1;
		}
		else {
			dp(i) = sin(2.*M_PI*vec(i))/(2.*M_PI)/vec(i);
		}
	}
	return dp;
}


ArrayXd laplacian(const ArrayXd & vec, const int & xd, const int & yd, const int & zd) {
	// initialize and zero the edge values
	ArrayXd lap = ArrayXd::Zero((vec.size()));
	// find laplacian for inner grid
	for (int z = 1; z < zd-1; z++) {
		for (int y = 1; y < yd-1; y++) {
			for (int x = 1; x < xd-1; x++) {
				lap(z*yd*xd + y*xd + x) = 	vec(z*yd*xd + y*xd + x+1) + vec(z*yd*xd + y*xd + x-1)  // x
												+	vec(z*yd*xd + (y+1)*xd + x) + vec(z*yd*xd + (y-1)*xd + x)  // y
												+	vec((z+1)*yd*xd + y*xd + x) + vec((z-1)*yd*xd + y*xd + x)  // z
												-	6*vec(z*yd*xd + y*xd + x)	; // center
			}
		}
	}
	return lap;
}

ArrayXd laplacian2d(const ArrayXd & vec, const int & xd, const int & yd) {
    // initialize and zero the edge values
    ArrayXd lap = ArrayXd::Zero((vec.size()));
    // find laplacian for inner grid
    for (int y = 1; y < yd-1; y++) {
        for (int x = 1; x < xd-1; x++) {
            lap(y*xd + x) =     vec(y*xd + x+1) + vec(y*xd + x-1)  // x
                                            +    vec((y+1)*xd + x) + vec((y-1)*xd + x)  // y
                                            -    4*vec(y*xd + x)    ; // center
        }
    }
    return lap;
}
	
// this might leak memory idk
//LSMLIB_REAL * laplacian(const LSMLIB_REAL * const vec, const int & xd, const int & yd, const int & zd) {
//	// initialize and zero the edge values
//	LSMLIB_REAL * lap = (LSMLIB_REAL*) malloc(xd*yd*zd*sizeof(LSMLIB_REAL));
//	// find laplacian for inner grid
//	for (int z = 1; z < zd-1; z++) {
//		for (int y = 1; y < yd-1; y++) {
//			for (int x = 1; x < xd-1; x++) {
//				lap[z*yd*xd + y*xd + x] = 	vec[z*yd*xd + y*xd + x+1] + vec[z*yd*xd + y*xd + x-1]  // x
//												+	vec[z*yd*xd + (y+1)*xd + x] + vec[z*yd*xd + (y-1)*xd + x]  // y
//												+	vec[(z+1)*yd*xd + y*xd + x] + vec[(z-1)*yd*xd + y*xd + x]  // y 
//												-	6*vec[z*yd*xd + y*xd + x]	; // center
//			}
//		}
//	}	
//	return lap;
//}

ArrayXd laplacian(const LSMLIB_REAL * const vec, const int & xd, const int & yd, const int & zd) {
	// initialize and zero the edge values
	ArrayXd lap = ArrayXd::Zero(xd*yd*zd);
	// find laplacian for inner grid
	for (int z = 1; z < zd-1; z++) {
		for (int y = 1; y < yd-1; y++) {
			for (int x = 1; x < xd-1; x++) {
				lap(z*yd*xd + y*xd + x) = 	vec[z*yd*xd + y*xd + x+1] + vec[z*yd*xd + y*xd + x-1]  // x
												+	vec[z*yd*xd + (y+1)*xd + x] + vec[z*yd*xd + (y-1)*xd + x]  // y
												+	vec[(z+1)*yd*xd + y*xd + x] + vec[(z-1)*yd*xd + y*xd + x]  // y 
												-	6*vec[z*yd*xd + y*xd + x]	; // center
			}
		}
	}	
	return lap;
}




#endif

