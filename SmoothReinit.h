/*
 *  smoothReinit.h
 *
 *  Created on: September 24, 2014
 *      Author: Reid
 */
#ifndef SMOOTHREINIT_H
#define SMOOTHREINIT_H

#include "definitions.h"
#include "levelsetElementaryFunctions.h"


// this doesn't work as well as i'd like it but at least it works (don't use it)
void smoothReinit(ArrayXd & lset, const int & xd, const int & yd, const int & zd) {
	
	size_t ngridpts = xd*yd*zd;
//	ArrayXd temp(ngridpts);
	Vector3i dims(xd, yd, zd);
	ArrayXd lap(ngridpts);
	ArrayXd gmag(ngridpts);
	LSMLIB_REAL * phi;
	LSMLIB_REAL * distance_function;
	LSMLIB_REAL * mask;
	
	// set up the variables
	phi = (LSMLIB_REAL*) malloc(ngridpts*sizeof(LSMLIB_REAL));
	distance_function = (LSMLIB_REAL*) malloc(ngridpts*sizeof(LSMLIB_REAL));
	mask = (LSMLIB_REAL*) malloc(ngridpts*sizeof(LSMLIB_REAL));
	int spatial_derivative_order = 2;
	int grid_dims[3];
	LSMLIB_REAL dx[3];
	for (size_t j = 0; j < 3; j++) {
		grid_dims[j] = dims(j);
		dx[j] = 1.;
	}
	for (size_t j = 0; j < ngridpts; j++) {
		phi[j] = lset(j);
		distance_function[j] = 0;
		mask[j] = 1;
	}
	
//	// smooth slightly
//	for (size_t j = 0; j < 6; j++) {
//		lap = laplacian(phi, dims(0), dims(1), dims(2));
//		for (size_t k = 0; k < ngridpts; k++) {
//			phi[k] += .025*lap(k);
//		}
//	}
	// reinit
	computeDistanceFunction3d(distance_function,phi,mask,spatial_derivative_order,grid_dims,dx);
	
	// do another pass
	// make the output distance function the new input uninitialized level set
	free(phi);
	phi = distance_function;
	distance_function = (LSMLIB_REAL*) malloc(ngridpts*sizeof(LSMLIB_REAL));

	// smoothed sharpen and reinit again
//	for (size_t j = 0; j < 5; j++) {
//		lap = laplacian(phi, dims(0), dims(1), dims(2));
//		lap = gaussFilter(lap,dims(0), dims(1), dims(2));
////		for(size_t k = 0; k < (size_t)ngridpts; k++) {
////			gmag(k) = phi[k];
////		}
////		gmag = gradmag(gmag, dims(0), dims(1), dims(2));
//		for (size_t k = 0; k < (size_t)ngridpts; k++) {
//			phi[k] -= .005*double((lap(k) > 0) - (lap(k) < 0));
//		}
//	}
	
	// hack to make the gradient better and the grains slightly larger
	for (size_t j = 0; j < ngridpts; j++) {
		phi[j] += -1;
	}
	computeDistanceFunction3d(distance_function,phi,mask,spatial_derivative_order,grid_dims,dx);
	for (size_t j = 0; j < ngridpts; j++) {
		lset(j) = distance_function[j] + 1;
	}
	
	delete[] phi;  
	delete[] distance_function;
	delete[] mask;

//	return outlset;
}

// use this method
void reinit(ArrayXd & lset, const int & xd, const int & yd, const int & zd) {
	
	size_t ngridpts = xd*yd*zd;
//	ArrayXd temp(ngridpts);
	Vector3i dims(xd, yd, zd);
	ArrayXd lap(ngridpts);
	ArrayXd gmag(ngridpts);
	LSMLIB_REAL * phi;
	LSMLIB_REAL * distance_function;
	LSMLIB_REAL * mask;
	
	
	// set up the variables
	phi = (LSMLIB_REAL*) malloc(ngridpts*sizeof(LSMLIB_REAL));
	distance_function = (LSMLIB_REAL*) malloc(ngridpts*sizeof(LSMLIB_REAL));
	mask = (LSMLIB_REAL*) malloc(ngridpts*sizeof(LSMLIB_REAL));
	int spatial_derivative_order = 1;
	int grid_dims[3];
	LSMLIB_REAL dx[3];
	for (size_t j = 0; j < 3; j++) {
		grid_dims[j] = dims(j);
		dx[j] = 1.;
	}
	for (size_t j = 0; j < ngridpts; j++) {
		phi[j] = lset(j);
		distance_function[j] = 0;
		mask[j] = 1;
	}
	
    //SMOOTHER PARTICLES TRY (NEW)
    // smooth slightly
    for (size_t j = 0; j < 6; j++) {
        lap = laplacian(phi, dims(0), dims(1), dims(2));
        for (size_t k = 0; k < ngridpts; k++) {
            phi[k] += .025*lap(k);  //Default: 0.025
        }
    }
	
	// reinit
	computeDistanceFunction3d(distance_function,phi,mask,spatial_derivative_order,grid_dims,dx);
	
	for (size_t j = 0; j < ngridpts; j++) {
		lset(j) = distance_function[j];
	}
	
	delete[] phi;  
	delete[] distance_function;
	delete[] mask;

//	return outlset;
}


// use this method
void reinit2D(ArrayXd & lset, const int & xd, const int & yd) {
    
    size_t ngridpts = xd*yd;
//    ArrayXd temp(ngridpts);
    Vector2i dims(xd, yd);
    ArrayXd lap(ngridpts);
    ArrayXd phicopy(ngridpts);
    ArrayXd gmag(ngridpts);
    LSMLIB_REAL * phi;
    LSMLIB_REAL * distance_function;
    LSMLIB_REAL * mask;
    
    
    // set up the variables
    phi = (LSMLIB_REAL*) malloc(ngridpts*sizeof(LSMLIB_REAL));
    distance_function = (LSMLIB_REAL*) malloc(ngridpts*sizeof(LSMLIB_REAL));
    mask = (LSMLIB_REAL*) malloc(ngridpts*sizeof(LSMLIB_REAL));
    int spatial_derivative_order = 1;
    int grid_dims[2];
    LSMLIB_REAL dx[2];
    for (size_t j = 0; j < 2; j++) {
        grid_dims[j] = dims(j);
        dx[j] = 1.;
    }
    for (size_t j = 0; j < ngridpts; j++) {
        phi[j] = lset(j);
        phicopy(j) = lset(j);
        distance_function[j] = 0;
        mask[j] = 1;
    }
    
    //SMOOTHER PARTICLES TRY (NEW) ????? TODO
    // smooth slightly
    for (size_t j = 0; j < 4; j++) {
        //lap = laplacian2d(phi, dims(0), dims(1));
        lap = laplacian2d(phicopy, dims(0), dims(1));
        for (size_t k = 0; k < ngridpts; k++) {
            phi[k] += .025*lap(k);  //Default: 0.025
        }
    }
    
    // reinit
    computeDistanceFunction2d(distance_function,phi,mask,spatial_derivative_order,grid_dims,dx);
    
    for (size_t j = 0; j < ngridpts; j++) {
        lset(j) = distance_function[j];
    }
    
    delete[] phi;
    delete[] distance_function;
    delete[] mask;

//    return outlset;
}

void smooth(ArrayXd & lset, const int & xd, const int & yd, const int & zd) {
	
	double shift = 1;
	
	size_t ngridpts = xd*yd*zd;
//	ArrayXd temp(ngridpts);
	Vector3i dims(xd, yd, zd);
	ArrayXd lap(ngridpts);
	ArrayXd gmag(ngridpts);
	LSMLIB_REAL * phi;
	LSMLIB_REAL * distance_function;
	LSMLIB_REAL * mask;
	
	
	// set up the variables
	phi = (LSMLIB_REAL*) malloc(ngridpts*sizeof(LSMLIB_REAL));
	distance_function = (LSMLIB_REAL*) malloc(ngridpts*sizeof(LSMLIB_REAL));
	mask = (LSMLIB_REAL*) malloc(ngridpts*sizeof(LSMLIB_REAL));
	int spatial_derivative_order = 1;
	int grid_dims[3];
	LSMLIB_REAL dx[3];
	for (size_t j = 0; j < 3; j++) {
		grid_dims[j] = dims(j);
		dx[j] = 1.;
	}
	for (size_t j = 0; j < ngridpts; j++) {
		phi[j] = lset(j)-shift;
		distance_function[j] = 0;
		mask[j] = 1;
	}
	
	// reinit
	computeDistanceFunction3d(distance_function,phi,mask,spatial_derivative_order,grid_dims,dx);
	
	for (size_t j = 0; j < ngridpts; j++) {
		lset(j) = distance_function[j]+shift;
	}
	
	delete[] phi;  
	delete[] distance_function;
	delete[] mask;
//	return outlset;
}



#endif
