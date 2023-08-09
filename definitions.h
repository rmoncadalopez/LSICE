/*
 * definitions.h
 *
 *  Created on: July 14, 2014
 *      Author: Reid
 */

#ifndef DEFINITIONS_H_
#define DEFINITIONS_H_

// eigen doesn't support vectorization for window$
// #ifdef _WIN32
// 	#define EIGEN_DONT_ALIGN_STATICALLY
// #endif

// C libraries
#include <float.h>		// DBL_MAX
#include <math.h> 		// sin, cos
#include <stdio.h>      // printf, fgets
#include <stdlib.h>     // atoi
#include <omp.h>		// OpenMP
#include <time.h>		// timex
#include <mpi.h>		// MPI

#include <sys/stat.h>

// C++ libraries
#include <string>		// strings
using std::string;

#include <random>

#include <iostream>		// terminal output
#include <cmath>
using std::cout;
using std::endl;

using std::sqrt;
using std::sin;
using std::cos;
//using std::isnan;

#include <iterator> //distance
using std::distance;

#include <list>
using std::list;

#include <fstream>		// file io
using std::ifstream;
using std::getline;
using std::istringstream;

using std::stringstream;

#include <algorithm>	// math stuff
using std::min;
using std::max;
using std::find;
using std::unique;

#include <vector>  // standard vector
//#include <array>
using std::vector;
//using std::array;

//#include <experimental/filesystem>  //copy files from input to output
//namespace fs = std::experimental::filesystem;

#include <ctime>
using std::clock_t;

#include<tr1/array>
using std::tr1::array;

#include <iomanip>

#include <Eigen/Core>    //MODIFICATION due to address in Rml Mac!!!!!!!!!!!!!!
//#include <eigen3/Eigen/Core>
using Eigen::Vector2d;
using Eigen::Vector2i;
using Eigen::Vector3i;
using Eigen::Vector4i;
using Eigen::Vector3d;
using Eigen::Vector4d;
using Eigen::Matrix2d;
using Eigen::Matrix;

using Eigen::ArrayXd;   //For LSModification
using Eigen::ArrayXi;

using Eigen::MatrixXd;
using Eigen::VectorXd;

#include <Eigen/Sparse>  //For Manipulation of Sparse Finite Difference
#include <Eigen/IterativeLinearSolvers> 	//For Classic iterative CG	
//#include <eigen3/Eigen/Sparse>  //For Manipulation of Sparse Finite Difference
//#include <eigen3/Eigen/IterativeLinearSolvers> 	//For Classic iterative CG

typedef Matrix<double, 8, 1> Vector8d;
typedef Matrix<double, 6, 1> Vector6d;

typedef vector<vector<double> > Matrixx;  //Object to generate Matrixx[i][j]
typedef vector<double> Row;


#include <Eigen/StdVector>
//#include <eigen3/Eigen/StdVector>

//#include <eigen3/Eigen/Eigen/Eigenvalues>  //For LSModification
//using Eigen::EigenSolver;


// smooth heaviside function of the negative of the input
vector<double> sHeaviOpp(const vector<double> & vec, double eps = 1.5) {
	vector<double> h(vec.size());
	for (size_t i = 0; i < vec.size(); i++) {
		if (-vec[i] > eps) {
			h[i] = 1;
		}
		else if (-vec[i] < -eps) {
			h[i] = 0;
		}
		else {
			h[i] = 0.5*(1. + -vec[i]/eps + sin(M_PI*-vec[i]/eps)/M_PI);
		}
	}
	return h;
}


// LSMLIB stuff
#include "LSMLIB_config.h"
#include "FMM_Core.h"
#include "lsm_fast_marching_method.h"


#include <bits/stdc++.h> 


//Constants for Controlling Simulation

//Function to define simulation parameters from text file and/or commandline

// size_t start_test = 100;
// double therm_coeff = 0.005;
// double melt_mass = 0.06;
// size_t break_step = 500;
// int break_prob = 99;

// size_t melt_step, break_step;
// double therm_coeff, melt_mass;
// int break_prob;

// //Import text
// char tempfname[300];
// sprintf(tempfname, ("./Input/SeaIceOcean/SimInput.dat").c_str());  //DIR
// std::string file_inputs = tempfname;


// //Choose Case
// size_t caseNo = 4;

// //Run Input Function
// readSimInputFile(file_inputs, caseNo, melt_step, therm_coeff, melt_mass, break_step, break_prob);


// #define START_TEMP melt_step //Can also be used as Step if melt Step is Constant
// #define THERM_COEFF therm_coeff //0.200 //0.005 for 100 steps //High speed makes unstable grains
// #define MELT_MASS melt_mass //0.02 eg for mass //Divided Break and Melt //10000 for point size  //SHOULD IT BE LARGER FOR FASTER SLOPE????
// #define BREAK_STEP break_step 
// #define BREAK_PROB break_prob // 10 for 100 steps


#endif /* DEFINITIONS_H_ */
